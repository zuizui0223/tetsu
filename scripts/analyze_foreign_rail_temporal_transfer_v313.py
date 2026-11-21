"""Retrospective foreign-rail frequency versus 2014->2017 recorded transitions.

Original site/taxon data stay in memory; only study/block/habitat aggregates
are exported. Design was frozen in commit 865cac4 before target effects.
"""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import importlib.metadata
import platform
import json
import sys
import numpy as np
from analyze_incremental_foreign_rail_v304 import (
    logit, metrics, fit_predict, mean_frequencies, crossfit_training_frequencies,
    local_source,
)

ROOT = Path(__file__).resolve().parents[1]
HABITATS = ("Track", "Interior Plot 1", "Interior Plot 2")
HABITAT_WEIGHTS = np.array([1., .5, .5])
MODELS = ("M0", "M1", "M1b", "M2")
PENALTIES = (.01, .001, .1)
COMPARISONS = {
    "foreign_common": ("M0", "M1"),
    "foreign_all_habitat_slopes": ("M1", "M2"),
    "foreign_track_given_interior_difference": ("M1b", "M2"),
}


def prepare_data(species, plots, meta, donors):
    if len(donors) != 4:
        raise ValueError("The fixed foreign score requires all four studies")
    known = [t for t in species if t]
    if len(known) != len(set(known)):
        raise ValueError("Duplicate binomial mapping")
    support = Counter(t for d in donors for t, f in d.items() if f > 0)
    foreign = mean_frequencies(donors)
    baseline = np.zeros(len(species), bool)
    for (_, year, _), values in plots.items():
        values = np.asarray(values, float)
        if len(values) != len(species):
            raise ValueError("Inconsistent source taxon columns")
        if year == 2014:
            baseline |= np.isfinite(values) & (values > 0)
    chosen = [i for i, t in enumerate(species) if t and baseline[i] and support[t] >= 2]
    sites = sorted((s for s in meta if all(
        (s, y, h) in plots and np.isfinite(plots[s, y, h]).all()
        for y in (2014, 2017) for h in HABITATS)), key=lambda s: (meta[s]["lat"], s))
    if len(sites) < 5 or not chosen:
        raise ValueError("Insufficient complete pairs or baseline eligible taxa")
    abundance = np.array([[plots[s, 2014, h][chosen] for h in HABITATS] for s in sites], float)
    later = np.array([[plots[s, 2017, h][chosen] for h in HABITATS] for s in sites], float)
    if np.any(abundance < 0) or np.any(later < 0):
        raise ValueError("Negative abundance is not a nonrecord")
    blocks = np.empty(len(sites), int)
    for b, indices in enumerate(np.array_split(np.arange(len(sites)), 5)):
        blocks[indices] = b
    taxa = [species[i] for i in chosen]
    return dict(
        sites=sites, taxa=taxa, abundance0=abundance, presence1=(later > 0).astype(float),
        foreign=logit([foreign[t] for t in taxa]), blocks=blocks,
        altitude=np.array([meta[s]["altitude"] for s in sites]),
        station=np.array([meta[s]["station"] for s in sites], float),
        audit={
            "source_taxon_columns": len(species), "mapped_binomials": len(known),
            "foreign_two_study_supported_before_baseline_restriction":
                sum(bool(t) and support[t] >= 2 for t in species),
            "baseline_recorded_binomials": sum(bool(t) and bool(baseline[i]) for i,t in enumerate(species)),
            "eligible_baseline_foreign_supported_taxa": len(taxa),
            "baseline_excluded_foreign_supported_taxa": sum(bool(t) and support[t] >= 2 and not baseline[i]
                                                            for i,t in enumerate(species)),
            "source_transects": len(meta), "complete_paired_transects": len(sites),
            "station_transects": sum(bool(meta[s]["station"]) for s in sites),
            "between_transects": sum(not meta[s]["station"] for s in sites),
            "block_sizes": [int(np.sum(blocks == b)) for b in range(5)],
        })


def baseline_frequencies(abundance0, train_sites):
    return crossfit_training_frequencies((np.asarray(abundance0) > 0).astype(float), train_sites)


def make_rows(c, endpoint, train_sites):
    if endpoint not in ("persistence", "new_record"):
        raise ValueError("Unknown endpoint")
    presence = c["abundance0"] > 0
    risk = presence if endpoint == "persistence" else ~presence
    index = np.argwhere(risk)
    site, habitat, taxon = index.T
    local_frequency = logit(baseline_frequencies(c["abundance0"], train_sites))
    neighbors = (presence.sum(axis=1, keepdims=True) - presence) / 2
    onehot = np.eye(3)[habitat]
    columns = [
        (habitat == 1).astype(float), (habitat == 2).astype(float),
        c["altitude"][site], c["station"][site],
    ]
    feature_names = ["interior1", "interior2", "altitude", "station"]
    for name, values in (("baseline_frequency", local_frequency), ("other_habitat_presence", neighbors)):
        for h, label in enumerate(HABITATS):
            columns.append(values[site, habitat, taxon] * onehot[:,h])
            feature_names.append(name + ":" + label)
    if endpoint == "persistence":
        for h, label in enumerate(HABITATS):
            columns.append(np.log1p(c["abundance0"][site, habitat, taxon]) * onehot[:,h])
            feature_names.append("baseline_abundance:" + label)
    nbase = len(columns)
    f = c["foreign"][taxon]
    columns.extend([f, f * np.array([0.,1.,-1.])[habitat], f * (habitat == 0)])
    feature_names.extend(["foreign", "foreign_interior_contrast", "foreign_track"])
    return dict(x=np.column_stack(columns), y=c["presence1"][site,habitat,taxon],
                index=index, base_columns=nbase, feature_names=feature_names)


def habitat_weights(habitat):
    habitat = np.asarray(habitat)
    counts = np.bincount(habitat, minlength=3)
    if np.any(counts == 0):
        raise ValueError("Training risk set missing a habitat")
    weights = HABITAT_WEIGHTS[habitat] / counts[habitat]
    return weights / weights.sum()


def balanced_average(values):
    if any(v is None for v in values):
        return None
    return float(np.dot(HABITAT_WEIGHTS, values) / HABITAT_WEIGHTS.sum())


def foreign_slopes(beta, sd):
    raw = np.asarray(beta)[1:] / np.asarray(sd)
    common, interior, track = raw[-3:]
    return {
        "Track": float(common+track),
        "Interior Plot 1": float(common+interior),
        "Interior Plot 2": float(common-interior),
        "Track_minus_interior_mean": float(track),
        "Track_minus_Interior_Plot_1": float(track-interior),
        "Track_minus_Interior_Plot_2": float(track+interior),
    }


def forecast_fold(c, endpoint, block, penalty):
    train_sites = c["blocks"] != block
    rows = make_rows(c, endpoint, train_sites)
    train = train_sites[rows["index"][:,0]]
    test = ~train
    if len(np.unique(rows["y"][train])) != 2:
        raise ValueError("Training risk set lacks both outcomes")
    w = habitat_weights(rows["index"][train,1])
    predictions = {}
    slopes = {}
    nbase = rows["base_columns"]
    for model, count in zip(MODELS, (nbase, nbase+1, nbase+2, nbase+3)):
        x = rows["x"][train,:count]
        p, beta = fit_predict(x, rows["y"][train], w, rows["x"][test,:count], penalty)
        predictions[model] = p
        if model == "M2":
            # Same weighted scaling used internally by v304.fit_predict.
            wn = w / w.sum()
            mean = np.sum(x*wn[:,None], axis=0)
            sd = np.sqrt(np.sum((x-mean)**2*wn[:,None], axis=0))
            sd[sd < 1e-10] = 1.
            slopes = foreign_slopes(beta, sd)
    return dict(index=rows["index"][test], y=rows["y"][test],
                predictions=predictions, foreign_slopes=slopes)


def score_habitats(y, p, habitats):
    scores = {}
    for h, label in enumerate(HABITATS):
        chosen = habitats == h
        n = int(chosen.sum())
        scores[label] = dict(n=n, events=int(y[chosen].sum()),
                             event_rate=float(y[chosen].mean()) if n else None,
                             **(metrics(y[chosen],p[chosen]) if n else {"log_loss":None,"brier":None}))
    return {
        "by_habitat": scores,
        "balanced": {key:balanced_average([scores[h][key] for h in HABITATS])
                     for key in ("log_loss","brier")},
    }


def summarize(values):
    available = [v for v in values if v is not None]
    n = len(available)
    leave = [float(np.mean(available[:j]+available[j+1:])) for j in range(n)] if n > 1 else []
    return {
        "mean": float(np.mean(available)) if n else None,
        "valid_blocks": n, "positive_blocks": sum(v > 0 for v in available),
        "block_values": values,
        "leave_one_block_mean_range": [min(leave),max(leave)] if leave else None,
    }


def analyze_endpoint(c, endpoint, penalty):
    folds = []
    for b in range(5):
        fit = forecast_fold(c, endpoint, b, penalty)
        scores = {m:score_habitats(fit["y"],p,fit["index"][:,1])
                  for m,p in fit["predictions"].items()}
        folds.append(dict(block=b, transects=int(np.sum(c["blocks"]==b)),
                          models=scores, foreign_slopes=fit["foreign_slopes"]))
    models = {}
    for model in MODELS:
        models[model] = {}
        for scope in ("balanced",*HABITATS):
            models[model][scope] = {}
            for metric in ("log_loss","brier"):
                values = [(f["models"][model]["balanced"][metric] if scope=="balanced"
                           else f["models"][model]["by_habitat"][scope][metric]) for f in folds]
                models[model][scope][metric] = summarize(values)
    comparisons = {}
    for name,(left,right) in COMPARISONS.items():
        comparisons[name] = {}
        for scope in ("balanced",*HABITATS):
            comparisons[name][scope] = {}
            for metric in ("log_loss","brier"):
                a = models[left][scope][metric]["block_values"]
                b = models[right][scope][metric]["block_values"]
                comparisons[name][scope][metric] = summarize([
                    x-y if x is not None and y is not None else None for x,y in zip(a,b)])
    counts = {}
    for label in HABITATS:
        cells = [f["models"]["M0"]["by_habitat"][label] for f in folds]
        n, events = sum(s["n"] for s in cells), sum(s["events"] for s in cells)
        counts[label] = dict(n=n, events=events, event_rate=events/n if n else None)
    return dict(endpoint=endpoint,penalty=penalty,counts=counts,models=models,
                comparisons=comparisons,
                foreign_slopes={key:summarize([f["foreign_slopes"][key] for f in folds])
                                for key in folds[0]["foreign_slopes"]},folds=folds)


def read_sources(data):
    # These existing parsers validate original source bytes/tables internally.
    import validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 as eu
    from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
    from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping
    from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
    sys.path.insert(0,str(ROOT/"research/validation_20261117"))
    from rail_frequency_sources import names, frequencies
    fixed = json.loads((ROOT/"research/validation_20261117/country_validation.json").read_text())
    expected = {
        "cz1.xlsx":"df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06",
        "cz2.xlsx":"00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2",
        "parma.pdf":eu.PUBLISHER_HASH["italy_original_pdf"],
        "hrba.pdf":eu.PUBLISHER_HASH["hr_ba_original_pdf"],
        "Matrix_plot.xlsx":fixed["sha256"]["Matrix_plot.xlsx"],
    }
    actual_hashes = {}
    for file in ("cz1.xlsx","cz2.xlsx","parma.pdf","hrba.pdf","Matrix_plot.xlsx"):
        actual_hashes[file] = hashlib.sha256((data/file).read_bytes()).hexdigest()
        if actual_hashes[file] != expected[file]:
            raise ValueError("Source byte mismatch: "+file)
    donors, study_audits = [], {}
    loaded = [
        ("CZ39",39,eu.publisher_czech(data/"cz1.xlsx",data/"cz2.xlsx")),
        ("IT28",28,eu.publisher_italy(data/"parma.pdf")),
        ("HR_BA15",15,eu.publisher_hr_ba(data/"hrba.pdf")),
    ]
    for tag,n,result in loaded:
        sites = result[0]
        if len(sites) != n:
            raise ValueError("Unexpected donor site count: "+tag)
        donors.append(frequencies([names(s["source_names"]) for s in sites]))
        study_audits[tag] = {"sites":n, "binomials":len(donors[-1])}
    xlsx = list((data/"si").rglob("*.xlsx"))
    if len(xlsx) != 1:
        raise ValueError("Expected one Slovenian source workbook")
    digest = hashlib.sha256(xlsx[0].read_bytes()).hexdigest()
    if digest != "26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f":
        raise ValueError("Slovenian workbook hash mismatch")
    jsons = [p for p in (data/"si").rglob("*.json") if "Slovenia_2026_RAIL_92sites_plant_supplement_source_v102" in p.name]
    if len(jsons) != 1:
        raise ValueError("Expected original Slovenian source audit")
    original = extract_original_xlsx(xlsx[0],jsons[0])
    sites,_,_ = source_sites_and_guilds(original)
    code,_ = original_slovenian_EPPO_codebook_binomial_mapping(original)
    if len(sites) != 89:
        raise ValueError("Unexpected Slovenian QC site count")
    donors.append(frequencies([names([code.get(t) for t in s["all_original_taxa"]]) for s in sites]))
    study_audits["SI89"] = {"sites":89,"binomials":len(donors[-1])}
    actual_hashes["Slovenia_source.xlsx"] = digest
    actual_hashes["Slovenia_source_audit.json"] = hashlib.sha256(jsons[0].read_bytes()).hexdigest()
    species, plots, meta, k_audit = local_source(data, {"eu":eu})
    c = prepare_data(species,plots,meta,donors)
    return c, {"source_sha256":actual_hashes,"donor_studies":study_audits,
               "kashmir_original_source":k_audit, "cohort":c["audit"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    c, audit = read_sources(args.data.resolve())
    result = {
        "version":"v313",
        "design_commit":"865cac4a834916fa82b346de466696f07bac31b3",
        "interpretation":"Retrospective recorded transitions; not causal transport, survival or historical prospective forecast.",
        "validation_unit":"Five contiguous latitude blocks in one Kashmir corridor; not five independent studies.",
        "source_audit":audit,
        "measurement":"Original species abundance values; no assumption of percent cover.",
        "runtime":{"python":platform.python_version(), **{name:importlib.metadata.version(name)
                   for name in ("numpy","scipy","pandas","openpyxl","PyMuPDF")}},
        "panels":[],
    }
    print(json.dumps(audit["cohort"]),flush=True)
    for endpoint in ("persistence","new_record"):
        for penalty in PENALTIES:
            result["panels"].append(analyze_endpoint(c,endpoint,penalty))
            print(f"Completed {endpoint} lambda={penalty}",flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()

