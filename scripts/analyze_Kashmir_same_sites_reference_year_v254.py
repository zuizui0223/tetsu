#!/usr/bin/env python3
"""v25.4: Same-site reference-specific temporal railway flora beta contrasts.

Uses publisher-authentic Dryad Matrix_plot.xlsx, checks publisher SHA256,
reads XLSX with Python stdlib ZIP+XML, and does NOT redistribute source flora.
Station/between is a spatial stratum, NOT observed train traffic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
import zipfile
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET

DRYAD_SHA256 = "529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474"
AGGREGATE_GIT_BLOB = "75dd78dd2b66a67543c248ff686c3a9dc016ea73"
NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
ZONES = ("Track", "Interior Plot 1", "Interior Plot 2")
YEARS = (2014, 2017)
GROUPS = ("station", "between")
REFERENCES = ("Interior Plot 1", "Interior Plot 2")


def blob_sha1(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\x00" + raw).hexdigest()


def excel_col_index(address: str) -> int:
    m = re.match(r"^([A-Z]+)[0-9]+$", address)
    if m is None:
        raise ValueError(f"Not an XLSX cell reference: {address!r}")
    idx = 0
    for char in m.group(1):
        idx = idx * 26 + (ord(char) - ord("A") + 1)
    return idx - 1


def parse_xlsx_matrix(path: Path) -> tuple[dict[tuple[str, int, str], set[int] | None], dict]:
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != DRYAD_SHA256:
        raise ValueError("Original Dryad Matrix_plot.xlsx SHA-256 failed closed")
    if not data.startswith(b"PK\x03\x04"):
        raise ValueError("Original workbook is not XLSX")

    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{{{NS}}}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{NS}}}t")))
        book = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = book.findall(f"{{{NS}}}sheets/{{{NS}}}sheet")
        if not sheets:
            raise ValueError("No original source sheet")
        relationships = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        sheet_rid = sheets[0].get(f"{{{REL_NS}}}id")
        target = None
        for rel in relationships.findall(f"{{{PKG_NS}}}Relationship"):
            if rel.get("Id") == sheet_rid and "worksheet" in (rel.get("Type") or ""):
                target = rel.get("Target")
                break
        if not target:
            raise ValueError("Original source sheet relation missing")
        target = "xl/" + target.lstrip("/").removeprefix("xl/")
        if target not in z.namelist():
            raise ValueError("Original source worksheet target missing")
        root = ET.fromstring(z.read(target))
        row_nodes = root.findall(f"{{{NS}}}sheetData/{{{NS}}}row")
        rows = []
        for r in row_nodes:
            values: dict[int, object] = {}
            for cell in r.findall(f"{{{NS}}}c"):
                col = excel_col_index(cell.get("r", ""))
                ctype = cell.get("t", "")
                val = cell.find(f"{{{NS}}}v")
                text = val.text if val is not None else None
                if ctype == "s":
                    if text is not None:
                        pos = int(text)
                        if pos >= len(shared):
                            raise ValueError("Shared-string index exceeds original source count")
                        values[col] = shared[pos]
                elif ctype == "inlineStr":
                    values[col] = "".join(x.text or "" for x in cell.iter(f"{{{NS}}}t"))
                elif ctype == "b":
                    values[col] = text == "1" if text is not None else None
                elif ctype in ("n", "", "e"):
                    if text is not None:
                        try:
                            values[col] = float(text)
                        except ValueError:
                            values[col] = text
                elif ctype == "str":
                    values[col] = text
                else:
                    raise ValueError(f"Unknown original source XLSX cell type {ctype}")
            rows.append(values)
    if len(rows) != 249:
        raise ValueError(f"Expected 249 XLSX row nodes incl. header, got {len(rows)}")
    header = rows[0]
    if sorted(header) != list(range(129)):
        raise ValueError("Expected first row with 6 metadata and 123 species columns")
    species_headers = [str(header[c]) for c in range(6, 129)]
    if len(set(species_headers)) != len(species_headers):
        raise ValueError("Duplicated original source species column headers")

    plots = {}
    missing = 0
    all_sites = set()
    for values in rows[1:]:
        site = " ".join(str(values.get(0, "")).split())
        try:
            year = int(values[4])
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError("Source year not an integer") from exc
        zone = " ".join(str(values.get(5, "")).split())
        if not site or zone not in ("Track", "Embankment", *REFERENCES) or year not in YEARS:
            raise ValueError(f"Unexpected source site/year/zone {site}/{year}/{zone}")
        key = (site, year, zone)
        if key in plots:
            raise ValueError("Duplicate original source site/year/zone")
        all_sites.add(site)
        this = set()
        complete = True
        for index in range(6, 129):
            v = values.get(index)
            if v is None:
                missing += 1
                complete = False
                continue
            if not isinstance(v, (float, int)) or not math.isfinite(v) or v < 0:
                raise ValueError("Unexpected source species-by-plot cover value")
            if v > 0:
                this.add(index - 6)
        plots[key] = this if complete and this else None
    if len(all_sites) != 31 or len(plots) != 248:
        raise ValueError("Expected 31 spatial sites x 2 years x 4 zones")
    if missing != 4:
        raise ValueError(f"Expected four source missing species values; found {missing}")
    stations = sum(not s.startswith("Between ") for s in all_sites)
    if stations != 16:
        raise ValueError(f"Expected 16 stations (other 15 between); got {stations}")
    return plots, {
        "original_file_sha256": DRYAD_SHA256,
        "original_file_size": len(data),
        "all_spatial_sites": len(all_sites),
        "station_sites_in_source": stations,
        "between_sites_in_source": len(all_sites)-stations,
        "plot_year_zone_rows": len(plots),
        "source_taxon_columns": len(species_headers),
        "source_missing_species_cells": missing,
        "original_sheet": sheets[0].get("name")
    }


def matching_sites(plots: dict, group: str, zones: tuple[str, ...]) -> list[str]:
    all_sites = {k[0] for k in plots}
    selected = []
    for site in sorted(all_sites):
        if ("between" if site.startswith("Between ") else "station") != group:
            continue
        if all(plots.get((site, y, zone)) is not None for y in YEARS for zone in zones):
            selected.append(site)
    return selected


def jaccard(a: set[int], b: set[int]) -> float:
    union = len(a | b)
    if union == 0:
        raise ValueError("Jaccard undefined for two empty species sets")
    return 1.0 - len(a & b) / union


def pair_distance(plots: dict, cohort: dict[str, list[str]], *, seed: int | None = None) -> dict:
    """Prepare observed / once-rarefied station-pair distances by year and zone."""
    rng = random.Random(seed) if seed is not None else None
    result = {}
    for group, sites in cohort.items():
        result[group] = {}
        for year in YEARS:
            sets = {}
            for site in sites:
                this = {z: plots[site, year, z] for z in ZONES}
                if any(x is None or len(x) == 0 for x in this.values()):
                    raise ValueError("Missing site entering common-cohort analysis")
                target = min(len(x) for x in this.values())
                for zone, species in this.items():
                    sets[site, zone] = (frozenset(rng.sample(sorted(species), target)) if rng is not None
                                        else frozenset(species))
            result[group][year] = {}
            for zone in ZONES:
                pairs = {}
                for i, a in enumerate(sites):
                    for b in sites[i+1:]:
                        pairs[a, b] = jaccard(sets[a, zone], sets[b, zone])
                result[group][year][zone] = pairs
    return result


def averaged_pair_distances(plots: dict, cohort: dict[str, list[str]], draws: int = 199) -> dict:
    acc = None
    for i in range(draws):
        curr = pair_distance(plots, cohort, seed=254019 + 107 * i)
        if acc is None:
            acc = curr
            for group in GROUPS:
                for year in YEARS:
                    for zone in ZONES:
                        acc[group][year][zone] = {key: 0.0 for key in curr[group][year][zone]}
        for group in GROUPS:
            for year in YEARS:
                for zone in ZONES:
                    for key, val in curr[group][year][zone].items():
                        acc[group][year][zone][key] += val / draws
    return acc


def mean_pair(d: dict[tuple[str, str], float], site_weights: dict[str, int] | None = None) -> float:
    num = 0.0
    denom = 0.0
    for (a, b), value in d.items():
        wt = site_weights[a] * site_weights[b] if site_weights is not None else 1
        num += wt * value
        denom += wt
    if denom == 0:
        raise ValueError("No distinct original site pair in resampled data")
    return num / denom


def calc(pair: dict, weights: dict[str, dict[str, int]] | None = None) -> dict:
    means = {}
    for g in GROUPS:
        means[g] = {}
        for y in YEARS:
            means[g][str(y)] = {
                zone: mean_pair(pair[g][y][zone], weights.get(g) if weights else None)
                for zone in ZONES
            }
    by_ref = {}
    for ref in REFERENCES:
        changes = {}
        for g in GROUPS:
            r = means[g]
            changes[g] = (r["2017"]["Track"] - r["2017"][ref]) - (r["2014"]["Track"] - r["2014"][ref])
        by_ref[ref] = {
            "station_change": changes["station"],
            "between_change": changes["between"],
            "station_minus_between_change": changes["station"] - changes["between"],
        }
    actual = (by_ref["Interior Plot 1"]["station_minus_between_change"] -
              by_ref["Interior Plot 2"]["station_minus_between_change"])
    # Track cancels exactly when the same sites and zone sampling are used.
    expected = ((means["station"]["2017"]["Interior Plot 2"] - means["station"]["2014"]["Interior Plot 2"])-
                (means["station"]["2017"]["Interior Plot 1"] - means["station"]["2014"]["Interior Plot 1"])-
                (means["between"]["2017"]["Interior Plot 2"] - means["between"]["2014"]["Interior Plot 2"])+
                (means["between"]["2017"]["Interior Plot 1"] - means["between"]["2014"]["Interior Plot 1"]))
    if abs(actual - expected) > 1e-11:
        raise AssertionError("Track cancellation invariant failed")
    return {"means": means, "by_reference": by_ref, "reference_sensitivity_station_minus_between": actual,
            "track_cancels_in_reference_sensitivity": True}


def percentile(values: list[float], p: float) -> float:
    x = sorted(values)
    k = (len(x)-1) * p
    lo = int(math.floor(k)); hi = int(math.ceil(k))
    return x[lo] + (k-lo)*(x[hi]-x[lo])


def sample_weights(sites: list[str], rng: random.Random) -> dict[str, int]:
    c = Counter(rng.choices(sites, k=len(sites)))
    return {s: c.get(s, 0) for s in sites}


def resample_summary(pair: dict, cohort: dict[str,list[str]], draws: int = 1999, seed: int = 254001) -> dict:
    rng = random.Random(seed)
    by_ref = {ref: [] for ref in REFERENCES}
    sensitivity = []
    skipped = 0
    for _ in range(draws):
        weights = {g: sample_weights(cohort[g], rng) for g in GROUPS}
        try:
            out = calc(pair, weights)
        except ValueError:
            skipped += 1
            continue
        for ref in REFERENCES:
            by_ref[ref].append(out["by_reference"][ref]["station_minus_between_change"])
        sensitivity.append(out["reference_sensitivity_station_minus_between"])
    def finish(xs):
        if len(xs) < draws * 0.95:
            raise ValueError("Too many degenerate spatial resamples; insufficient support")
        return {"draws":len(xs),"percentile_025_technical_not_field_CI":percentile(xs,.025),
                "percentile_975_technical_not_field_CI":percentile(xs,.975),
                "fraction_negative":sum(x < 0 for x in xs)/len(xs)}
    jackknife = {}
    for g in GROUPS:
        vals = {ref:[] for ref in REFERENCES}
        for omitted in cohort[g]:
            weights = {gg:{s: int(gg!=g or s!=omitted) for s in cohort[gg]} for gg in GROUPS}
            output = calc(pair, weights)
            for ref in REFERENCES:
                vals[ref].append(output["by_reference"][ref]["station_minus_between_change"])
        jackknife[g] = {ref: {"min":min(vals[ref]),"max":max(vals[ref])} for ref in REFERENCES}
    return {"method":"site-cluster multiplier bootstrap: spatial sites resampled as paired multi-zone multiyear clusters; same-site pairs excluded by U-statistic weighting",
            "bootstrap_requested":draws,"bootstrap_skipped_degenerate":skipped,
            "by_reference":{ref:finish(xs) for ref,xs in by_ref.items()},
            "reference_sensitivity":finish(sensitivity),
            "leave_one_original_site_out_by_group":jackknife,
            "percentiles_are_technical_site_resampling_not_randomized_field_CI":True}



def within_site_temporal_turnover(plots: dict, cohort: dict[str,list[str]],
                                  draws: int = 199, bootstrap: int = 1999) -> dict:
    """2014-vs-2017 taxonomic turnover within each *original spatial site*.

    Raw analysis uses all recorded taxa. Technical rarefaction reduces all six
    site-year-zone lists to a common minimum richness at each site.
    """
    rng = random.Random(254099)
    score = {g: {} for g in GROUPS}
    for g in GROUPS:
        for site in cohort[g]:
            raw = {z: jaccard(plots[site,2014,z], plots[site,2017,z]) for z in ZONES}
            min_n = min(len(plots[site,y,z]) for y in YEARS for z in ZONES)
            if min_n < 1:
                raise ValueError("No positive species in common site/year/zone")
            rare = {z: 0.0 for z in ZONES}
            for _ in range(draws):
                for z in ZONES:
                    aa=set(rng.sample(sorted(plots[site,2014,z]), min_n))
                    bb=set(rng.sample(sorted(plots[site,2017,z]), min_n))
                    rare[z]+=jaccard(aa,bb)/draws
            score[g][site] = {"raw":raw,"rarefied":rare}

    def evaluate(mode: str, weights: dict[str,dict[str,int]] | None = None) -> dict:
        means = {}
        for g in GROUPS:
            w=weights[g] if weights else {s:1 for s in cohort[g]}
            denominator=sum(w.values())
            if denominator==0:raise ValueError("Empty turnover site bootstrap")
            means[g]={z:sum(w[s]*score[g][s][mode][z] for s in cohort[g])/denominator for z in ZONES}
        out={}
        for ref in REFERENCES:
            out[ref]={
                "station_track_minus_reference":means["station"]["Track"]-means["station"][ref],
                "between_track_minus_reference":means["between"]["Track"]-means["between"][ref],
                "station_minus_between_track_reference_turnover_difference":
                  (means["station"]["Track"]-means["station"][ref])-
                  (means["between"]["Track"]-means["between"][ref])
            }
        return {"by_group_zone":means,"by_reference":out}

    out={}
    for mode,seed in (("raw",254711),("rarefied",254712)):
        observed=evaluate(mode)
        rand=random.Random(seed)
        draws_ref={ref:[] for ref in REFERENCES}
        for _ in range(bootstrap):
            w={g:sample_weights(cohort[g],rand) for g in GROUPS}
            b=evaluate(mode,w)
            for ref in REFERENCES:
                draws_ref[ref].append(b["by_reference"][ref]["station_minus_between_track_reference_turnover_difference"])
        out[mode]={"observation":observed,
                   "site_bootstrap_technical":{
                       ref:{"draws":len(v),"q025_not_field_CI":percentile(v,.025),
                            "q975_not_field_CI":percentile(v,.975),
                            "fraction_negative":sum(x<0 for x in v)/len(v)}
                       for ref,v in draws_ref.items()}}
    return {"interpretation":"Within-site Jaccard dissimilarity between original 2014 and 2017 surveys, unlike between-site beta convergence",
            "rarefaction":"all six site-year-zone species sets rarefied to minimum richness within that site; %d independently sampled draws"%draws,
            "not_traffic_or_heritable_adaptation":True,
            "results":out}

def pairwise_reproduce_previous(plots: dict, verified: dict) -> dict:
    matched = []
    for group in GROUPS:
        for zone, ref in (("Track", "Interior Plot 2"), ("Track", "Interior Plot 1"),
                          ("Embankment", "Interior Plot 2")):
            sites = matching_sites(plots, group, (zone, ref))
            for year in YEARS:
                def avg_distance(x: str) -> float:
                    return sum(jaccard(plots[s,year,x],plots[t,year,x])
                               for i,s in enumerate(sites) for t in sites[i+1:]) / (len(sites)*(len(sites)-1)/2)
                delta = avg_distance(zone) - avg_distance(ref)
                prior = next((v for v in verified["all_species"] if v["group"] == group
                              and v["zones"] == [zone,ref] and v["year"] == year),None)
                if prior is None or len(sites) != prior["n"] or abs(delta-prior["all_species_distance_delta"]) > 1e-10:
                    raise ValueError(f"Earlier source v25.3 all-species raw aggregate not reproduced: {group}/{zone}/{ref}/{year}")
                matched.append({"group":group,"rail_zone":zone,"reference_zone":ref,
                                "year":year,"source_site_count":len(sites),"nonrarefied_delta_exactly_reproduced":delta})
    return {"twelve_source_aggregates_reproduced":len(matched),"all_passed_1e_minus_10":True,
            "same_original_sites_used_in_both_years_by_construction":True,
            "strata":matched}


def audit(plots: dict, verified: dict, source_meta: dict, bootstrap: int = 1999, rarefaction: int = 199) -> dict:
    reproduced = pairwise_reproduce_previous(plots,verified)
    sites = {g: matching_sites(plots,g,ZONES) for g in GROUPS}
    if min(len(x) for x in sites.values()) < 5:
        raise ValueError("At least five original spatial sites per stratum needed for common-cohort analysis")
    observed_pairs = pair_distance(plots,sites)
    observed = calc(observed_pairs)
    rare_pairs = averaged_pair_distances(plots,sites,draws=rarefaction)
    rare = calc(rare_pairs)
    return {
        "version":"v25.4",
        "source":source_meta,
        "prior_derived_aggregate_git_blob_sha1":AGGREGATE_GIT_BLOB,
        "source_reproduction":reproduced,
        "same_site_common_3_zone_cohort":{g:len(sites[g]) for g in GROUPS},
        "sampling_unit":"original spatial site; same retained site IDs in 2014 and 2017, across Track and both Interiors",
        "nonrarefied":observed,
        "rarefaction":"all three zones within each site-year independently rarefied to minimum within-site-year richness; 199 draw pairwise distances averaged; same sites & cohorts",
        "rarefied":rare,
        "cluster_bootstrap_raw":resample_summary(observed_pairs,sites,draws=bootstrap,seed=254001),
        "cluster_bootstrap_rarefied":resample_summary(rare_pairs,sites,draws=bootstrap,seed=254002),
        "within_site_taxonomic_turnover_2014_to_2017":within_site_temporal_turnover(plots,sites,draws=rarefaction,bootstrap=bootstrap),
        "key_causal_limits":["rail traffic exposure not observed", "station/between not randomized", "sites may vary in management, region and climate", "technical bootstrap is not field experiment", "source XLSX/plot identities intentionally not republished"]
    }


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",required=True,type=Path)
    parser.add_argument("--prior-aggregate",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--bootstrap",type=int,default=1999)
    parser.add_argument("--rarefaction",type=int,default=199)
    args=parser.parse_args()
    raw=args.prior_aggregate.read_bytes()
    if blob_sha1(raw)!=AGGREGATE_GIT_BLOB:
        raise SystemExit("The prior independent source-derived aggregate changed: fail closed")
    plots,meta=parse_xlsx_matrix(args.source)
    report=audit(plots,json.loads(raw),meta,bootstrap=args.bootstrap,rarefaction=args.rarefaction)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"STATUS":"source_verified","common_site_cohort":report["same_site_common_3_zone_cohort"],
                      "raw":{k:v["station_minus_between_change"] for k,v in report["nonrarefied"]["by_reference"].items()},
                      "rarefied":{k:v["station_minus_between_change"] for k,v in report["rarefied"]["by_reference"].items()}},ensure_ascii=False))


if __name__=="__main__":
    main()
