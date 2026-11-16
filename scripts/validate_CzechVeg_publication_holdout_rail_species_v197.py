#!/usr/bin/env python3
"""v19.7: honest original-study holdout of Czech railway-text species filters.

The v17.4 whole-source rank was selected using the same CzechVeg study
collection that supplied its leave-one-author-out checks. Here, 36 original
source *publications* are allocated to four folds using plot opportunity
counts ONLY (no plant identities or outcome signs). For each fold, a new
positive/negative species rank is learned from the other publications, then
scored on untouched publication blocks at exactly the original same-study x
phytoregion x decade x area matching grain.

Study holdout is internal prediction, NOT external railway-versus-natural
validation or evidence of train dispersal, trait evolution, or a field p
value. Original sources are downloaded in ephemeral Actions workspaces and
validated by the v17.4 publisher MD5/SHA checks; only aggregate summaries
are written to GitHub Actions artifacts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    verify_file, plot_labels, HEADERS_NAME, SPECIES_NAME,
)
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import header_source_context
from analyze_CzechVeg_rail_road_same_study_trait_v123 import source_study_meta
from analyze_original_CzechVeg_independent_railroad_plot_baseline_39rail_stations_v174 import (
    station_source, source_vegetation_scan, retain_blocks, standardized_taxon_panel,
    MIN_ORIGINAL_STATION_OCCURRENCES,
    MIN_BACKGROUND_MATCHED_PLOT_SPECIES_INCIDENCES,
    MIN_DISTINCT_ORIGINAL_BOTANICAL_STUDIES_WITH_TAXON,
)

K_FOLDS = 4
TOP_N = 20
LOCKED_FROM_V174 = (
    "Echium vulgare",
    "Poa compressa",
    "Artemisia vulgaris",
    "Conyza canadensis",
)
SCORE = "CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence"
RAIL_N = "CzechVeg_same_study_region_decade_area_source_N_RAIL_text_plot_with_this_taxon"
ROAD_N = "CzechVeg_same_study_region_decade_area_source_N_ROAD_text_plot_with_this_taxon"
PAPER_N = "CzechVeg_same_study_source_N_distinct_publications_with_taxon"
STATION_N = "station_original_2020_22_n_of_39"


def split_source_studies(blocks, nfold=K_FOLDS):
    """Assign whole source studies, balancing only pre-outcome survey opportunity."""
    weights = Counter()
    for key, arms in blocks.items():
        weights[key[0]] += min(len(arms["rail"]), len(arms["road"]))
    if len(weights) < nfold * 2:
        raise ValueError("Too few separate original source studies for study-blocked CV")
    cap = math.ceil(len(weights) / nfold)
    fold_studies = [set() for _ in range(nfold)]
    masses = [0] * nfold
    order = sorted(weights, key=lambda name: (
        -weights[name], hashlib.sha256(name.encode("utf-8")).hexdigest()
    ))
    for name in order:
        possible = [i for i in range(nfold) if len(fold_studies[i]) < cap]
        j = min(possible, key=lambda i: (masses[i], len(fold_studies[i]), i))
        fold_studies[j].add(name)
        masses[j] += weights[name]
    assert set.union(*fold_studies) == set(weights)
    assert sum(len(q) for q in fold_studies) == len(weights)
    return tuple(frozenset(s) for s in fold_studies), masses


def supported(taxon_result):
    """Minimum original study and plot support, irrespective the effect sign."""
    x = taxon_result or {}
    return (
        x.get(SCORE) is not None
        and x.get(STATION_N, 0) >= MIN_ORIGINAL_STATION_OCCURRENCES
        and x.get(RAIL_N, 0) + x.get(ROAD_N, 0)
        >= MIN_BACKGROUND_MATCHED_PLOT_SPECIES_INCIDENCES
        and x.get(PAPER_N, 0) >= MIN_DISTINCT_ORIGINAL_BOTANICAL_STUDIES_WITH_TAXON
    )


def rank_train(panel, sign, n=TOP_N):
    if sign not in ("positive", "negative"):
        raise ValueError("Must explicitly choose positive or negative direction")
    eligible = [(name, x[SCORE]) for name, x in panel.items() if supported(x)]
    eligible = [(name, score) for name, score in eligible
                if (score > 0 if sign == "positive" else score < 0)]
    # Rank is entirely training-only: no heldout taxon presence or score used.
    eligible.sort(key=lambda item: (
        -item[1] if sign == "positive" else item[1], item[0]
    ))
    return [name for name, _ in eligible[:n]]


def heldout_score(panel, name):
    val = panel.get(name)
    if not supported(val):
        return None
    return val[SCORE]


def validate_rank(train, test, direction, n=TOP_N):
    selected = rank_train(train, direction, n=n)
    valid = [heldout_score(test, name) for name in selected]
    valid = [score for score in valid if score is not None]
    favorable = sum(score > 0 if direction == "positive" else score < 0
                    for score in valid)
    return {
        "training_selected_n": len(selected),
        "heldout_evaluable_n": len(valid),
        "heldout_missing_min_source_support_n": len(selected) - len(valid),
        "heldout_direction_replicated_n": favorable,
        "heldout_opposite_or_zero_n": len(valid) - favorable,
        "heldout_direction_replicated_fraction_DESCRIPTIVE_ONLY": (
            favorable / len(valid) if valid else None
        ),
        "heldout_median_weighted_plot_prevalence_difference": (
            round(statistics.median(valid), 7) if valid else None
        ),
        "selection_never_used_heldout_outcomes": True,
    }


def fold_diagnostic(blocks, plot_taxa, station, background, denom):
    groups, opportunity = split_source_studies(blocks)
    full_studies = set.union(*(set(x) for x in groups))
    folds = []
    for i, unseen in enumerate(groups):
        train_blocks = {k: v for k, v in blocks.items() if k[0] not in unseen}
        test_blocks = {k: v for k, v in blocks.items() if k[0] in unseen}
        assert train_blocks and test_blocks
        assert {k[0] for k in train_blocks}.isdisjoint({k[0] for k in test_blocks})
        assert ({k[0] for k in train_blocks} | {k[0] for k in test_blocks}) == full_studies
        train, tq = standardized_taxon_panel(train_blocks, plot_taxa, station, background, denom)
        test, vq = standardized_taxon_panel(test_blocks, plot_taxa, station, background, denom)
        frozen = {}
        for taxon in LOCKED_FROM_V174:
            tr = train.get(taxon)
            te = test.get(taxon)
            frozen[taxon] = {
                "training_effect": round(tr[SCORE], 7) if supported(tr) else None,
                "heldout_effect": round(te[SCORE], 7) if supported(te) else None,
                "heldout_source_studies_with_this_taxon": te.get(PAPER_N, 0) if te else 0,
                "heldout_taxon_observed_plots": (
                    te.get(RAIL_N, 0) + te.get(ROAD_N, 0) if te else 0
                ),
                "heldout_positive_if_evaluable": (
                    te[SCORE] > 0 if supported(te) else None
                ),
                "previously_selected_taxon_NOT_external_confirmatory_hypothesis": True,
            }
        folds.append({
            "fold": i + 1,
            "study_n_train": tq["n_original_source_field_studies"],
            "study_n_heldout": vq["n_original_source_field_studies"],
            "study_opportunity_mass_heldout": opportunity[i],
            "train_original_blocks_n": len(train_blocks),
            "heldout_original_blocks_n": len(test_blocks),
            "heldout_original_RAIL_text_plot_n": vq["n_original_rail_text_plots"],
            "heldout_original_ROAD_text_plot_n": vq["n_original_road_text_plots"],
            "train_selected_positive": validate_rank(train, test, "positive"),
            "train_selected_negative_placebo": validate_rank(train, test, "negative"),
            "locked_v174_taxa": frozen,
        })
    return folds


def summarize(folds):
    ans = {}
    for key in ("train_selected_positive", "train_selected_negative_placebo"):
        evaluations = sum(f[key]["heldout_evaluable_n"] for f in folds)
        concordant = sum(f[key]["heldout_direction_replicated_n"] for f in folds)
        ans[key] = {
            "across_folds_total_evaluable_taxon_fold_pairs_NONINDEPENDENT": evaluations,
            "direction_replicated_taxon_fold_pairs": concordant,
            "fraction_descriptive_NOT_field_p": concordant / evaluations if evaluations else None,
            "per_fold_fractions": [
                f[key]["heldout_direction_replicated_fraction_DESCRIPTIVE_ONLY"]
                for f in folds
            ],
        }
    ans["locked_prior_v174_taxa"] = {}
    for taxon in LOCKED_FROM_V174:
        v = [f["locked_v174_taxa"][taxon]["heldout_effect"] for f in folds]
        found = [x for x in v if x is not None]
        ans["locked_prior_v174_taxa"][taxon] = {
            "heldout_supported_folds": len(found),
            "heldout_positive_folds": sum(x > 0 for x in found),
            "heldout_effects_in_fold_order_MISSING_NULL": v,
            "a_fold_is_NOT_an_independent_external_railway_study": True,
        }
    return ans


def run(zenodo, headers, species, s1, s2):
    original_headers = verify_file(zenodo, HEADERS_NAME, headers)
    original_species = verify_file(zenodo, SPECIES_NAME, species)
    station, station_qa = station_source(s1, s2)
    plotmeta, label_qa = plot_labels(headers)
    yearmeta, year_qa = header_source_context(headers, plotmeta)
    meta, study_qa = source_study_meta(headers, yearmeta)
    candidate = {
        pid for pid, row in meta.items()
        if row["broad"] in ("rail", "road") or row["strict"] in ("rail", "road")
    }
    plot_taxa, background, scan, denom = source_vegetation_scan(
        species, plotmeta, candidate, set(station)
    )
    blocks, match_qa = retain_blocks(meta, plot_taxa, "broad", None)
    n_rail = sum(len(v["rail"]) for v in blocks.values())
    n_road = sum(len(v["road"]) for v in blocks.values())
    n_studies = len({k[0] for k in blocks})
    if (len(blocks), n_studies, n_rail, n_road) != (44, 36, 229, 462):
        raise ValueError("Original full matched frame changed; fail closed")
    folds = fold_diagnostic(blocks, plot_taxa, station, background, denom)
    return {
        "version": "v19.7_original_publication_holdout",
        "original_2026_CzechVeg_MD5": original_headers,
        "original_2026_CzechVeg_species_MD5": original_species,
        "source_Preslia_2024_39station_verified": station_qa,
        "original_plant_rows_source_counts": scan,
        "publication_region_decade_plotarea_source_check": {
            "matched_blocks": len(blocks),
            "original_publications": n_studies,
            "rail_text_plots": n_rail,
            "road_text_plots": n_road,
            "source_year_region_study_gates": match_qa,
            "study_qa_has_orig_study_reference": bool(study_qa),
            "plot_label_source_checked": bool(label_qa),
            "plot_year_source_checked": bool(year_qa),
        },
        "CV_folds_selected_by_study_plot_opportunity_only": folds,
        "summary": summarize(folds),
        "no_heldout_study_leakage": True,
        "not_out_of_source_external_replication": True,
        "no_rail_transport_selection_adaptation_causal_effect": True,
        "no_valid_field_p_or_confidence_interval": True,
        "source_labels_are_historical_RAIL_ROAD_text_not_true_ballast": True,
        "prior_v174_top_taxa_selected_from_entire_CzechVeg_not_confirmatory": True,
    }


def main():
    p = argparse.ArgumentParser()
    for arg in ("zenodo-record", "source-CzechVeg-headers", "source-CzechVeg-species",
                "Preslia-original-S1", "Preslia-original-S2", "output"):
        p.add_argument("--" + arg, required=True)
    a = p.parse_args()
    obj = run(
        json.loads(Path(a.zenodo_record).read_text(encoding="utf-8")),
        a.source_CzechVeg_headers, a.source_CzechVeg_species,
        a.Preslia_original_S1, a.Preslia_original_S2,
    )
    target = Path(a.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "version": obj["version"],
        "original_source": obj["publication_region_decade_plotarea_source_check"],
        "heldout_fold_source_summary": [
            {k: v for k, v in fold.items() if k in (
                "fold", "study_n_train", "study_n_heldout",
                "train_selected_positive", "train_selected_negative_placebo",
                "locked_v174_taxa"
            )} for fold in obj["CV_folds_selected_by_study_plot_opportunity_only"]
        ],
        "descriptive_heldout_results": obj["summary"],
        "no_independent_field_p": True,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
