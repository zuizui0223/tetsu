#!/usr/bin/env python3
"""v20.0: falsify broad railway-text functional prediction with label/era gates.

Original-source CzechVeg-Open (not the stratified national sample), and original
Pladias traits, are downloaded and verified by the existing v19.8 pipeline.
The primary check is a strict HABITAT-field-only RAIL/ROAD classification
instead of the broad HABITAT-or-LOCALITY classification. Each classification
is additionally repeated in observation-year cohorts, never publication year.

All matched blocks require original publication x phytogeographical region x
survey decade x measured plot-area band support. Heldout folds are complete
original publication titles, assigned using sampling opportunity alone BEFORE
evaluating plant community outcomes. Underpowered cohorts are reported as
NOT IDENTIFIED, never a zero effect. Country=only Czechia; RAIL/ROAD textual
labels do not prove active ballast, management, seed carriage, or adaptation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import random
import statistics
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import plot_labels
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import header_source_context
from analyze_CzechVeg_rail_road_same_study_trait_v123 import (
    common_original_study_blocks, source_study_meta,
)
from validate_CzechVeg_publication_holdout_rail_species_v197 import split_source_studies
from validate_CzechVeg_source_publication_holdout_guild_vs_taxon_v198 import (
    MATCHING, PANEL_THRESHOLDS, frame_records, source_blocks_and_traits,
)
from validate_CzechVeg_heldout_life_diaspore_trait_label_placebo_v199 import (
    MODELS, source_heldout_auc,
)

LABELS = ("broad", "strict")
ERAS = ("all_years", "before_2000", "year_2000_onward", "year_2010_onward")
MIN_PUBLICATIONS = 8
MIN_RAIL_PLOTS = 30
MIN_ROAD_PLOTS = 30
MIN_REGIONS = 3
N_STUDY_BOOTSTRAPS = 2000
SEED = 20261116


def year_in_era(year, era):
    if not isinstance(year, int):
        raise TypeError("Survey year must be an original integer year")
    if era == "all_years":
        return True
    if era == "before_2000":
        return year < 2000
    if era == "year_2000_onward":
        return year >= 2000
    if era == "year_2010_onward":
        return year >= 2010
    raise ValueError("Unregistered original survey-year cohort")


def matched_label_era_blocks(meta, plants, label, era):
    if label not in LABELS or era not in ERAS:
        raise ValueError("Non-preregistered source label/era sensitivity")
    eligible = {
        pid: rec for pid, rec in meta.items()
        if pid in plants and rec[label] in ("rail", "road")
        and year_in_era(rec["source_year_int"], era)
    }
    blocks, quality = common_original_study_blocks(eligible, plants, label, MATCHING)
    return blocks, {
        "candidate_original_plots_with_nonempty_vascular_records_and_year":
            len(eligible),
        "candidate_original_rail_plots": sum(
            rec[label] == "rail" for rec in eligible.values()),
        "candidate_original_road_plots": sum(
            rec[label] == "road" for rec in eligible.values()),
        "source_matching": quality,
    }


def study_bootstrap_descriptive_by_paper(values, n=N_STUDY_BOOTSTRAPS, seed=SEED):
    """Equal-paper resampling diagnostics, NOT causal or independent-country CI."""
    v = list(values)
    if len(v) < MIN_PUBLICATIONS:
        raise ValueError("Descriptive bootstrap requires original study support gate")
    rng = random.Random(seed)
    a = []
    for _ in range(n):
        a.append(statistics.mean(v[rng.randrange(len(v))] for _ in v))
    a.sort()
    return {
        "mean_original_publication_equal": statistics.mean(v),
        "percentile_2_5_study_resampling_NOT_independent_field_CI":
            a[int(.025 * n)],
        "percentile_97_5_study_resampling_NOT_independent_field_CI":
            a[min(n-1, int(.975 * n))],
        "n_original_studies": len(v),
        "n_descriptive_study_resamples": n,
    }


def run_panel(blocks, plants, panel):
    # The original plot richness/trait-known eligibility rule is determined
    # without fitting the railway-vs-road outcome. Fold assignment is based on
    # pre-eligibility publication opportunity, NOT on fitted trait effects.
    selected, sites, original, frame = frame_records(blocks, plants, panel)
    papers = {key[0] for key in selected}
    rail_n = sum(len(v["rail"]) for v in selected.values())
    road_n = sum(len(v["road"]) for v in selected.values())
    regions = {str(key[1]) for key in selected}
    support = {
        **frame,
        "n_retained_original_publications": len(papers),
        "n_retained_original_phytoregions": len(regions),
        "n_retained_source_rail_plots": rail_n,
        "n_retained_source_road_plots": road_n,
        "minimum_original_publications": MIN_PUBLICATIONS,
        "minimum_original_plots_per_arm": min(MIN_RAIL_PLOTS, MIN_ROAD_PLOTS),
        "minimum_original_phytoregions": MIN_REGIONS,
    }
    reasons = []
    if len(papers) < MIN_PUBLICATIONS:
        reasons.append("under_8_original_publications")
    if rail_n < MIN_RAIL_PLOTS:
        reasons.append("under_30_rail_text_source_plots")
    if road_n < MIN_ROAD_PLOTS:
        reasons.append("under_30_road_text_source_plots")
    if len(regions) < MIN_REGIONS:
        reasons.append("under_3_source_phytogeographic_regions")
    if len({k[0] for k in blocks}) < MIN_PUBLICATIONS:
        reasons.append("cannot_split_4_folds_with_2_source_studies_minimum")
    if reasons:
        return {"identified": False, "source_support": support,
                "reasons": reasons,
                "underpowered_not_zero_effect": True}
    folds, masses = split_source_studies(blocks)
    res = {mode: source_heldout_auc(
        selected, sites, original, folds, mode) for mode in MODELS}
    if any(not v["identified"] for v in res.values()):
        return {"identified": False, "source_support": support,
                "reasons": ["original_paper_folds_not_evaluable"],
                "underpowered_not_zero_effect": True}
    paper_ids = set(res["effort"]["original_source_paper_AUC_by_hash"])
    if any(set(v["original_source_paper_AUC_by_hash"]) != paper_ids
           for v in res.values()):
        raise AssertionError("Original papers differ between source models")
    if len(paper_ids) != len(papers):
        raise AssertionError("Source paper holdout silently lost a publication")
    perpaper = {
        p: {mode: res[mode]["original_source_paper_AUC_by_hash"][p]
            for mode in MODELS}
        for p in sorted(paper_ids)
    }
    differences = {}
    for high, low in (
        ("lifeform", "effort"),
        ("dispersal", "effort"),
        ("combined", "effort"),
        ("combined", "lifeform"),
    ):
        name = high + "_minus_" + low
        differences[name] = study_bootstrap_descriptive_by_paper(
            [a[high] - a[low] for a in perpaper.values()],
            seed=SEED + len(perpaper) * 37 + len(name))
    for mode in ("lifeform", "combined"):
        name = mode + "_minus_constant_chance_auc_0_5"
        differences[name] = study_bootstrap_descriptive_by_paper(
            [a[mode] - .5 for a in perpaper.values()],
            seed=SEED + len(perpaper) * 41 + len(name))
    return {
        "identified": True,
        "source_support": support,
        "study_split_opportunity_mass": masses,
        "publication_holdout_fold_size": [len(x) for x in folds],
        "original_publication_equal_holdout_AUC": {
            mode: res[mode]["paper_balanced_heldout_AUC"] for mode in MODELS
        },
        "untrained_constant_score_baseline_AUC": .5,
        "descriptive_study_resampling_paired_AUC": differences,
        "per_study_AUC_only_hashed_original_study_titles": perpaper,
        "n_original_publications_lifeform_auc_gt_chance": sum(
            v["lifeform"] > .5 for v in perpaper.values()),
        "n_original_publications_lifeform_gt_effort": sum(
            v["lifeform"] > v["effort"] for v in perpaper.values()),
        "same_original_publications_for_all_models": True,
        "publication_fold_assignment_pre_model_using_only_plot_opportunity": True,
        "source_text_railway_not_verified_ballast": True,
        "bootstrap_does_not_correct_geographic_or_study_selection_bias": True,
    }


def assert_original_broad_source_frame(rebuilt, frozen):
    if set(rebuilt) != set(frozen):
        raise AssertionError("v19.9 original matched source blocks drifted")
    for key in frozen:
        for arm in ("rail", "road"):
            if set(rebuilt[key][arm]) != set(frozen[key][arm]):
                raise AssertionError("v19.9 source broad actual plot frame drifted")


def run(zenodo, headers, species, diaspore, lifeform):
    broad_original, plants, source = source_blocks_and_traits(
        zenodo, headers, species, diaspore, lifeform)
    plotmeta, label_census = plot_labels(headers)
    years, year_census = header_source_context(headers, plotmeta)
    meta, study_census = source_study_meta(headers, years)
    original_label_cross = Counter(
        (rec["strict"] or "none", rec["broad"] or "none")
        for pid, rec in meta.items() if pid in plants
    )
    out = {}
    for label in LABELS:
        out[label] = {}
        for era in ERAS:
            blocks, candidate = matched_label_era_blocks(
                meta, plants, label, era)
            if label == "broad" and era == "all_years":
                assert_original_broad_source_frame(blocks, broad_original)
            out[label][era] = {
                "original_label_and_survey_year_source_opportunity": candidate,
                "original_publication_holdout": {
                    panel: run_panel(blocks, plants, panel)
                    for panel in PANEL_THRESHOLDS
                },
            }
    # Fail closed if the old broad/all-years sample is not exactly the
    # published source-reproduced v19.9. This is a reproduction check only.
    broad = out["broad"]["all_years"]["original_publication_holdout"]
    expected = {
        "all_vascular_min5": {"effort": .4332111,
                              "lifeform": .6003338,
                              "combined": .6035221},
        "both_trait50_min5": {"effort": .4465063,
                              "lifeform": .5899882,
                              "combined": .6098787},
    }
    for panel, reference in expected.items():
        result = broad[panel]
        if not result["identified"]:
            raise AssertionError("Known v19.9 all-years source baseline lost support")
        observed = result["original_publication_equal_holdout_AUC"]
        for key, value in reference.items():
            if abs(observed[key] - value) > .00002:
                raise AssertionError("True v19.9 source AUC reproduction drift: "
                                     + panel + "/" + key)
    return {
        "version": "v20.0_CzechVeg_original_strict_HABITAT_and_survey_year_holdout",
        "full_original_117739_headers_2510340_species_and_Pladias_source_QA":
            source,
        "source_labels_original_source_census": label_census,
        "source_year_quality": year_census,
        "source_original_study_metadata_quality": study_census,
        "original_plots_with_actual_vascular_records_strict_x_broad_label_counts":
            {a + "_vs_" + b: n for (a,b),n in sorted(original_label_cross.items())},
        "v19.9_broad_all_years_original_baseline_reproduced": True,
        "results": out,
        "country": "Czechia_only",
        "source_text_not_active_rail_ballast_or_source_rand_sample": True,
        "source_pladias_taxon_trait_data_not_redistributed": True,
        "era_cutoffs_use_original_YEAR_not_reference_publication_year": True,
        "potential_Pladias_study_and_geo_dependence_not_causally_removed": True,
        "no_inference_about_global_homogenization_or_evolved_adaptation": True,
    }


def main():
    parser = argparse.ArgumentParser()
    for field in ("zenodo-record", "headers", "species", "diaspore-xlsx",
                  "lifeform-xlsx", "output"):
        parser.add_argument("--" + field, required=True)
    args = parser.parse_args()
    result = run(
        json.loads(Path(args.zenodo_record).read_text(encoding="utf8")),
        args.headers, args.species, args.diaspore_xlsx, args.lifeform_xlsx)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2),
                    encoding="utf8")
    print(json.dumps({
        "version": result["version"],
        "source_reproduction_pass": result[
            "v19.9_broad_all_years_original_baseline_reproduced"],
        "label_cross_counts": result[
            "original_plots_with_actual_vascular_records_strict_x_broad_label_counts"],
        "results": {
            label: {
                era: {
                    "source_opportunity": v[
                        "original_label_and_survey_year_source_opportunity"],
                    "panels": {
                        name: {
                            "identified": panel["identified"],
                            "support": panel["source_support"],
                            "reasons": panel.get("reasons"),
                            "auc": panel.get(
                                "original_publication_equal_holdout_AUC"),
                            "study_paired": panel.get(
                                "descriptive_study_resampling_paired_AUC"),
                        } for name, panel in v[
                            "original_publication_holdout"].items()
                    },
                } for era, v in eras.items()
            } for label, eras in result["results"].items()
        },
    }, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
