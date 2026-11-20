#!/usr/bin/env python3
"""v30.1: Does a stature-elevation slope trained on OTHER genuine mountain
roads predict vegetation change along a NEW road with one calibration site?

This is MIREN mountain ROADS, not railways. It borrows no response from any
held-out road beyond the single, elevation-chosen anchor site. No genetic
adaptation, train dispersal, or railway-specific effect is identified.
The already source-verified v30.0 pipeline runs unchanged and is checked
against its committed six-panel frozen summary BEFORE accepting v30.1.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from statistics import mean
import numpy as np
import evaluate_MIREN_native_elevation_heldout_stature_predictability_v300 as source

TECHNICAL_SEEDS = (301201, 301202)
DRAWS_PER_CHAIN = 199
ZONES = ("roadside", "interior", "roadside_minus_interior")
ANCHORS = ("lowest", "median")
MIN_ROADS_FOR_DESCRIPTIVE_REGION = 2
MIN_REGIONS_FOR_MACRO = 3

def _road_xy(road: dict, T: np.ndarray, zone: str) -> tuple[np.ndarray, np.ndarray]:
    """The original physical site CWM and measured elevations; never export."""
    x = np.asarray(road["elev_m"], dtype=float) / 1000.0
    W = np.asarray(road["matrices"][zone], dtype=float)
    T = np.asarray(T, dtype=float)
    if (len(x) < 6 or np.ptp(x) < .15 - 1e-12 or W.ndim != 2 or
        T.ndim != 2 or T.shape[1] != 2 or W.shape != (len(x), len(T)) or
        not np.isfinite(W).all() or not np.isfinite(T).all() or
        not np.isfinite(x).all() or
        not np.allclose(W.sum(axis=1), 1, atol=1e-10)):
        raise ValueError("v30.1 requires authentic complete six-site source Road")
    return x, W @ T[:, 0]

def _within_road_slope(x: np.ndarray, y: np.ndarray) -> float:
    xc = x - x.mean()
    xx = float(xc @ xc)
    if xx <= 1e-12:
        raise ValueError("No genuine altitude spread")
    return float(xc @ (y-y.mean()) / xx)

def _anchor_index(x: np.ndarray, choice: str) -> int:
    if choice == "lowest":
        return int(np.argmin(x))
    if choice == "median":
        return int(np.argmin(np.abs(x-np.median(x))))
    raise ValueError("Prespecified anchors are lowest or median altitude")

def road_out_calibrated_transfer(
    roads: list[dict], T: np.ndarray, anchor: str = "lowest"
) -> dict:
    """Leave every entire physical road out of slope TRAINING.

    Fit slope on OTHER roads only, giving every training road equal weight.
    One held-out road site, selected using elevations only, supplies an
    intercept anchor. Its other original sites are wholly unseen outcomes.
    Compare transferred slope to the SAME anchor-as-constant baseline.
    Each source road, not each site, contributes one evaluation weight.
    """
    if len(roads) < MIN_ROADS_FOR_DESCRIPTIVE_REGION:
        raise ValueError("At least 2 genuine original roads for out-road CV")
    if anchor not in ANCHORS:
        raise ValueError("Unexpected anchor rule")
    by_zone = {}
    for zone in ("roadside", "interior"):
        xy = [_road_xy(r, T, zone) for r in roads]
        slopes = [_within_road_slope(x, y) for x, y in xy]
        gains, positives = [], 0
        for i, (x, y) in enumerate(xy):
            # The excluded road's slope and all its target y are NOT used
            # to fit b. Calibration uses only the preselected single site.
            b_train = mean(s for j, s in enumerate(slopes) if j != i)
            a = _anchor_index(x, anchor)
            take = np.arange(len(x)) != a
            reference = float(y[a])
            pred0 = np.full(int(take.sum()), reference, dtype=float)
            pred1 = reference + b_train * (x[take]-x[a])
            rmse0 = float(np.sqrt(np.mean((y[take]-pred0)**2)))
            rmse1 = float(np.sqrt(np.mean((y[take]-pred1)**2)))
            gain = rmse0 - rmse1
            if not all(map(math.isfinite, (b_train, rmse0, rmse1, gain))):
                raise ValueError("Non-finite original road heldout prediction")
            gains.append(gain)
            positives += gain > 0
        by_zone[zone] = {
            "gain_zlog_height_RMSE": mean(gains),
            "positive_original_heldout_roads": positives,
            "original_heldout_road_count": len(roads),
            "one_original_site_per_road_used_only_for_intercept": True,
            "target_sites_not_used_in_training_or_calibration": True,
        }
    by_zone["roadside_minus_interior"] = {
        "gain_zlog_height_RMSE":
            by_zone["roadside"]["gain_zlog_height_RMSE"] -
            by_zone["interior"]["gain_zlog_height_RMSE"],
        "original_heldout_road_count": len(roads),
    }
    return by_zone

def conditional_trait_label_null(
    roads: list[dict], T: np.ndarray, groups: list[np.ndarray],
    tree_idx: int, region_idx: int, mode: str, draws: int = DRAWS_PER_CHAIN
) -> dict:
    if draws != DRAWS_PER_CHAIN:
        raise ValueError("Predeclared 2 chains x 199 shuffles required")
    result = {}
    for anchor in ANCHORS:
        obs = road_out_calibrated_transfer(roads, T, anchor)
        placebo = {z: [] for z in ZONES}
        for seed in TECHNICAL_SEEDS:
            rng = np.random.default_rng(
                seed + 17011*tree_idx + 5009*region_idx +
                sum(map(ord, mode+anchor)))
            for _ in range(draws):
                shuffled = source.permute_paired_traits(T, groups, rng)
                d = road_out_calibrated_transfer(roads, shuffled, anchor)
                for z in ZONES:
                    placebo[z].append(d[z]["gain_zlog_height_RMSE"])
        zone_effects = {}
        for z in ZONES:
            actual = obs[z]["gain_zlog_height_RMSE"]
            null = mean(placebo[z])
            zone_effects[z] = {
                "real_out_road_calibrated_slope_gain": actual,
                "phylo75_x_occupancy_joint_trait_label_null_gain": null,
                "real_minus_joint_trait_label_null_gain": actual-null,
                "technical_draws_NOT_ecological_field_replicates": len(placebo[z]),
                "positive_heldout_original_roads":
                    obs[z].get("positive_original_heldout_roads"),
                "true_heldout_original_roads":
                    obs[z]["original_heldout_road_count"],
            }
        result[anchor] = zone_effects
    return {
        "status": "SOURCE_COMPUTABLE_ROAD_OUT_WITH_ONE_FIXED_ANCHOR",
        "original_road_count": len(roads),
        "same_original_community_sites_traits_elevations_as_v300": True,
        "definition": "TRAIN slopes on all OTHER roads, intercept calibrates from exactly one elevation-chosen ORIGINAL site of heldout road; remaining sites genuinely held out",
        "anchor_effects": result,
    }

def region_equal_aggregate(regions: dict) -> dict:
    """Macro averaging is over eligible source regions, never field sites."""
    if len(regions) < MIN_REGIONS_FOR_MACRO:
        return {
            "status": "NOT_IDENTIFIABLE_THREE_SOURCE_REGIONS",
            "source_regions": len(regions),
            "effects": None,
        }
    out = {}
    for anchor in ANCHORS:
        out[anchor] = {}
        for zone in ZONES:
            raw = {r: z["anchor_effects"][anchor][zone][
                    "real_out_road_calibrated_slope_gain"]
                   for r, z in regions.items()}
            exc = {r: z["anchor_effects"][anchor][zone][
                    "real_minus_joint_trait_label_null_gain"]
                   for r, z in regions.items()}
            leave = [mean(v for r, v in exc.items() if r != omitted)
                     for omitted in exc]
            out[anchor][zone] = {
                "true_observed_region_equal_gain": mean(raw.values()),
                "region_equal_excess_beyond_joint_label_null": mean(exc.values()),
                "n_regions_positive_real_prediction": sum(v > 0 for v in raw.values()),
                "n_regions_positive_conditional_excess": sum(v > 0 for v in exc.values()),
                "source_regions": len(regions),
                "region_delete_one_excess_range_NOT_global_95CI":
                    [min(leave), max(leave)],
                "region_codes_and_excess_only": exc,
            }
    return {
        "status": "ELIGIBLE_DESCRIPTIVE_REGION_EQUAL_NOT_STUDY_META_ANALYSIS",
        "source_regions": len(regions),
        "effects": out,
    }

def verify_frozen_original_v300(out: dict, reference: dict) -> None:
    if (out["version"] != "v30.0" or reference["version"] != "v30.0" or
        not out["original_v299_six_tree_by_support_panel_true_native_altitude_effects_reproduced_to_1e_10"]):
        raise ValueError("Original v30.0/v29.9 source pipeline not reproduced")
    for frozen in reference["source_tree_by_field_support_panel"]:
        tree, kind = frozen["original_tree"], frozen["source_trait_support"]
        panel = ("primary_5species_25percent" if kind == "primary" else
                 "strict_10species_50percent")
        raw = out["region_equal_source_environmental_heldout_predictability"][
            tree][panel]["observed_original_site_heldout_elevation_gain"][
            "source_region_equal_effects"]
        if raw is None:
            raise ValueError("Original v30.0 supported region disappeared")
        old = {
            "roadside": (
                "real_observed_original_roadside_site_elevation_LOO_RMSE_gain",
                "original_source_region_equal_mean_heldout_Elevation_predictive_gain"),
            "roadside_excess": (
                "roadside_gain_excess_vs_joint_phylo75_occ_label_null",
                "original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null"),
            "interior_excess": (
                "inland_gain_excess_vs_joint_phylo75_occ_label_null",
                "original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null")
        }
        for label, (fk, rk) in old.items():
            zone = "interior" if label == "interior_excess" else "roadside"
            if label == "roadside_excess":
                key = rk
            else:
                key = rk
            if abs(float(raw[zone][key]) - float(frozen[fk])) > 1e-10:
                raise ValueError("Original first party v30.0 frozen outcome drift")
        if raw["roadside"]["source_original_region_codes_evaluated"] != frozen[
            "original_source_region_codes_with_true_heldout_prediction"]:
            raise ValueError("Original v30.0 region cohort drift")

def run(args: argparse.Namespace) -> dict:
    original_hook = source.source_trait_relabel_predictability_null
    def add_prespecified_road_transfer(
        roads, T, groups, tree_idx, region_idx, mode,
        draws=source.DRAWS_PER_CHAIN
    ):
        base = original_hook(roads, T, groups, tree_idx, region_idx, mode, draws)
        base["v301_road_out_calibrated_transfer"] = conditional_trait_label_null(
            roads, T, groups, tree_idx, region_idx, mode, draws)
        return base
    source.source_trait_relabel_predictability_null = add_prespecified_road_transfer
    try:
        original = source.run(
            args.miren, args.height, args.seed, args.tree_dir,
            args.original_v299_summary, args.output, args.draws)
    finally:
        source.source_trait_relabel_predictability_null = original_hook
    verify_frozen_original_v300(
        original, json.loads(args.original_v300_frozen_summary.read_text()))
    summary = {}
    for tree, data in original["original_first_party_tree_source_and_original_region_support"].items():
        summary[tree] = {}
        for panel, info in data["panels"].items():
            regions = {}
            stable = {}
            roads_total = 0
            sites_total = 0
            for region, item in info["original_source_region_support_and_predictivity"].items():
                base = item.get("new_predictive_heldout_effects")
                if base is None:
                    continue
                value = base["v301_road_out_calibrated_transfer"]
                regions[region] = value
                roads_total += value["original_road_count"]
                sites_total += item["source_CV_site_support"][
                    "original_heldout_sites_on_eligible_roads"]
                year = item.get("original_road_single_field_year_only_predictive_sensitivity")
                if year is not None:
                    stable[region] = year["v301_road_out_calibrated_transfer"]
            summary[tree][panel] = {
                "source_regions_with_at_least_two_original_roads": len(regions),
                "original_road_count": roads_total,
                "original_paired_field_sites_used": sites_total,
                "primary": region_equal_aggregate(regions),
                "single_original_survey_year_roads_sensitivity":
                    region_equal_aggregate(stable),
                "no_original_plot_species_GPS_or_elevations_exported": True,
            }
    result = {
        "version": "v30.1",
        "status": "SOURCE_RECALCULATED_AFTER_INDEPENDENT_v300_FROZEN_AUDIT",
        "independent_original_botanical_programmes_added": 0,
        "genuine_MIREN_mountain_ROADS_not_railway": True,
        "n_technical_label_assignments_NOT_field_p_values": 398,
        "original_MIREN_v300_six_panels_source_observed_gain_reproduced_to_1e_10": True,
        "one_heldout_road_anchor_height_USED_for_intercept_not_a_cold_start": True,
        "train_slope_only_on_other_physically_distinct_original_roads": True,
        "two_anchors_a_priori": list(ANCHORS),
        "source_tree_panel_region_aggregates": summary,
        "no_global_transfer_claim_and_no_CAUSAL_railway_claim": True,
    }
    args.v301_output.parent.mkdir(parents=True, exist_ok=True)
    args.v301_output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    print(json.dumps({
        "version": result["version"],
        "results": {t: {
            p: {
                "original_regions": x["source_regions_with_at_least_two_original_roads"],
                "original_roads": x["original_road_count"],
                "low_anchor_roadside_real_gain":
                    (x["primary"]["effects"] or {}).get("lowest", {}).get(
                        "roadside", {}).get("true_observed_region_equal_gain"),
                "low_anchor_roadside_excess":
                    (x["primary"]["effects"] or {}).get("lowest", {}).get(
                        "roadside", {}).get("region_equal_excess_beyond_joint_label_null"),
            } for p,x in panels.items()
        } for t,panels in summary.items()},
    }, ensure_ascii=False, indent=2))
    return result

def main():
    p = argparse.ArgumentParser()
    for name in ("miren", "height", "seed", "tree-dir", "original-v299-summary",
                 "original-v300-frozen-summary", "output", "v301-output"):
        p.add_argument("--"+name, required=True, type=Path)
    p.add_argument("--draws", type=int, default=DRAWS_PER_CHAIN)
    a = p.parse_args()
    if a.draws != DRAWS_PER_CHAIN:
        raise ValueError("Frozen v30.1 conditional null draws changed")
    run(a)

if __name__ == "__main__":
    main()
