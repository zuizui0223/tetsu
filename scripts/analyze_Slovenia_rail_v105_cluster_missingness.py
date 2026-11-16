#!/usr/bin/env python3
"""v10.5 railway-site/segment sensitivity of the original Slovenian 2026 field data.

This is a descriptive secondary analysis, not a treatment effect or a
test of seed carriage, train slipstream, mowing, plasticity or selection.
Original source: Šabić, Grašič & Rovanšek (2026), DOI
10.1007/s44473-026-00155-z; EXACT provenance and same-side design
are inherited unchanged from the successfully verified v10.3/v10.4.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_wind_diaspore_vs_annual_v104 import (
    WIND_STRA, WIND_SENS, SOURCE_DOI, SOURCE_SHA256,
    extract_original_xlsx, validate_original_keys, pair_source,
    measure, location_level, run as run_v104
)

BOOTSTRAP_DRAWS = 6000
BOOTSTRAP_SEED = 20261114
MEASURES = {
    "strict_annual": "source_strict_annual_BS_minus_BT",
    "strict_therophyte": "source_strict_therophyte_BS_minus_BT",
    "strict_wind_diaspore": WIND_STRA + "_delta",
    "wind_plus_Lycopodium_spores": WIND_SENS + "_delta",
    "cover_original_units_not_percent": "source_cover_BS_minus_BT",
}


def percentile(xs, p):
    if not xs:
        return None
    x = sorted(xs)
    at = (len(x) - 1) * p
    lo, hi = int(math.floor(at)), int(math.ceil(at))
    return x[lo] * (hi - at) + x[hi] * (at - lo) if lo != hi else x[lo]


def segment_cluster_diagnostic(rows, key, draws=BOOTSTRAP_DRAWS,
                               seed=BOOTSTRAP_SEED):
    """Resample railway SECTIONS (not track sides or individual locations).

    Descriptive uncertainty conditional on the nonrandom Slovenian sample,
    NOT a confidence interval for all railways or a causal p-value.
    """
    by = defaultdict(list)
    for r in rows:
        v = r.get(key)
        if v is None:
            continue
        if not math.isfinite(float(v)):
            raise ValueError("Nonfinite field location measure")
        by[str(r["section_number"])].append(float(v))
    ids = sorted(by)
    if len(ids) < 2:
        return {"n_source_segments": len(ids), "not_identifiable": True}
    site_mean = mean([v for k in ids for v in by[k]])
    segment_mean = mean([mean(by[k]) for k in ids])
    leave_site = [mean([v for k in ids if k != omitted for v in by[k]])
                  for omitted in ids]
    leave_segment = [mean([mean(by[k]) for k in ids if k != omitted])
                     for omitted in ids]
    rng = random.Random(seed)
    site_reps, seg_reps = [], []
    for _ in range(draws):
        chosen = [rng.choice(ids) for _ in ids]
        site_reps.append(mean([v for k in chosen for v in by[k]]))
        seg_reps.append(mean([mean(by[k]) for k in chosen]))
    return {
        "n_original_field_locations": sum(len(by[k]) for k in ids),
        "n_source_segments": len(ids),
        "n_cluster_bootstrap_draws": draws,
        "cluster_unit": "original_source_railway_section",
        "equal_field_location_mean_BS_minus_BT": round(site_mean, 6),
        "equal_rail_section_mean_BS_minus_BT": round(segment_mean, 6),
        "equal_field_location_section_cluster_bootstrap95": [
            round(percentile(site_reps, p), 6) for p in (.025, .975)],
        "equal_section_cluster_bootstrap95": [
            round(percentile(seg_reps, p), 6) for p in (.025, .975)],
        "leave_one_section_out_equal_location_range": [
            round(min(leave_site), 6), round(max(leave_site), 6)],
        "leave_one_section_out_equal_section_range": [
            round(min(leave_segment), 6), round(max(leave_segment), 6)],
        "n_leave_one_section_out_positive_equal_location": sum(x > 0 for x in leave_site),
        "n_leave_one_section_out_negative_equal_location": sum(x < 0 for x in leave_site),
        "n_sections_with_negative_mean": sum(mean(by[k]) < 0 for k in ids),
        "n_sections_with_positive_mean": sum(mean(by[k]) > 0 for k in ids),
        "nonprobability_original_one_country_sample_cannot_support_global_CI": True,
    }


def pair_missingness_bounds(row):
    """Sharp extreme all-taxa wind bounds when every missing strategy is unknown."""
    result = {}
    for zone in ("BS", "BT"):
        known = row["n_original_taxa_known_strategy_" + zone]
        missing = row["n_original_taxa_missing_strategy_" + zone]
        if not isinstance(known, int) or not isinstance(missing, int) or known < 0 or missing < 0:
            raise ValueError("Invalid original categorical trait counts")
        total = known + missing
        if total <= 0:
            raise ValueError("Original habitat zone has no plant taxa")
        p = row["wind_and_diaspore_group_shares_BS_minus_BT"]
        # Numerator below is an original wind fraction among known-strategy
        # taxa. Its counts are separately recovered from original trait data,
        # NEVER treating unknown-strategy plants as known nonwind.
        result[zone] = {"known": known, "missing": missing,
                        "coverage": known / total}
    return result


def paired_fraction_extreme_bounds(known_wind_BS, known_BS, missing_BS,
                                   known_wind_BT, known_BT, missing_BT):
    """BS-BT sharp bounds across unknown strategy status at the two sides."""
    for x in (known_wind_BS, known_BS, missing_BS,
              known_wind_BT, known_BT, missing_BT):
        if not isinstance(x, int) or x < 0:
            raise ValueError("Counts must be nonnegative integers")
    if (known_wind_BS > known_BS or known_wind_BT > known_BT or
        known_BS + missing_BS == 0 or known_BT + missing_BT == 0):
        raise ValueError("No valid source trait denominator")
    bs_total = known_BS + missing_BS
    bt_total = known_BT + missing_BT
    lo = known_wind_BS / bs_total - (known_wind_BT + missing_BT) / bt_total
    hi = (known_wind_BS + missing_BS) / bs_total - known_wind_BT / bt_total
    return lo, hi


def original_pair_wind_count(row, zone):
    known = row["n_original_taxa_known_strategy_" + zone]
    # v10.4 retains the source-known category fraction, but only the
    # *difference*. Reconstruct an integer numerator from the ORIGINAL
    # species×zone set instead of assuming a value for missing categories.
    if "_wind_original_counts" not in row:
        raise ValueError("Original exact wind counts not attached")
    wind = row["_wind_original_counts"][zone]
    if not isinstance(wind, int) or wind < 0 or wind > known:
        raise ValueError("Wind numerator disagrees with known strategy count")
    return wind


def exact_pair_annotation(paired, species, traits):
    """Join original rows by their prevalidated location & true track side."""
    for r in paired:
        lid, side = r["location_id"], r["track_side"]
        n = {}
        for zone in ("BS", "BT"):
            codes = species[(lid, side + zone)]
            n[zone] = sum(traits[c].get("Diaspore.strategy")
                          in ("Epilobium", "Phragmites") for c in codes)
            if n[zone] > r["n_original_taxa_known_strategy_" + zone]:
                raise ValueError("Source-strategy wind numerator exceeds known taxa")
        r["_wind_original_counts"] = n
    return paired


def section_and_missingness_report(paired, locations):
    # Sites, not sides, carry final equal-site weight.
    by = defaultdict(list)
    threshold_results = {}
    for r in paired:
        m = pair_missingness_bounds(r)
        bs, bt = m["BS"], m["BT"]
        low, high = paired_fraction_extreme_bounds(
            original_pair_wind_count(r, "BS"), bs["known"], bs["missing"],
            original_pair_wind_count(r, "BT"), bt["known"], bt["missing"])
        by[r["location_id"]].append({
            "coverage_delta": bs["coverage"] - bt["coverage"],
            "min_strategy_coverage": min(bs["coverage"], bt["coverage"]),
            "wind_bound_low": low,
            "wind_bound_high": high,
            "strict_wind_delta": r["wind_and_diaspore_group_shares_BS_minus_BT"][WIND_STRA],
        })
    if len(by) != 89 or sum(len(v) for v in by.values()) != 173:
        raise ValueError("Source site/side inventory drifted from v10.4")
    site = []
    for r in locations:
        site_sides = by[r["location_id"]]
        site.append({
            "section_number": r["section_number"],
            "coverage_delta": mean(p["coverage_delta"] for p in site_sides),
            "wind_bound_low": mean(p["wind_bound_low"] for p in site_sides),
            "wind_bound_high": mean(p["wind_bound_high"] for p in site_sides),
        })
    lo = mean(r["wind_bound_low"] for r in site)
    hi = mean(r["wind_bound_high"] for r in site)
    for threshold in (0.0, 0.5, 0.75, 0.9, 1.0):
        chosen = {}
        for r in locations:
            vs = [x["strict_wind_delta"] for x in by[r["location_id"]]
                  if x["min_strategy_coverage"] >= threshold
                  and x["strict_wind_delta"] is not None]
            if vs:
                chosen[r["location_id"]] = {
                    "section_number": r["section_number"],
                    "strict_wind_delta": mean(vs)}
        cells = list(chosen.values())
        threshold_results[str(threshold)] = {
            "n_original_locations_with_qualifying_sides": len(cells),
            "n_rail_sections": len(set(x["section_number"] for x in cells)),
            "equal_location_mean_BS_minus_BT": round(
                mean(x["strict_wind_delta"] for x in cells), 6) if cells else None,
            "number_of_original_side_pairs_qualifying": sum(
                x["min_strategy_coverage"] >= threshold
                and x["strict_wind_delta"] is not None
                for items in by.values() for x in items),
            "threshold_selected_without_replacing_missing_with_nonwind": True,
        }
    return {
        "source_Diaspore_strategy_coverage_BS_minus_BT": segment_cluster_diagnostic(
            site, "coverage_delta"),
        "fraction_of_original_89_locations_measured": 1,
        "extreme_missingness_wind_BS_BT_mean_sharp_bounds_all_taxa": [
            round(lo, 6), round(hi, 6)],
        "mean_sign_is_identified_under_every_unknown_trait_assignment":
            "negative" if hi < 0 else ("positive" if lo > 0 else "not_identified"),
        "strategy_complete_case_sensitivity_original_sites": threshold_results,
        "unknown_trait_values_never_imputed_as_nonwind": True,
    }


def verify_against_locked_v104(original, previous103, previous104):
    if (previous104.get("complete") is not True or
        previous104.get("original_supplement_sha256") != SOURCE_SHA256 or
        previous104.get("original_source_DOI") != SOURCE_DOI):
        raise ValueError("v10.4 verified original-source receipt not supplied")
    actual = run_v104(original, previous103)
    if (actual["n_original_source_rail_field_locations"] != 89 or
        actual["n_original_BS_vs_BT_same_side_pairs"] != 173):
        raise ValueError("Unexpected original published field pairing")
    for measure in (WIND_STRA, WIND_SENS):
        got = actual["observed_within_location_BS_minus_BT_group_differences"][measure]
        ref = previous104["observed_within_location_BS_minus_BT_group_differences"][measure]
        for key in ("n_independent_actual_RAIL_field_locations",
                    "mean_BS_minus_BT_fraction", "n_locations_BT_larger",
                    "n_locations_BS_larger"):
            if got[key] != ref[key]:
                raise ValueError("v10.5 did not reproduce v10.4 source " + measure + "/" + key)
    return actual


def run(original, prev103, prev104):
    verify_against_locked_v104(original, prev103, prev104)
    geos, zones, species, traits, _, provenance = validate_original_keys(original)
    original_pairs = pair_source(geos, zones, species, traits)
    paired = exact_pair_annotation(measure(original_pairs, species, traits),
                                   species, traits)
    locations = location_level(paired)
    if len(locations) != 89 or len(set(x["section_number"] for x in locations)) != 15:
        raise ValueError("Original field location/railway segment source keys changed")
    n_per_segment = defaultdict(int)
    for r in locations:
        n_per_segment[str(r["section_number"])] += 1
    eligible = {k for k, n in n_per_segment.items() if n >= 3}
    if len(eligible) != 11:
        raise ValueError("Original 11 source railway segments with >=3 locations changed")
    report = {
        "complete": True, "version": "v10.5",
        "source_original_DOI": SOURCE_DOI,
        "original_supplement_sha256": SOURCE_SHA256,
        "n_source_original_field_locations": len(locations),
        "n_source_original_same_side_BS_BT_pairs": len(paired),
        "n_original_railway_sections": len(n_per_segment),
        "n_railway_sections_with_3plus_sites": len(eligible),
        "section_cluster_bootstrap_diagnostics": {},
        "original_trait_missingness_sensitivity": section_and_missingness_report(paired, locations),
        "limitation": "Original Slovenian single-country nonrandom cross-sectional 2022-23 railway plots; no train-carried seed, mowing, slipstream, direct within-species response, heritability, road controls or causal intervention. Section bootstrap describes only variation in these sampled segments.",
    }
    for label, key in MEASURES.items():
        report["section_cluster_bootstrap_diagnostics"][label] = {
            "all_15_source_sections": segment_cluster_diagnostic(locations, key),
            "only_11_sections_with_at_least_three_field_locations":
                segment_cluster_diagnostic(
                    [r for r in locations if str(r["section_number"]) in eligible],
                    key),
        }
    return report


def main():
    p = argparse.ArgumentParser()
    for name in ("original-xlsx", "source-audit102", "verified103", "verified104", "output"):
        p.add_argument("--" + name, required=True)
    args = p.parse_args()
    original = extract_original_xlsx(args.original_xlsx, args.source_audit102)
    prev103 = json.loads(Path(args.verified103).read_text(encoding="utf-8"))
    prev104 = json.loads(Path(args.verified104).read_text(encoding="utf-8"))
    result = run(original, prev103, prev104)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
