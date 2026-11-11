#!/usr/bin/env python3
"""v7.3: actual France GBIF original 81-file GRIIS event source vs old
speciesKey indexed facets; GIFT growth form; frozen 50km RAIL/ROAD paths.

This is a SOURCE-definition diagnostic, not measured train dispersal,
freight activity, cutting/wind selection, plasticity or evolution.
The 120 within-block cell pairs are dependent, with no pairwise p-value.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from analyze_macro_50km_rail_road_aliens_v34 import _routes, partial_rank_r
from analyze_macro_rail_assemblages_v22 import (
    incidence_probabilities, expected_rarefied_overlap)

BLOCK = "fr_burgundy"
PERIOD = "2018_2025"
GROUPS = ("herbaceous", "woody")
SOURCES = ("original_event_all", "original_event_precise5km", "indexed_facet_2026")
STRATEGIES = ("nearest_in_cell", "nearest_largest_component_in_cell")
K = 20
MIN_GROUP_TAXA = 3
MIN_CELLS = 8
MIN_PAIRS = 20
EXPECTED_ALL_PLANT_BBOX_ROWS = 697577
EXPECTED_GRIIS_EVENT_COUNT = 24279
EXPECTED_GRIIS_5KM_COUNT = 18991


def col_growth_map(crosswalk, gift):
    if (crosswalk.get("countries", {}).get("FR") is None or
            gift.get("complete") is not True or
            gift.get("source", {}).get("GIFT_meta_ID") != "1.2.1"):
        raise ValueError("Official France GRIIS-COL / vetted GIFT source missing")
    fr = crosswalk["countries"]["FR"]
    old_to_col = fr["old_GRIIS_to_COL_strict_accepted_map"]
    old_to_growth = gift["growth_form_by_source_gbif_species_key"]
    by_col = defaultdict(set)
    for old, col in old_to_col.items():
        if old in old_to_growth:
            if old_to_growth[old] not in GROUPS:
                raise ValueError("Unvetted GIFT growth form")
            by_col[col].add(old_to_growth[old])
    return ({c: next(iter(groups)) for c, groups in by_col.items() if len(groups) == 1},
            {c: list(groups) for c, groups in by_col.items() if len(groups) > 1},
            old_to_col)


def assign_groups(counts, growth):
    result = {g: {} for g in GROUPS}
    unknown = 0
    if not isinstance(counts, dict):
        raise ValueError("Missing source species-key counts")
    for taxon, n in counts.items():
        if type(n) is not int or n <= 0:
            raise ValueError("Invalid original source record count")
        group = growth.get(taxon)
        if group is None:
            unknown += n
        else:
            result[group][taxon] = result[group].get(taxon, 0) + n
    return result, unknown


def build_cells(events, crosswalk, gift, flora):
    src = events.get("source_pre_GRIIS_taxon_outcome_81_file_selection", {})
    if (events.get("frozen_geographical_block") != BLOCK or
        events.get("source_period_analyzed") != PERIOD or
        events.get("n_original_equal_area_50km_cells") != 16 or
        events.get("NO_RAIL_effect_or_train_dispersal_result") is not True or
        events.get("source_archive_NOT_complete_regional_original_events_until_A5_partition_gap_resolved") is not True or
        src.get("n_all_81_bbox_rows_EXPECTED") != EXPECTED_ALL_PLANT_BBOX_ROWS or
        src.get("n_source_files_in_conservative_A5_geography_candidate_union") != 81 or
        events.get("n_country_GRIIS_strict_COL_event_species") != EXPECTED_GRIIS_EVENT_COUNT or
        events.get("n_country_GRIIS_strict_COL_event_species_known_uncertainty_le5000m") != EXPECTED_GRIIS_5KM_COUNT):
        raise ValueError("Actual frozen original France GBIF source provenance differs")
    if (flora.get("complete") is not True or flora.get("source_errors") or
        flora.get("regions", {}).get(BLOCK, {}).get("country") != "FR"):
        raise ValueError("Original France GRIIS facet source incomplete")
    growth, ambiguous, old_to_col = col_growth_map(crosswalk, gift)
    native_cells = events["by_period"][PERIOD]["cells"]
    indexed_cells = flora["regions"][BLOCK]["cells"]
    if set(native_cells) != set(indexed_cells) or len(native_cells) != 16:
        raise ValueError("Original event/legacy indexed cell identifiers do not match")
    cells = {}
    for cid, ev in native_cells.items():
        old = indexed_cells[cid]
        if not old.get("source_pre_outcome_both_period_recording_eligible"):
            raise ValueError("Original same 50km cell not pre-outcome GBIF eligible")
        old_time = old.get("periods", {}).get(PERIOD, {})
        if not old_time.get("full_facet_source_pagination_completed"):
            raise ValueError("Original indexed speciesKey source not completely paginated")
        legacy_counts = old_time["observed_alien_specieskey_counts"]
        harmonized = defaultdict(int)
        for legacy_key, n in legacy_counts.items():
            if legacy_key in old_to_col:
                harmonized[old_to_col[legacy_key]] += n
        variants = {
            "original_event_all": ev["COL_strict_country_GRIIS_speciesKey_event_counts"],
            "original_event_precise5km": ev["COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts"],
            "indexed_facet_2026": harmonized,
        }
        cells[cid] = {"variants": {}, "unknown": {},
                      "indexed_all_plantae_source_records":old_time["indexed_all_plantae_source_records"]}
        for key, counts in variants.items():
            groups, unknown = assign_groups(counts, growth)
            cells[cid]["variants"][key] = groups
            cells[cid]["unknown"][key] = unknown
        if sum(map(sum, (v.values() for v in cells[cid]["variants"]["original_event_all"].values()))) + cells[cid]["unknown"]["original_event_all"] != ev["n_GRIIS_strict_COL_event_species"]:
            raise ValueError("Original all-event strict COL record conservation failed")
        if sum(map(sum, (v.values() for v in cells[cid]["variants"]["original_event_precise5km"].values()))) + cells[cid]["unknown"]["original_event_precise5km"] != ev["n_GRIIS_strict_COL_event_species_known_uncertainty_le5000m"]:
            raise ValueError("Original ≤5km source record conservation failed")
    return cells, len(growth), len(ambiguous)


def source_summaries(cells):
    summaries = {}
    for source in SOURCES:
        ns = {g: 0 for g in GROUPS}
        species = {g: set() for g in GROUPS}
        unknown = 0
        for v in cells.values():
            unknown += v["unknown"][source]
            for group in GROUPS:
                x = v["variants"][source][group]
                ns[group] += sum(x.values())
                species[group].update(x)
        summaries[source] = {
            "count_granularity": "unique original event-species source" if source != "indexed_facet_2026" else "indexed facet speciesKey-record source",
            "record_counts_by_growth_form": {**ns, "unknown": unknown},
            "distinct_COL_species_categories_by_growth_form":{g: len(species[g]) for g in GROUPS},
            "percentage_classified": round(100 * (ns["herbaceous"]+ns["woody"]) /
                                           (sum(ns.values())+unknown), 3)
        }
    return summaries


def compare_anchor(cells, rail, road, anchor):
    railpairs = _routes(rail, BLOCK, anchor, "pair_routes")
    roadpairs = _routes(road, BLOCK, anchor, "pair_paths")
    if set(railpairs) != set(roadpairs):
        raise ValueError("Physical rail/road pair universes differ")
    rows = []
    for a, b in combinations(sorted(cells), 2):
        r, d = railpairs[(a, b)], roadpairs[(a, b)]
        ra, rd = r.get("cell_centroid_snap_augmented_rail_path_km"), d.get("centroid_snap_augmented_road_path_km")
        geo = r.get("source_cell_centroid_geodesic_km")
        ratios = [r.get("rail_to_geodesic_ratio"), d.get("road_to_geodesic_ratio")]
        if any(type(t) not in (int, float) or not math.isfinite(t) or t <= 0
               for t in [ra, rd, geo, *ratios]):
            continue
        if min(ratios) < .95:
            continue
        n1 = cells[a]["indexed_all_plantae_source_records"]
        n2 = cells[b]["indexed_all_plantae_source_records"]
        if type(n1) is not int or type(n2) is not int or min(n1,n2) <= 0:
            raise ValueError("Missing indexed source observation effort")
        y = {}
        valid = True
        for source in SOURCES:
            y[source] = {}
            for group in GROUPS:
                aa, bb = (cells[c]["variants"][source][group] for c in (a,b))
                if min(len(aa),len(bb)) < MIN_GROUP_TAXA or min(sum(aa.values()),sum(bb.values())) < K:
                    valid = False
                    continue
                pa = incidence_probabilities(aa, k=K)
                pb = incidence_probabilities(bb, k=K)
                similar = expected_rarefied_overlap(pa, pb)
                if similar is None:
                    valid = False
                    continue
                raw = len(set(aa) & set(bb)) / len(set(aa) | set(bb))
                y[source][group] = {"rarefied20": similar, "raw_jaccard": raw}
        if valid:
            rows.append({
                "a":a,"b":b,
                "x":math.log(ra/rd),
                "covars":[math.log(geo),math.log(rd),abs(math.log(n1/n2)),(math.log(n1)+math.log(n2))/2],
                "responses":y
            })
    n_cells = len({c for p in rows for c in (p["a"], p["b"])})
    report = {"n_same_pair_cohort":len(rows), "n_distinct_source_50km_cells":n_cells,
              "independent_replicates_are_region_blocks_not_cell_pairs":True,
              "responses":{}}
    for source in SOURCES:
        report["responses"][source] = {}
        for group in GROUPS:
            report["responses"][source][group] = {}
            for metric in ("rarefied20","raw_jaccard"):
                scores = [p["responses"][source][group][metric] for p in rows]
                r = (partial_rank_r([p["x"] for p in rows], scores,
                                    [p["covars"] for p in rows])
                     if n_cells >= MIN_CELLS and len(rows) >= MIN_PAIRS else None)
                report["responses"][source][group][metric] = {"partial_rho":r,
                                "n_dependent_cell_pairs":len(rows), "n_distinct_cells":n_cells}
    return report


def run(events, crosswalk, gift, flora, rail, road):
    cells,nmapped,nambig = build_cells(events, crosswalk, gift, flora)
    result = {
        "complete": True, "source_geography": BLOCK, "period": PERIOD,
        "country_independent_macro_replicate_count": 1,
        "data_sources": "actual 81-source-object GBIF Parquet GRIIS COL original events vs independently archived old GBIF facets; GIFT v7.0 and 2026 Overture",
        "n_COL_species_growth_form_mapped": nmapped, "n_ambiguous_growth_forms_excluded":nambig,
        "source_summaries":source_summaries(cells),
        "anchor_sensitivity":{s:compare_anchor(cells,rail,road,s) for s in STRATEGIES},
        "causal_seed_transport_or_measured_stress_evolution":False,
        "limitations":[
            "This is one geographical 200 km block; 120 shared-cell comparisons are NOT independent replicates, no inferential pairwise p-values.",
            "81 development Parquet objects represent a bounded nonexhaustive source; A5 partition polygon/centroid mismatch remains unresolved.",
            "Indexed all-Plantae facet effort covariates have different provenance from original Parquet event records and do not fully control observer/dataset effects.",
            "Mapped rail/road 2026 physical routes are not observed active historical freight or seed flux.",
            "GIFT species growth form does not measure wind/mowing exposure, regeneration capacity or selection; no plasticity or local adaptation is identified.",
            "Changing similarity metric from record-rarefied to raw Jaccard may change ecological association signs."]}
    return result


def main():
    ap = argparse.ArgumentParser()
    for key in ("events","crosswalk","gift","flora","rail","road","output"):
        ap.add_argument("--"+key,required=True)
    args = ap.parse_args()
    load = lambda k: json.loads(Path(getattr(args,k)).read_text(encoding="utf-8"))
    out = run(*(load(k) for k in ("events","crosswalk","gift","flora","rail","road")))
    dst=Path(args.output); dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete": True,
        "original_trait_coverage":out["source_summaries"]["original_event_all"],
        "primary_nearest_in_cell":out["anchor_sensitivity"]["nearest_in_cell"],
        "no_causal_railway_dispersal_or_stress_adaptation":True},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
