#!/usr/bin/env python3
"""v5.7 actual leave-largest-GBIF-dataset-out GRIIS flora / railway paths.

Run a TRUE indexed source control, stronger than merely adjusting the
all-plants datasetKey compositional-overlap covariate: for EACH original
50km source cell×period, source v5.6 queries the full GRIIS-listed
introduced speciesKey RECORD counts, the current MOST ABUNDANT all-
Plantae datasetKey-only counts, and the complement (all other datasets).

This analysis compares FULL vs WITHOUT_TOP on the EXACT same cells
and EXACT same RAIL/ROAD mapped source pair universe at each K.
K=20,40,80 RECORD rarefaction for BOTH periods; cells require >=5
introduced indexed keys and >=K source records in both current source
FULL and WITHOUT_TOP. A K-specific cohort may shrink with K;
never compare coefficients across depths as identical populations.

Source covariates: physical geographic straight line, ROAD shortest
path, original GBIF all-Plantae record-count scale/balance, and current
all-Plantae GBIF datasetKey source overlap. Both anchor strategies
reported. No pairwise p-values; independent geography is six selected
country-blocks at most. Dataset UUIDs are not field surveys, historical
train freight service or actual propagule transfers.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json,math
from pathlib import Path
from statistics import median

from analyze_6country_50km_network_plants_v44 import (
    BLOCKS_SIX,PERIODS,STRATEGIES,combine_source_regions,
    MIN_SOURCE_CELLS,MIN_VALID_PAIRS
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_macro_50km_rail_road_aliens_v34 import _routes,partial_rank_r
from analyze_dataset_overlap_adjusted_rail_plants_v51 import (
    merge_dataset_chunks,dataset_overlap
)
from analyze_rarefaction_depth_6country_v47 import (
    stable_probabilities,expected_overlap_ratio
)

DEPTHS=(20,40,80)
COMPARISONS=("full","without_dominant_dataset")
METRIC_KEYS={
    "full":"full_source_current_country_GRIIS_aliens",
    "without_dominant_dataset":"remaining_datasets_country_GRIIS_aliens",
}
NO_P_VALUE=True


def merge_species_sources(chunks):
    out={}
    if not chunks:
        raise ValueError("Missing real country-source GRIIS top-dataset facets")
    for r in chunks:
        if r.get("complete") is not True or r.get("source_errors"):
            raise ValueError("Source GRIIS species partition incomplete")
        block=r.get("source_geo_block")
        if block not in BLOCKS_SIX or block in out:
            raise ValueError("Unexpected/duplicate source country block")
        if r.get("source_country")!=BLOCKS[block][0]:
            raise ValueError("Griis dataset partition source country mismatch")
        cells=r.get("source_cells")
        if not isinstance(cells,dict) or len(cells)!=16:
            raise ValueError("Missing complete 50km source cells")
        out[block]=r
    if set(out)!=set(BLOCKS_SIX):
        raise ValueError("Source did not cover exactly six pre-frozen country regions")
    return {"complete":True,"regions":out}


def records_at_k(row,metric,k):
    original=row[METRIC_KEYS[metric]]
    if not isinstance(original,dict) or any(not isinstance(key,str) or
        isinstance(n,bool) or not isinstance(n,int) or n<=0
        for key,n in original.items()):
        raise ValueError("Current indexed GBIF source species-key counts malformed")
    return len(original)>=5 and sum(original.values())>=k


def source_eligible_cell_pairs(block,flora,splits,compo,k):
    name=BLOCKS[block][0]
    original=flora["regions"][block]
    if original["country"]!=name or splits["regions"][block]["source_country"]!=name:
        raise ValueError("Source country mismatch in GBIF GRIIS leaveout")
    frag=splits["regions"][block]["source_cells"]
    data=compo["regions"][block]["cells"]
    if not(set(frag)==set(data)==set(original["cells"])):
        raise ValueError("Source 50km geographic cell identities do not match")
    selected={}
    allK={}
    for cid,old in original["cells"].items():
        if not old.get("source_pre_outcome_both_period_recording_eligible"):
            continue
        if not all(old["periods"].get(p,{}).get("alien_observed_minimum_5_species_and20_records") is True
                   for p in PERIODS):
            # Original source used 'alien_observed_minimum_5_species_and_20_records';
            # protect that exact historical key rather than silently
            # turning all cell sources into missing.
            if not all(old["periods"].get(p,{}).get(
                "alien_observed_minimum_5_species_and_20_records") is True for p in PERIODS):
                continue
        part=frag[cid]
        if part.get("country")!=name or part.get("original_two_period_country_GRIIS_eligible") is not True:
            raise ValueError("Current indexed source expected country-GRIIS original eligibility")
        comp=data[cid]
        if comp.get("country")!=name or not comp.get("original_GRIIS_two_period_source_eligible"):
            raise ValueError("Original datasetKey source eligibility mismatch")
        if not all(comp["periods"].get(p,{}).get("source_profile_valid_for_recording_bias_covariates") for p in PERIODS):
            continue
        if not all(part.get("periods",{}).get(p,{}).get(
            "dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys")
            for p in PERIODS):
            raise ValueError("GBIF dataset-specific source partition missing")
        good_all=all(records_at_k(part["periods"][p],"full",k) for p in PERIODS)
        good_other=all(records_at_k(part["periods"][p],"without_dominant_dataset",k)
                       for p in PERIODS)
        if good_all:allK[cid]=True
        if good_all and good_other:
            selected[cid]=part["periods"]
    return selected,len(allK)


def compute_block(block,flora,rail,road,splits,compo):
    result={"country":BLOCKS[block][0],"depths":{}}
    for k in DEPTHS:
        cells,n_allK=source_eligible_cell_pairs(block,flora,splits,compo,k)
        by_period={}
        for period in PERIODS:
            by_anchor={}
            profile={
                t:{cid:stable_probabilities(cells[cid][period][METRIC_KEYS[t]],k)
                   for cid in cells}
                for t in COMPARISONS
            }
            for strategy in STRATEGIES:
                rs=_routes(rail,block,strategy,"pair_routes")
                rd=_routes(road,block,strategy,"pair_paths")
                if set(rs)!=set(rd):raise ValueError("Rail/road source geometric pair universe mismatch")
                pairs=[]
                for a,b in combinations(sorted(cells),2):
                    r=rs[(a,b)];d=rd[(a,b)]
                    railkm=r.get("cell_centroid_snap_augmented_rail_path_km")
                    roadkm=d.get("centroid_snap_augmented_road_path_km")
                    geo=r.get("source_cell_centroid_geodesic_km")
                    rratio=r.get("rail_to_geodesic_ratio")
                    dratio=d.get("road_to_geodesic_ratio")
                    if not all(isinstance(x,(float,int)) and not isinstance(x,bool) and
                               math.isfinite(x) and x>0 for x in
                               (railkm,roadkm,geo,rratio,dratio)):
                        continue
                    if rratio<.95 or dratio<.95:continue
                    olda=flora["regions"][block]["cells"][a]["periods"][period]
                    oldb=flora["regions"][block]["cells"][b]["periods"][period]
                    n1=olda["indexed_all_plantae_source_records"]
                    n2=oldb["indexed_all_plantae_source_records"]
                    if n1<=0 or n2<=0:raise ValueError("Plant observation opportunity source denominator zero")
                    da=compo["regions"][block]["cells"][a]["periods"][period]
                    db=compo["regions"][block]["cells"][b]["periods"][period]
                    y={
                        name:expected_overlap_ratio(profile[name][a],profile[name][b])
                        for name in COMPARISONS
                    }
                    if any(v is None or not math.isfinite(v) for v in y.values()):
                        continue
                    pairs.append({
                        "x":math.log(railkm/roadkm),
                        "source_covariates":[math.log(geo),math.log(roadkm),
                            abs(math.log(n1/n2)),(math.log(n1)+math.log(n2))/2,
                            dataset_overlap(da,db)],
                        "responses":y,
                    })
                values={}
                for metric in COMPARISONS:
                    x=[p["x"] for p in pairs]
                    y=[p["responses"][metric] for p in pairs]
                    z=[p["source_covariates"] for p in pairs]
                    values[metric]=partial_rank_r(x,y,z) if pairs else None
                summary={
                    "n_source_cells_with_full_and_without_top_dataset_K_records_in_BOTH_periods":len(cells),
                    "n_source_cells_with_K_records_all_datasets_before_top_source_exclusion":n_allK,
                    "n_matched_RAIL_ROAD_independently_sourced_pairs_NOT_independent":len(pairs),
                    "partial_rank_same_pair_flora_full_and_after_dominant_dataset_removal":values,
                    "NO_pvalue_or_causal_seed_transport":True,
                }
                by_anchor[strategy]=summary
                if strategy=="nearest_in_cell":
                    print(json.dumps({
                        "block":block,"country":BLOCKS[block][0],
                        "period":period,"K":k,
                        "source_cells_all":n_allK,
                        "source_cells_after_dominant_dataset_removed":len(cells),
                        "matched_nonindependent_pairs":len(pairs),
                        "rho_all_datasets":values["full"],
                        "rho_without_dominant_dataset":values["without_dominant_dataset"],
                    },ensure_ascii=False),flush=True)
            by_period[period]=by_anchor
        result["depths"][str(k)]=by_period
    return result


def summarize_country(blocks,k,period,strategy,metric):
    vals={};audit={}
    for block,row in blocks.items():
        d=row["depths"][str(k)][period][strategy]
        n=d["n_source_cells_with_full_and_without_top_dataset_K_records_in_BOTH_periods"]
        pairs=d["n_matched_RAIL_ROAD_independently_sourced_pairs_NOT_independent"]
        value=d["partial_rank_same_pair_flora_full_and_after_dominant_dataset_removal"][metric]
        qualifies=n>=MIN_SOURCE_CELLS and pairs>=MIN_VALID_PAIRS and value is not None
        audit[block]={"n_source_50km_cells":n,
                      "n_dependent_cell_pair_rows":pairs,
                      "qualifies":qualifies}
        if qualifies:vals[block]=float(value)
    ordered=list(vals.values())
    loo={block:round(median(v for other,v in vals.items() if other!=block),6)
         for block in vals if len(vals)>1}
    return {
        "n_source_supported_independent_country_blocks":len(vals),
        "hypothesis_negative_country_directions":sum(v<0 for v in ordered),
        "opposite_positive_country_directions":sum(v>0 for v in ordered),
        "country_equal_median_partial_rank":round(median(ordered),6) if ordered else None,
        "regions":vals,"region_source_support":audit,
        "leave_one_country_out_medians":loo,
        "NO_pairwise_significance_or_causal_seed_claim":True,
    }


def run(flora,rail,road,source,ds):
    if not all(v.get("complete") is True for v in (flora,rail,road,source,ds)):
        raise ValueError("GBIF GRIIS/RAIL/ROAD source acquisition incomplete")
    if flora.get("source_errors"):raise ValueError("Original v3.6 GRIIS source errors")
    for v in (rail,road,source,ds):
        if set(v.get("regions",{}))!=set(BLOCKS_SIX):
            raise ValueError("Original six-country geographic source blocks not preserved")
    out={"utc":datetime.now(timezone.utc).isoformat(),
         "question":"Does observed source GRIIS-introduced plant similarity association persist after REMOVING the dominant GBIF all-plant source dataset in each cell-period?",
         "source_data_subtraction_is_not_observer_effort_or_true_plant_absence":True,
         "source_top_datasets_selected_from_all_plants_NOT_from_effect_signs":True,
         "regions":{},"equal_country_summaries":{}}
    for block in BLOCKS_SIX:
        out["regions"][block]=compute_block(block,flora,rail,road,source,ds)
    for k in DEPTHS:
        out["equal_country_summaries"][str(k)]={}
        for p in PERIODS:
            out["equal_country_summaries"][str(k)][p]={}
            for strategy in STRATEGIES:
                out["equal_country_summaries"][str(k)][p][strategy]={
                    metric:summarize_country(out["regions"],k,p,strategy,metric)
                    for metric in COMPARISONS
                }
    out["limits"]=[
        "Both FULL and WITHOUT_TOP outcomes are computed on EXACTLY the same K-specific source cells and matched physical RAIL/ROAD route pairs, preventing a cohort difference within each comparison.",
        "Top dataset from all-Plantae marginal source; not established to represent the same top INTRODUCED-plant source in every country.",
        "Current indexed speciesKey source facets are unfiltered status, with marginal DETECTED/NOT_DETECTED source quality gates, not individual occurrence statuses.",
        "GBIF speciesKey indexed records are not unique plants/field plots or survey absences; multiple dataset UUIDs do not imply independent confirmation.",
        "K20, K40 and K80 use different eligible cell populations; do not compare depth-specific country summaries as a monotonic ecological change.",
        "One country regional block represents one geographic result, not 120 shared-cell pair observations; no independent pairwise p-values.",
        "Physical 2026 Overture networks do not establish active historical freight transfer or rail-side floral trait selection; no causal train seed movement inferable.",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--flora",required=True)
    p.add_argument("--rail-sources",nargs="+",required=True)
    p.add_argument("--road-sources",nargs="+",required=True)
    p.add_argument("--species-sources",nargs="+",required=True)
    p.add_argument("--dataset-sources",nargs="+",required=True)
    p.add_argument("--output",default="results/audit/GBIF_sixcountry_topdataset_leaveout_flora_rail_v57.json")
    a=p.parse_args()
    read=lambda paths:[json.loads(Path(x).read_text(encoding="utf-8")) for x in paths]
    flora=json.loads(Path(a.flora).read_text(encoding="utf-8"))
    rail=combine_source_regions(read(a.rail_sources))
    road=combine_source_regions(read(a.road_sources))
    source=merge_species_sources(read(a.species_sources))
    ds=merge_dataset_chunks(read(a.dataset_sources))
    result=run(flora,rail,road,source,ds)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    for k in DEPTHS:
        s=result["equal_country_summaries"][str(k)]["2018_2025"]["nearest_in_cell"]
        print(json.dumps({"recent_K":k,
            "full_source_on_strict_same_pairs":s["full"],
            "WITHOUT_dominant_source_on_same_pairs":s["without_dominant_dataset"],
        },ensure_ascii=False),flush=True)
    print(json.dumps({"complete":True,"country_blocks":len(result["regions"]),
                      "NO_causal_train_dispersal":True},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
