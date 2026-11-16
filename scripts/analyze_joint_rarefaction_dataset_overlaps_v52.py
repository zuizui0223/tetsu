#!/usr/bin/env python3
"""v5.2: six-country fixed-K80 source cohort joint GBIF recording-bias diagnosis.

Combine TWO documented influences on source flora/physical rail correlation:
(1) records-per-cell rarefaction depth K=20,40,80
(2) all-Plantae GBIF DATASET_KEY source composition overlap
on exactly the SAME K>=80 source-eligible country cell pairs, 2010-17
and 2018-25, both predeclared Overture rail/road source anchors.

This is a descriptive, posthoc source robustness check, not evidence
of seed dispersal. No independent cell-pair p values. Source dataset
overlap is not survey-level effort correction. GBIF original row DOI
and verified active freight/rail vegetation remain unmeasured.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json,math
from pathlib import Path
from statistics import median

from analyze_6country_50km_network_plants_v44 import (
    BLOCKS_SIX,PERIODS,STRATEGIES,combine_source_regions,MIN_SOURCE_CELLS,MIN_VALID_PAIRS
)
from analyze_macro_50km_rail_road_aliens_v34 import partial_rank_r
from analyze_dataset_overlap_adjusted_rail_plants_v51 import (
    valid_dataset_profiles,dataset_overlap,merge_dataset_chunks
)
from analyze_rarefaction_depth_6country_v47 import (
    stable_probabilities,expected_overlap_ratio,network_source_pairs
)

DEPTHS=(20,40,80)
FIXED_RECORD_COHORT=80


def block_source(block,flora,rail,road,datasets):
    original,excluded=valid_dataset_profiles(block,flora,datasets)
    fixed={
        cid:obs for cid,obs in original.items()
        if all(sum(obs[p]["counts"].values())>=FIXED_RECORD_COHORT for p in PERIODS)
    }
    out={
        "country":flora["regions"][block]["country"],
        "n_cells_original_GRIIS_and_dataset_source_eligible":len(original),
        "n_cells_fixed_at_K80_both_periods":len(fixed),
        "excluded_dataset_quality_cells":len(excluded),
        "periods":{}
    }
    for period in PERIODS:
        out["periods"][period]={}
        for strategy in STRATEGIES:
            routes=network_source_pairs(block,rail,road,strategy)
            records=[]
            for a,b in combinations(sorted(fixed),2):
                route=routes.get((a,b))
                if route is None:continue
                left=fixed[a][period];right=fixed[b][period]
                n1=left["all_plant_indexed_records"]
                n2=right["all_plant_indexed_records"]
                if n1<=0 or n2<=0:raise ValueError("Source all-plants record effort invalid")
                d1=datasets["regions"][block]["cells"][a]["periods"][period]
                d2=datasets["regions"][block]["cells"][b]["periods"][period]
                records.append({
                    "a":a,"b":b,"network":route,
                    "baseline":[route["geo"],route["road"],
                                abs(math.log(n1/n2)),
                                (math.log(n1)+math.log(n2))/2],
                    "source_dataset_overlap":dataset_overlap(d1,d2),
                    "dataset_top1_fraction_gap":abs(
                        d1["top_datasetKey_source_share"]-d2["top_datasetKey_source_share"])
                })
            perdepth={}
            for k in DEPTHS:
                probs={cid:stable_probabilities(cell[period]["counts"],k)
                       for cid,cell in fixed.items()}
                yy=[];xx=[];zz=[];aa=[];bb=[]
                for r in records:
                    y=expected_overlap_ratio(probs[r["a"]],probs[r["b"]])
                    if y is None or not math.isfinite(y):continue
                    yy.append(y)
                    xx.append(r["network"]["X"])
                    zz.append(r["baseline"])
                    aa.append(r["source_dataset_overlap"])
                    bb.append(r["dataset_top1_fraction_gap"])
                coefficients={
                    "baseline_rho":partial_rank_r(xx,yy,zz) if xx else None,
                    "dataset_overlap_adjusted_rho":partial_rank_r(
                        xx,yy,[z+[p] for z,p in zip(zz,aa)]) if xx else None,
                    "dataset_overlap_plus_top1_share_gap_adjusted_rho":partial_rank_r(
                        xx,yy,[z+[p,q] for z,p,q in zip(zz,aa,bb)]) if xx else None
                }
                perdepth[str(k)]={
                    "same_fixed_K80_cell_cohort":True,
                    "n_qualified_50km_source_cells":len(fixed),
                    "n_nonindependent_source_cell_pairs":len(yy),
                    "coefficients":coefficients,
                    "no_railway_seed_transport_or_pairwise_p_value":True
                }
            out["periods"][period][strategy]=perdepth
        print(json.dumps({"block":block,"country":out["country"],"period":period,
            "fixed_K80_cells":len(fixed),
            "K20_baseline":out["periods"][period]["nearest_in_cell"]["20"]["coefficients"]["baseline_rho"],
            "K20_dataset_adjusted":out["periods"][period]["nearest_in_cell"]["20"]["coefficients"]["dataset_overlap_adjusted_rho"],
            "K80_baseline":out["periods"][period]["nearest_in_cell"]["80"]["coefficients"]["baseline_rho"],
            "K80_dataset_adjusted":out["periods"][period]["nearest_in_cell"]["80"]["coefficients"]["dataset_overlap_adjusted_rho"]
        },ensure_ascii=False),flush=True)
    return out


def aggregate(regions,period,strategy,k,field):
    values={};coverage={}
    for block,r in regions.items():
        row=r["periods"][period][strategy][str(k)]
        n=row["n_qualified_50km_source_cells"]; npair=row["n_nonindependent_source_cell_pairs"]
        rho=row["coefficients"][field]
        valid=(n>=MIN_SOURCE_CELLS and npair>=MIN_VALID_PAIRS and isinstance(rho,(float,int))
               and not isinstance(rho,bool) and math.isfinite(rho) and -1<=rho<=1)
        if valid:values[block]=float(rho)
        coverage[block]={"n_cells":n,"n_dependent_cell_pairs":npair,"supported":valid}
    vals=list(values.values())
    loo={b:round(median(v for c,v in values.items() if c!=b),6)
         for b in values if len(values)>1}
    return {
        "n_independent_country_blocks":len(vals),
        "n_negative_hypothesis_compatible_rail_directions":sum(v<0 for v in vals),
        "n_positive_opposite_directions":sum(v>0 for v in vals),
        "median_country_equal_partial_rank":round(median(vals),6) if vals else None,
        "country_partial_rhos":values,
        "leave_one_country_out_medians":loo,
        "source_support":coverage,
        "not_cell_pair_p_value_or_causal_train_seed_transfer":True
    }


def run(flora,rail,road,datasets):
    if not all(s.get("complete") is True for s in (flora,rail,road,datasets)):
        raise ValueError("One original source not complete")
    if flora.get("source_errors"):
        raise ValueError("Original GRIIS source error")
    for x in (rail,road,datasets):
        if set(x.get("regions",{}))!=set(BLOCKS_SIX):
            raise ValueError("Six-country source geography changed")
    if not set(BLOCKS_SIX).issubset(flora.get("regions",{})):
        raise ValueError("Original source GRIIS country frame missing")
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "study":"Fixed K80 source cell cohort, K20/40/80 indexed-record rarefaction and independent GBIF dataset recording composition overlap",
        "posthoc_noncausal_method_sensitivity":True,
        "country_blocks":list(BLOCKS_SIX),
        "n_frozen_country_blocks":6,
        "results":{},"equal_country_blocks":{}
    }
    for block in BLOCKS_SIX:
        out["results"][block]=block_source(block,flora,rail,road,datasets)
    for p in PERIODS:
        out["equal_country_blocks"][p]={}
        for strat in STRATEGIES:
            out["equal_country_blocks"][p][strat]={
                str(k):{
                    field:aggregate(out["results"],p,strat,k,field)
                    for field in (
                        "baseline_rho","dataset_overlap_adjusted_rho",
                        "dataset_overlap_plus_top1_share_gap_adjusted_rho"
                    )
                } for k in DEPTHS
            }
    out["limits"]=[
        "This is a POSTHOC RECORD and dataset-source sensitivity, not a pre-registered transport experiment.",
        "A GBIF DATASET_KEY is not an independent survey, nor a known observer or known detection probability.",
        "Dataset source overlaps are ALL Plantae indexed records, while response counts are GRIIS country-listed introduced plant keys.",
        "Same frozen 50km source pairs at K20/K40/K80 cannot eliminate time-variable recorder, climate, rail operation, or land cover confound.",
        "Only country-block summaries can be compared; pairwise source cell similarities share cells and have NO valid independent pair-level p-values.",
        "Overture physical 2026 rail/road source network cannot identify freight flows or actual railway-carried seed dispersal during GBIF observation years 2010–2025.",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--flora",required=True)
    p.add_argument("--rail-sources",nargs="+",required=True)
    p.add_argument("--road-sources",nargs="+",required=True)
    p.add_argument("--dataset-sources",nargs="+",required=True)
    p.add_argument("--output",default="results/audit/fixed80_dataset_source_rarefaction_6country_v52.json")
    a=p.parse_args()
    read=lambda paths:[json.loads(Path(q).read_text(encoding="utf-8")) for q in paths]
    flora=json.loads(Path(a.flora).read_text(encoding="utf-8"))
    rail=combine_source_regions(read(a.rail_sources))
    road=combine_source_regions(read(a.road_sources))
    ds=merge_dataset_chunks(read(a.dataset_sources))
    result=run(flora,rail,road,ds)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    for p in PERIODS:
        for k in DEPTHS:
            val=result["equal_country_blocks"][p]["nearest_in_cell"][str(k)]
            print(json.dumps({
                "period":p,"fixed80_rarefaction_K":k,
                "baseline":val["baseline_rho"]["median_country_equal_partial_rank"],
                "dataset_overlap_adjusted":val["dataset_overlap_adjusted_rho"]["median_country_equal_partial_rank"],
                "support":val["dataset_overlap_adjusted_rho"]["n_independent_country_blocks"],
                "hypothesis_negative_directions":val["dataset_overlap_adjusted_rho"]["n_negative_hypothesis_compatible_rail_directions"],
            },ensure_ascii=False),flush=True)
    print(json.dumps({"complete":True,"n_blocks":len(result["results"])}),flush=True)


if __name__=="__main__":
    main()
