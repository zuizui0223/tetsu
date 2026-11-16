#!/usr/bin/env python3
"""v5.1 REAL six-country data-source-adjusted observed flora vs RAIL/ROAD routes.

Join original source-frozen 50km 6-country GBIF speciesKey GRIIS-listed
introduced flora, original 2026 Overture RAIL and ROAD source routes and
new live GBIF FULL DATASET_KEY facet fractions per same cell×period.

Compare IDENTICAL pairwise geographic observations:
  baseline partial Spearman-style rank residualization as in v4.4
  baseline plus full indexed all-plant DATASET composition overlap
  baseline plus overlap AND per-cell max dataset dominance imbalance.

A dataset overlap metric is NOT proof of the same observer or detection
protocol or a causal correction. Dataset composition is ALL PLANT
records (not only country-GRIIS introduced taxa), while the outcome
is GRIIS-filtered indexed positive speciesKey record similarity.
Dataset source statuses are marginal DETECTED/NOT_DETECTED counts, and
the current GBIF indexed record source is not DOI-backed raw occurrences.

NO pairwise pseudo-independent p-values. One source-frozen 200km
country geographic block is the independent unit, not 120 pairwise rows.
All 6 country blocks, two time periods, two anchor strategies and
raw/Jaccard + 20-record rarefied outcomes must be reported.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json,math
from pathlib import Path
from statistics import median

from analyze_macro_50km_rail_road_aliens_v34 import (
    _eligible_both_period_cells,_routes,partial_rank_r,STRATEGIES
)
from analyze_6country_50km_network_plants_v44 import (
    BLOCKS_SIX, PERIODS, combine_source_regions,
    MIN_SOURCE_CELLS, MIN_VALID_PAIRS
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_macro_rail_assemblages_v22 import expected_rarefied_overlap

RESPONSES=("rarefied20_indexed_source_overlap","raw_indexed_source_species_Jaccard")


def merge_dataset_chunks(chunks):
    if not chunks or not isinstance(chunks,list):
        raise ValueError("Actual per-country GBIF dataset source missing")
    joined={}
    for chunk in chunks:
        if chunk.get("complete") is not True or chunk.get("errors"):
            raise ValueError("GBIF dataset source acquisition incomplete or errors")
        block=chunk.get("block")
        if block not in BLOCKS_SIX or block in joined:
            raise ValueError("Unpredeclared or duplicate source country block")
        if chunk.get("country")!=BLOCKS[block][0]:
            raise ValueError("GBIF dataset source country identity mismatch")
        cells=chunk.get("cells")
        if not isinstance(cells,dict) or len(cells)!=16:
            raise ValueError("Incomplete 16-cell original source geographic frame")
        joined[block]={"country":chunk["country"],"cells":cells,
                       "recording_source":chunk["source"]}
    if set(joined)!=set(BLOCKS_SIX):
        raise ValueError("Six predeclared country dataset sources not all retrieved")
    return {"complete":True,"regions":joined,
            "n_unique_country_blocks_with_true_source_dataset_composition":6}


def valid_dataset_profiles(block,flora,datasets):
    if datasets["regions"][block]["country"]!=BLOCKS[block][0]:
        raise ValueError("Actual GBIF dataset source country mismatches original GRIIS flora")
    plants=_eligible_both_period_cells(block,flora)
    raw=datasets["regions"][block]["cells"]
    if set(raw)!=set(flora["regions"][block]["cells"]):
        raise ValueError("Dataset and original GBIF GRIIS source cell geometry mismatch")
    included={}
    rejections={}
    for cid,source_flora in plants.items():
        data=raw.get(cid)
        if not data or data.get("country")!=BLOCKS[block][0] or not data.get("original_GRIIS_two_period_source_eligible"):
            raise ValueError("Actual GBIF datasets no longer match original source-eligible cell")
        okay=True
        for p in PERIODS:
            entry=data.get("periods",{}).get(p,{})
            if not entry.get("source_profile_valid_for_recording_bias_covariates"):
                okay=False
                break
            profile=entry.get("dataset_key_record_counts")
            total=entry.get("source_indexed_all_Plantae_records")
            if not isinstance(profile,dict) or not isinstance(total,int) or total<=0:
                raise ValueError("Missing actual complete indexed dataset source")
            if any(not isinstance(k,str) or not k or isinstance(v,bool) or not isinstance(v,int) or v<0
                   for k,v in profile.items()) or sum(profile.values())>total:
                raise ValueError("Malformed GBIF dataset proportions or record counts")
        if okay:
            included[cid]=source_flora
        else:
            rejections[cid]="GBIF dataset marginal unknown, non-detected >2%, or source drift >10%"
    return included,rejections


def dataset_overlap(p,q):
    """Lower bound on weighted shared dataset source fraction.

    p,q = complete datasetKey facets / full all-Plantae source count.
    Unattributed records do not become an artificial shared dataset.
    """
    a=p["dataset_key_record_counts"];b=q["dataset_key_record_counts"]
    n1=p["source_indexed_all_Plantae_records"]
    n2=q["source_indexed_all_Plantae_records"]
    if n1<=0 or n2<=0:
        raise ValueError("Zero indexed source dataset denominator")
    fraction=math.fsum(min(na/n1,b[k]/n2) for k,na in a.items() if k in b)
    if not math.isfinite(fraction) or fraction < -1e-12 or fraction>1+1e-8:
        raise ValueError("Unphysical shared source dataset composition fraction")
    return min(1.,max(0.,fraction))


def compute_country(block,flora,rail,road,datasets):
    included,excluded=valid_dataset_profiles(block,flora,datasets)
    outcome={
        "country":BLOCKS[block][0],
        "n_original_qualified_two_period_flora_cells":len(_eligible_both_period_cells(block,flora)),
        "n_true_dataset_quality_qualified_cells":len(included),
        "source_cells_excluded_due_to_dataset_marginal_quality":len(excluded),
        "periods":{}
    }
    for period in PERIODS:
        outcome["periods"][period]={}
        for strategy in STRATEGIES:
            rsource=_routes(rail,block,strategy,"pair_routes")
            dsource=_routes(road,block,strategy,"pair_paths")
            if set(rsource)!=set(dsource):
                raise ValueError("RAIL and ROAD source pair universe differs")
            rows=[]
            for a,b in combinations(sorted(included),2):
                k=(a,b)
                railrow=rsource[k];roadrow=dsource[k]
                rratio=railrow.get("rail_to_geodesic_ratio")
                dratio=roadrow.get("road_to_geodesic_ratio")
                railkm=railrow.get("cell_centroid_snap_augmented_rail_path_km")
                roadkm=roadrow.get("centroid_snap_augmented_road_path_km")
                geokm=railrow.get("source_cell_centroid_geodesic_km")
                if not all(isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
                           and x>0 for x in (rratio,dratio,railkm,roadkm,geokm)):
                    continue
                if min(rratio,dratio)<.95:
                    continue
                c1=included[a][period]
                c2=included[b][period]
                n1=c1["all_plant_indexed_records"]
                n2=c2["all_plant_indexed_records"]
                if n1<=0 or n2<=0:
                    raise ValueError("GBIF country flora record effort denominator invalid")
                raw_a=set(c1["counts"]);raw_b=set(c2["counts"])
                rawJ=len(raw_a&raw_b)/len(raw_a|raw_b)
                rare=expected_rarefied_overlap(c1["rarified_probs"],c2["rarified_probs"])
                if rare is None:
                    continue
                d1=datasets["regions"][block]["cells"][a]["periods"][period]
                d2=datasets["regions"][block]["cells"][b]["periods"][period]
                overlap=dataset_overlap(d1,d2)
                top1a=d1["top_datasetKey_source_share"];top1b=d2["top_datasetKey_source_share"]
                if not all(isinstance(x,(float,int)) and math.isfinite(x) and 0<=x<=1
                           for x in (top1a,top1b)):
                    raise ValueError("Missing source max-dataset dominance")
                rows.append({
                    "a":a,"b":b,
                    "rail_over_road":math.log(railkm/roadkm),
                    "baseline_cov":[math.log(geokm),math.log(roadkm),
                                    abs(math.log(n1/n2)),
                                    (math.log(n1)+math.log(n2))/2],
                    "dataset_overlap":overlap,
                    "same_top_dataset_key":d1["top_datasetKey_source_id"]==d2["top_datasetKey_source_id"],
                    "abs_difference_top_dataset_share":abs(top1a-top1b),
                    "rarefied20_indexed_source_overlap":rare,
                    "raw_indexed_source_species_Jaccard":rawJ
                })
            out={
                "n_distinct_source_cells_after_dataset_quality":len(included),
                "n_dependent_50km_cell_pair_rows_with_rail_road_and_dataset":len(rows),
                "median_shared_dataset_record_fraction":round(median(x["dataset_overlap"] for x in rows),6) if rows else None,
                "n_pairs_whose_largest_source_dataset_key_agrees":sum(x["same_top_dataset_key"] for x in rows),
                "n_different_largest_source_dataset_key_pairs":sum(not x["same_top_dataset_key"] for x in rows),
                "responses":{},
            }
            for resp in RESPONSES:
                x=[v["rail_over_road"] for v in rows]
                y=[v[resp] for v in rows]
                z=[v["baseline_cov"] for v in rows]
                overlap=[v["dataset_overlap"] for v in rows]
                topdiff=[v["abs_difference_top_dataset_share"] for v in rows]
                base_r=partial_rank_r(x,y,z)
                conditioned=partial_rank_r(x,y,[a+[b] for a,b in zip(z,overlap)])
                conditioned2=partial_rank_r(x,y,[a+[b,c] for a,b,c in zip(z,overlap,topdiff)])
                source_r=partial_rank_r(overlap,y,z)
                out["responses"][resp]={
                    "partial_rho_original_controls_on_IDENTICAL_dataset_supported_cell_pairs":base_r,
                    "partial_rho_add_allPlantae_dataset_composition_overlap":conditioned,
                    "partial_rho_add_overlap_and_max_dataset_dominance_difference":conditioned2,
                    "partial_rho_dataset_overlap_with_introduced_flora_similarity_given_base_controls":source_r,
                    "direction_negative_shorter_RAIL_route_coincides_with_more_similar_recorded_aliens":True,
                    "NO_causal_source_bias_adjustment_or_valid_pairwise_p_value":True
                }
            outcome["periods"][period][strategy]=out
            if strategy=="nearest_in_cell":
                main=out["responses"]["rarefied20_indexed_source_overlap"]
                print(json.dumps({
                    "block":block,"country":outcome["country"],"period":period,
                    "n_original_source_cells":outcome["n_original_qualified_two_period_flora_cells"],
                    "n_dataset_qualified_source_cells":len(included),
                    "n_true_source_pairs":len(rows),
                    "median_pair_dataset_record_overlap":out["median_shared_dataset_record_fraction"],
                    "rho_baseline_rarefied":main["partial_rho_original_controls_on_IDENTICAL_dataset_supported_cell_pairs"],
                    "rho_add_dataset_overlap_rarefied":main["partial_rho_add_allPlantae_dataset_composition_overlap"],
                    "rho_add_overlap_top1_rarefied":main["partial_rho_add_overlap_and_max_dataset_dominance_difference"],
                    "rho_dataset_source_overlap_with_flora":main["partial_rho_dataset_overlap_with_introduced_flora_similarity_given_base_controls"]
                },ensure_ascii=False),flush=True)
    return outcome


def country_equal_aggregate(regions,period,strategy,response,field):
    values={}
    audit={}
    for block,entry in regions.items():
        row=entry["periods"][period][strategy]
        value=row["responses"][response][field]
        nc=row["n_distinct_source_cells_after_dataset_quality"]
        npair=row["n_dependent_50km_cell_pair_rows_with_rail_road_and_dataset"]
        valid=(nc>=MIN_SOURCE_CELLS and npair>=MIN_VALID_PAIRS
               and isinstance(value,(int,float)) and not isinstance(value,bool)
               and math.isfinite(value) and -1<=value<=1)
        audit[block]={"source_cells":nc,"source_cell_pairs_NOT_independent":npair,"supported":valid}
        if valid:values[block]=float(value)
    vals=list(values.values())
    loo={k:round(median(z for n,z in values.items() if n!=k),6)
         for k in values if len(values)>1}
    return {
        "n_independent_qualifying_country_blocks":len(vals),
        "n_negative_railway_hypothesis_directions":sum(v<0 for v in vals),
        "n_positive_opposite_directions":sum(v>0 for v in vals),
        "median_country_equal_partial_rank":round(median(vals),6) if vals else None,
        "leave_one_country_out_medians":loo,
        "leave_country_out_median_sign_stable":bool(loo) and (
            all(x<0 for x in loo.values()) or all(x>0 for x in loo.values())),
        "regional_coefficients":values,
        "source_support":audit,
        "not_pooled_nonindependent_pair_pvalue":True,
    }


def analyze(flora,rail,road,datasets):
    if not all(v.get("complete") is True for v in (flora,rail,road,datasets)):
        raise ValueError("Incomplete independent real source data")
    if flora.get("source_errors"):
        raise ValueError("Original country GRIIS source errors")
    for sr in (rail,road,datasets):
        if set(sr.get("regions",{}))!=set(BLOCKS_SIX):
            raise ValueError("Source region universe does not match six original blocks")
    if not set(BLOCKS_SIX).issubset(flora.get("regions",{})):
        raise ValueError("Country-GRIIS source region missing")
    out={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "study":"POSTHOC source dataset-composition overlap vs GRIIS introduced flora similarity, conditional on physical RAIL/ROAD paths",
        "n_predeclared_country_macro_blocks":6,
        "source_dataset_marginal_ALL_Plantae_NOT_only_aliens":True,
        "not_observer_deduplicated_or_survey_adjusted":True,
        "regions":{},"country_equal_aggregates":{}
    }
    for block in BLOCKS_SIX:
        out["regions"][block]=compute_country(block,flora,rail,road,datasets)
    fields=(
        "partial_rho_original_controls_on_IDENTICAL_dataset_supported_cell_pairs",
        "partial_rho_add_allPlantae_dataset_composition_overlap",
        "partial_rho_add_overlap_and_max_dataset_dominance_difference",
        "partial_rho_dataset_overlap_with_introduced_flora_similarity_given_base_controls",
    )
    for p in PERIODS:
        out["country_equal_aggregates"][p]={}
        for strategy in STRATEGIES:
            out["country_equal_aggregates"][p][strategy]={
                response:{
                    field:country_equal_aggregate(out["regions"],p,strategy,response,field)
                    for field in fields
                } for response in RESPONSES
            }
    out["limitations"]=[
        "Dataset source overlap is a compositional proxy based on full ALL-PLANT recorded taxa, not country-listed introduced-only taxon surveys.",
        "Historical GBIF indexed flora taxon counts follow DETECTED marginal status while actual datasetKey facets include few NOT_DETECTED and unknown marginal records.",
        "A dataset ID is not an independent observer, site, effort, standardized plot or a justified causal detection adjustment.",
        "Observed taxon absences are unknown; no event-level GBIF DOI original record downloads or quality-controlled observer IDs.",
        "Cell pairs share sites and only one geographically selected source block per country; no false p-values or global causal conclusions.",
        "2026 source physical RAIL/ROAD network structure does not identify active historic freight or station propagule traffic in 2010-25.",
        "Adding posthoc covariates can create unstable rank residual coefficients; the result is a diagnostic rather than model selection for positive results.",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--flora",required=True)
    p.add_argument("--rail-sources",nargs="+",required=True)
    p.add_argument("--road-sources",nargs="+",required=True)
    p.add_argument("--dataset-sources",nargs="+",required=True)
    p.add_argument("--output",default="results/audit/GBIF_dataset_adjusted_6country_rail_plant_v51.json")
    args=p.parse_args()
    read=lambda paths:[json.loads(Path(x).read_text(encoding="utf-8")) for x in paths]
    flora=json.loads(Path(args.flora).read_text(encoding="utf-8"))
    rail=combine_source_regions(read(args.rail_sources))
    road=combine_source_regions(read(args.road_sources))
    data=merge_dataset_chunks(read(args.dataset_sources))
    result=analyze(flora,rail,road,data)
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    for p in PERIODS:
        a=result["country_equal_aggregates"][p]["nearest_in_cell"]["rarefied20_indexed_source_overlap"]
        print(json.dumps({"period":p,
            "baseline_identical_pairs":a["partial_rho_original_controls_on_IDENTICAL_dataset_supported_cell_pairs"],
            "dataset_overlap_adjusted":a["partial_rho_add_allPlantae_dataset_composition_overlap"],
            "dataset_overlap_and_dominance_adjusted":a["partial_rho_add_overlap_and_max_dataset_dominance_difference"]
        },ensure_ascii=False),flush=True)
    print(json.dumps({"complete":True,"n_blocks":len(result["regions"]),
                      "posthoc_source_data_quality_not_causal":True}),flush=True)


if __name__=="__main__":
    main()
