#!/usr/bin/env python3
"""v4.4 independent six-country 50km macro RAIL vs ROAD topology/plant survey audit.

After v3.4's two-country source result (US Midwest and France Burgundy)
had opposite recent-period signs, four other independent countries were
selected from the EXISTING frozen 10-block country-masked GRIIS data
on SOURCE COVERAGE ONLY, before reading their rail/plant route outcomes:
Canada Ontario, Britain Yorkshire, Germany Saxony, Australia Southeast.

All analysis inherits the exact original v3.4 rarity/record drift gates,
two source-cell network anchor strategies, two time periods and both
indexed species-key similarity definitions. No changing sample grid,
calipers or outcomes to chase a positive sign.

One COUNTRY-BLOCK is the summary replicate, not its 120 overlapping
pairwise cell comparisons. No frequentist p-value is computed. This is
still 2026 mapped routes vs 2010–25 GBIF-positive RECORDING patterns,
not train-operated seed transfer or verified plant colonization.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
from statistics import median

from analyze_macro_50km_rail_road_aliens_v34 import (
    analyze_one_block, STRATEGIES
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS

EXISTING_BLOCKS=("us_midwest","fr_burgundy")
INDEPENDENT_ADDITION=("ca_ontario","gb_yorkshire","de_saxony","au_southeast")
BLOCKS_SIX=EXISTING_BLOCKS+INDEPENDENT_ADDITION
PERIODS=("2010_2017","2018_2025")
RESPONSES=("source_expected_rarefied20_intersection_union_ratio",
           "source_observed_alien_jaccard")
MIN_SOURCE_CELLS=8
MIN_VALID_PAIRS=20
REPLICATION_GATE=5


def combine_source_regions(sources,allowed=BLOCKS_SIX):
    """Merge successful source regions without losing provenance or duplicates."""
    if not sources or not isinstance(sources,list):
        raise ValueError("No original source chunks")
    result={}
    source_count=0
    for chunk in sources:
        if chunk.get("complete") is not True or chunk.get("errors") or chunk.get("source_errors"):
            raise ValueError("Source region batch incomplete or errored; fail closed")
        regions=chunk.get("regions")
        if not isinstance(regions,dict) or not regions:
            raise ValueError("Missing source regions in chunk")
        for block,row in regions.items():
            if block not in allowed:
                raise ValueError("Non-predeclared source region included")
            if block in result:
                raise ValueError("Same country source region appeared twice")
            if not isinstance(row,dict) or "routes" not in row:
                raise ValueError("Missing actual network route source")
            expected_country=BLOCKS[block][0]
            given_country=row.get("country") or row.get("source",{}).get("country")
            if given_country!=expected_country:
                raise ValueError("Actual source country identification mismatch")
            result[block]=row
            source_count+=1
    if set(result)!=set(allowed):
        raise ValueError("Six-country network source incomplete: expected exact predeclared regions")
    return {"complete":True,"regions":result,
            "n_true_source_region_chunks":len(sources),
            "n_geographically_independent_blocks":source_count}


def aggregate(blocks,period,strategy,response):
    """Non-independent pairs distilled to ONE direction per source-country block."""
    vals={}
    audit={}
    for block,row in blocks.items():
        if period not in row["periods"] or strategy not in row["periods"][period]:
            raise ValueError("Original regional result strategy/period not defined")
        entry=row["periods"][period][strategy]
        nsource=row["n_source_qualified_two_period_plant_cells"]
        npairs=entry["n_plant_pair_routes_with_both_network_lengths"]
        ref=entry["responses"].get(response)
        if ref is None:
            raise ValueError("Unexpected species-key response source")
        rho=ref["descriptive_partial_rank_correlation"]
        supported=(isinstance(nsource,int) and nsource>=MIN_SOURCE_CELLS and
                   isinstance(npairs,int) and npairs>=MIN_VALID_PAIRS and
                   rho is not None and isinstance(rho,(float,int)) and
                   not isinstance(rho,bool) and math.isfinite(rho) and -1<=rho<=1)
        audit[block]={"source_qualified_cells":nsource,
                      "source_pairs_sharing_cells":npairs,
                      "rank_correlation_identified":rho is not None,
                      "regional_source_support":supported,
                      "country":BLOCKS[block][0]}
        if supported:vals[block]=float(rho)
    values=list(vals.values())
    med=median(values) if values else None
    loo={k:round(median(v for j,v in vals.items() if j!=k),6)
         for k in vals if len(vals)>1}
    out={
        "n_original_country_macro_blocks":len(blocks),
        "n_independent_regions_with_qualified_network_and_flora":len(vals),
        "negative_rho_expected_if_relatively_shorter_rail_paths_share_MORE_introduced_taxa":True,
        "n_blocks_in_hypothesis_compatible_negative_direction":sum(v<0 for v in values),
        "n_blocks_in_opposite_positive_direction":sum(v>0 for v in values),
        "n_blocks_exact_zero":sum(v==0 for v in values),
        "country_equal_weight_median_rho":round(med,6) if med is not None else None,
        "leave_one_independent_country_out_medians":loo,
        "leave_country_out_sign_stable":bool(loo) and (all(v<0 for v in loo.values()) or all(v>0 for v in loo.values())),
        "minimum_five_independent_regions_for_descriptive_generalization":len(vals)>=REPLICATION_GATE,
        "n_geographically_independent_country_blocks_not_pairwise_pseudoreplicates":len(vals),
        "regional_rhos":vals,
        "pre_outcome_source_support_audit":audit,
        "has_no_hypothesis_significance_p_value":True
    }
    return out


def analyze(flora,rail,road):
    if flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("Country-masked GRIIS source incomplete")
    if not set(BLOCKS_SIX).issubset(flora.get("regions",{})):
        raise ValueError("Frozen six block source GRIIS taxonomy missing")
    if rail.get("complete") is not True or road.get("complete") is not True:
        raise ValueError("Network source complete flag missing")
    if set(rail.get("regions",{}))!=set(BLOCKS_SIX) or set(road.get("regions",{}))!=set(BLOCKS_SIX):
        raise ValueError("Rail/road independent geographic source frames disagree")
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "research_stage":"EXPLORATORY six independent-country source blocks versus two-block baseline",
        "six_country_source_geographies":list(BLOCKS_SIX),
        "independent_additions_chosen_from_frozen_10block_8country_source_by_RECORD_COVERAGE_before_their_network_outcomes":list(INDEPENDENT_ADDITION),
        "original_test_blocks":list(EXISTING_BLOCKS),
        "response_definitions":list(RESPONSES),
        "observational_and_dependent_pairwise_no_p_values":True,
        "not_rail_freight_train_seed_dispersal_effect":True,
        "region_results":{},
        "country_equal_weight_aggregate":{},
    }
    for block in BLOCKS_SIX:
        before=flora["regions"][block].get("country")
        actual=BLOCKS[block][0]
        if before!=actual or rail["regions"][block].get("country")!=actual or road["regions"][block].get("source",{}).get("country")!=actual:
            raise ValueError("Actual source country mismatch")
        result=analyze_one_block(block,flora,rail,road)
        out["region_results"][block]=result
        for period in PERIODS:
            for strategy in STRATEGIES:
                r=result["periods"][period][strategy]
                x=r["responses"][RESPONSES[0]]["descriptive_partial_rank_correlation"]
                print(json.dumps({
                    "block":block,"country":actual,
                    "period":period,"anchor":strategy,
                    "n_source_cells":result["n_source_qualified_two_period_plant_cells"],
                    "n_paired_network_species_rows":r["n_plant_pair_routes_with_both_network_lengths"],
                    "source_rarefied20_partial_rho":x,
                    "source_jaccard_partial_rho":r["responses"][RESPONSES[1]]["descriptive_partial_rank_correlation"],
                },ensure_ascii=False),flush=True)
    for period in PERIODS:
        out["country_equal_weight_aggregate"][period]={}
        for strategy in STRATEGIES:
            out["country_equal_weight_aggregate"][period][strategy]={
                response:aggregate(out["region_results"],period,strategy,response)
                for response in RESPONSES
            }
    out["data_and_inference_limits"]=[
        "Source plant observations are positively recorded GBIF species-key facets, not event-deduplicated ecological plots or verified absence.",
        "Country GRIIS introduced status varies in taxonomic coverage and vintage, with nonmembers unknown rather than native.",
        "Country-area >=80% masks exclude some coastal/border cells but are not verified terrestrial ecoregion or vegetation masks.",
        "Rail and road Overture are 2026 physical mapped paths; not historical 2010-2017 or 2018-2025 active freight routes.",
        "Rail/road network anchoring to 50km cell centroid is a source accessibility proxy, not sampled plant locations.",
        "The rail/road source geometric paths use physical 75km-halo graphs, not actual vehicle traffic, turns or train schedules.",
        "Only one frozen broad geographical block per country in this six-country exploratory contrast, not independent replicated regions within each country.",
        "Partial rank association is purely descriptive, confounded by unknown climate, biome, land use, ports, observer effort and freight.",
        "No p-value and no false n from pairwise 120 observations sharing the same source 16 cells.",
        "Country selection was based on recorded source coverage from the existing ten-block frame, not a random global sample.",
    ]
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--flora",required=True)
    ap.add_argument("--rail-sources",nargs="+",required=True)
    ap.add_argument("--road-sources",nargs="+",required=True)
    ap.add_argument("--output",default="results/audit/6country_50km_rail_vs_road_plant_v44.json")
    args=ap.parse_args()
    flora=json.loads(Path(args.flora).read_text(encoding="utf-8"))
    rail=combine_source_regions([json.loads(Path(p).read_text(encoding="utf-8")) for p in args.rail_sources])
    road=combine_source_regions([json.loads(Path(p).read_text(encoding="utf-8")) for p in args.road_sources])
    result=analyze(flora,rail,road)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"n_country_source_blocks":len(result["region_results"]),
        "main_recent_rarefied":result["country_equal_weight_aggregate"]["2018_2025"]["nearest_in_cell"][RESPONSES[0]],
        "NO_train_seed_dispersal_inferred":result["not_rail_freight_train_seed_dispersal_effect"],
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
