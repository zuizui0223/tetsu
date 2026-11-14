#!/usr/bin/env python3
"""v4.6 six-independent-country macro GRIIS flora vs RAIL circuit structure.

Extends UNCHANGED 2010-17/2018-25 GBIF speciesKey source, two fixed source
anchor strategies and original rank controls from v3.8 to four NEW
country blocks preselected by botanical recording quality:
US Midwest, France Burgundy, CA Ontario, GB Yorkshire, DE Saxony,
AU Southeast. These six source-country blocks are independent
regional units, but source 50km cell pairs within a block are DEPENDENT.

No new tuning of edge resistances to plant associations; source
effective resistance = mapped rail km-equivalent circuit model, and
a ROAD shortest route is held as a control. Unlike v4.2, major ROAD
circuit Q is only observed in US/FR and NOT adjusted for all 6.
Do not treat raw source country-list incidence as actual plot absence
or a train-carried propagule mechanism.
"""
from __future__ import annotations
import argparse,json,math
from datetime import datetime,timezone
from pathlib import Path
from statistics import median

from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_6country_50km_network_plants_v44 import (
    BLOCKS_SIX,EXISTING_BLOCKS,INDEPENDENT_ADDITION,
    STRATEGIES,PERIODS,combine_source_regions,
    MIN_SOURCE_CELLS,MIN_VALID_PAIRS,REPLICATION_GATE
)
from analyze_macro_alien_circuit_increment_v38 import analyze_city

RESPONSES=("rarefied20_source_similarity","raw_observed_key_jaccard")
Q_FIELD="partial_rank_extra_path_redundancy_Q_at_fixed_shortest_rail"
SHORTEST_FIELD="partial_rank_shortest_rail_vs_road"
EFF_FIELD="partial_rank_Circuitscape_equivalent_effective_rail_vs_road"


def merge_circuit_source_chunks(chunks):
    if not chunks or not isinstance(chunks,list):
        raise ValueError("No physical source rail circuit artifacts")
    regions={}
    exclusions=[]
    for source in chunks:
        if source.get("complete") is not True or source.get("errors"):
            raise ValueError("Incomplete original 50km circuit source")
        rows=source.get("regions")
        if not isinstance(rows,dict) or not rows:
            raise ValueError("Missing source resistor regions")
        for block,record in rows.items():
            if block=="jp_tohoku":
                exclusions.append(block)
                continue
            if block not in BLOCKS_SIX:
                raise ValueError("Unfrozen network source blocks included")
            if block in regions:
                raise ValueError("Duplicate independent country circuit source")
            if record.get("country")!=BLOCKS[block][0] or record.get("actual_Julia_Circuitscape_executed") is not False:
                raise ValueError("Source country mismatch or unsupported real-Julia claim")
            if not isinstance(record.get("strategies"),dict) or set(record["strategies"])!=set(STRATEGIES):
                raise ValueError("Original physical network circuit anchor source changed")
            regions[block]=record
    if set(regions)!=set(BLOCKS_SIX):
        raise ValueError("Six country effective resistance source incomplete")
    return {"complete":True,"regions":regions,
            "predeclared_known_upstream_extra_JP_excluded":exclusions}


def blocks_summary(by_block,period,strategy,response,field):
    vals={}
    audit={}
    for name,v in by_block.items():
        source=v["periods"][period][strategy]
        r=source["responses"][response].get(field)
        n=source["n_dependent_cell_pair_rows"]
        nc=v["n_two_period_eligible_GRIIS_50km_cells"]
        qualifies=(isinstance(nc,int) and nc>=MIN_SOURCE_CELLS
                   and isinstance(n,int) and n>=MIN_VALID_PAIRS
                   and isinstance(r,(float,int)) and not isinstance(r,bool)
                   and math.isfinite(r) and -1<=r<=1)
        audit[name]={"n_valid_50km_source_cells":nc,
                     "n_source_pairwise_observations_sharing_cells":n,
                     "country":BLOCKS[name][0],
                     "partial_rank_identifiable":r is not None,
                     "block_support":qualifies}
        if qualifies:vals[name]=float(r)
    loo={name:round(median(v for block,v in vals.items() if block!=name),6)
         for name in vals if len(vals)>1}
    values=list(vals.values())
    return {
        "n_frozen_country_macro_blocks":len(by_block),
        "n_independent_country_blocks_with_positive_network_and_flora_source":len(vals),
        "positive_direction_if_Q_increases_observed_introduced_species_similarity":True,
        "n_blocks_in_railway_redundancy_hypothesis_positive_direction":sum(r>0 for r in values),
        "n_blocks_opposite_negative_direction":sum(r<0 for r in values),
        "n_blocks_zero":sum(r==0 for r in values),
        "country_equal_weight_median_descriptive_partial_rank":round(median(values),6) if values else None,
        "leave_one_country_out_medians":loo,
        "leave_one_out_direction_stable":bool(loo) and (all(v>0 for v in loo.values()) or all(v<0 for v in loo.values())),
        "at_least_5_independent_country_blocks_source_support":len(vals)>=REPLICATION_GATE,
        "regional_partial_rhos":vals,
        "source_quality_and_identification":audit,
        "NO_pairwise_p_value_or_propagule_causality":True
    }


def run(flora,rail,road,circuit):
    if not all(v.get("complete") is True for v in (flora,rail,road,circuit)):
        raise ValueError("Predefined source complete flags missing")
    if not set(BLOCKS_SIX).issubset(flora.get("regions",{})):
        raise ValueError("Original country-GRIIS 10block source missing requested six")
    for item in (rail,road,circuit):
        if set(item.get("regions",{}))!=set(BLOCKS_SIX):
            raise ValueError("Physical railway / road source has changed country geography")
    result={
        "utc":datetime.now(timezone.utc).isoformat(),
        "actual_Julia_Circuitscape_was_run_on_SYNTHETIC_TESTS_only":True,
        "actual_50km_source_networks_electrical_solver_SciPy":True,
        "study_kind":"Posthoc descriptive six-country GRIIS-introduced indexed flora vs physical RAIL circuit redundancy conditional on ROAD shortest routes",
        "predefined_country_blocks":list(BLOCKS_SIX),
        "new_country_blocks_selected_by_coverage_before_their_rail_outcomes":list(INDEPENDENT_ADDITION),
        "region_results":{},"country_equal_weight_results":{},
        "NO_railway_seed_dispersal_identified":True
    }
    for name in BLOCKS_SIX:
        v=analyze_city(name,flora,rail,road,circuit,allowed_blocks=BLOCKS_SIX)
        result["region_results"][name]=v
        for period in PERIODS:
            for anchor in STRATEGIES:
                row=v["periods"][period][anchor]
                response=row["responses"][RESPONSES[0]]
                print(json.dumps({
                    "block":name,"country":BLOCKS[name][0],
                    "period":period,"anchor":anchor,
                    "n_source_qualified_50km_cells":v["n_two_period_eligible_GRIIS_50km_cells"],
                    "n_GRIIS_flora_source_pairs_with_valid_circuit":row["n_dependent_cell_pair_rows"],
                    "partial_rho_Q_added_given_shortest_rail":response[Q_FIELD],
                    "partial_rho_Circuitscape_eff_over_ROAD":response[EFF_FIELD],
                    "partial_rho_shortest_rail_over_ROAD":response[SHORTEST_FIELD],
                    "raw_species_Jaccard_partial_rho_Q":row["responses"][RESPONSES[1]][Q_FIELD],
                },ensure_ascii=False),flush=True)
    for period in PERIODS:
        result["country_equal_weight_results"][period]={}
        for anchor in STRATEGIES:
            result["country_equal_weight_results"][period][anchor]={
                response:{
                    metric:blocks_summary(result["region_results"],period,anchor,response,metric)
                    for metric in (Q_FIELD,SHORTEST_FIELD,EFF_FIELD)
                } for response in RESPONSES
            }
    result["inference_limitations"]=[
        "Source physical edge resistance is railway segment length in km-equivalent units, not train traffic, active freight or measured seed movement.",
        "All circuit source values were derived in SciPy but Julia Circuitscape program only had independent small synthetic network parity.",
        "Only one 200km selected source area per country; six countries are not random global terrestrial ecoregion replicates.",
        "ROAD circuit Q exists only for earlier US/FR source, so road Q is NOT a common six-country confound control. The common road predictor here is Dijkstra physical path distance.",
        "GBIF taxon facets are indexed positive source records with no survey-verified species absence, no original DOI event deduplication, observer effort or coordinate uncertainty.",
        "Plant data are country GRIIS-listed introduced taxa only; unlisted flora cannot be treated as native.",
        "The 2026 rail map postdates 2010-17 and 2018-25 plant records; no treatment onset or train seed colonization chronology is identified.",
        "Cell source-pair rank correlations are dependent because 50km cells appear in many pairs; no per-pair p-values or fake 211 observations as independent regional replication.",
        "True climate/soil/ecoregion/urban/port/source reporting controls remain missing; across-country source-list completeness differs.",
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--flora",required=True)
    p.add_argument("--rail-sources",nargs="+",required=True)
    p.add_argument("--road-sources",nargs="+",required=True)
    p.add_argument("--circuit-sources",nargs="+",required=True)
    p.add_argument("--output",default="results/audit/6country_50km_GRIIS_circuit_increment_v46.json")
    a=p.parse_args()
    read=lambda items:[json.loads(Path(x).read_text(encoding="utf-8")) for x in items]
    flora=json.loads(Path(a.flora).read_text(encoding="utf-8"))
    rail=combine_source_regions(read(a.rail_sources))
    road=combine_source_regions(read(a.road_sources))
    circuit=merge_circuit_source_chunks(read(a.circuit_sources))
    result=run(flora,rail,road,circuit)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "source_blocks":len(result["region_results"]),
        "late_period_main_rarefied_Q":result["country_equal_weight_results"]["2018_2025"]["nearest_in_cell"]["rarefied20_source_similarity"][Q_FIELD],
        "no_railway_seed_mechanism_claim":result["NO_railway_seed_dispersal_identified"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
