#!/usr/bin/env python3
"""v2.7 noncausal partial road-covariate control for introduced plant × rail graph.

Combine 4 independently source-verified artifacts:
 GBIF full source taxon facets by geographic cell x period (v2.0)
 GRIIS country introduced species backbone keys (v2.4)
 Overture physical rail connector graph (v2.1)
 Overture clipped geodesic road/rail length density per source cell (v2.6)

Compare physical rail-connected vs mapped-unconnected cell pairs on
observed introduced species-key similarity, with **road density**
matching in addition to geographic distance and indexed recording effort.
Road kilometers are NOT network travel times; cannot equate matching with
full road-network confounding control. Country blocks, not cell pairs,
are the independent replication unit. No seed-transport causal effect.
"""
from __future__ import annotations
import argparse
import json
from datetime import datetime,timezone
from itertools import combinations
from math import isfinite, log
from pathlib import Path
from statistics import median

from analyze_griis_alien_rail_source_v25 import (
    valid_griis,cell_observations,ALIEN_RAREFACTION_RECORDS
)
from analyze_macro_rail_assemblages_v22 import (
    COUNTRIES,PERIODS,expected_rarefied_overlap,
    geographic_distance_km,jaccard
)

ROAD_PSEUDO_DENSITY=0.02
N_MATCH_MIN=3
N_INDEPENDENT_CITY_MIN=5


def _finite_nonnegative(x):
    return isinstance(x,(float,int)) and not isinstance(x,bool) and isfinite(x) and x>=0


def validate_transport_city(transport,plant,rail):
    if transport.get("country")!=plant.get("country") or transport.get("country")!=rail.get("country"):
        raise ValueError("Country source transport mismatch")
    cells=transport.get("cells")
    if not isinstance(cells,dict) or set(cells)!=set(plant.get("cells",{})):
        raise ValueError("Road covariate source cell frame mismatches GBIF")
    if set(cells)!=set(rail.get("cell_component_cover",{})):
        raise ValueError("Road covariate source cell frame mismatches mapped rail graph")
    metrics={}
    for cid,v in cells.items():
        km=v.get("road_km")
        railkm=v.get("rail_km")
        area=v.get("area_km2")
        if not all(_finite_nonnegative(x) for x in (km,railkm,area)) or area<=0:
            raise ValueError(f"Invalid real road and rail lengths or cell area: {cid}")
        road_density=km/area
        rail_density=railkm/area
        from_source=v.get("road_length_density_km_per_km2")
        if not _finite_nonnegative(from_source) or abs(road_density-from_source)>1e-5:
            raise ValueError("Road length per area inconsistent with source")
        metrics[cid]={"road_density":road_density,"rail_density":rail_density}
    return metrics


def build_city_pairs(plant,rail,transport,alien_keys,period):
    metrics=validate_transport_city(transport,plant,rail)
    taxa,_=cell_observations(plant,alien_keys,period)
    lookup={}
    for e in rail["cell_pair_connection_descriptors"]:
        k=tuple(sorted([e["a"],e["b"]]))
        if k in lookup or k[0]==k[1]:
            raise ValueError("Duplicate rail network source pair")
        lookup[k]=e.get("within_clipped_graph_connected")
    if len(lookup)!=36 or any(not isinstance(v,bool) for v in lookup.values()):
        raise ValueError("Incomplete rail connectivity source")
    rows=[]
    ids=sorted(taxa)
    for a,b in combinations(ids,2):
        p,q=taxa[a],taxa[b]
        conn=lookup[(a,b)]
        road_a,road_b=(metrics[x]["road_density"] for x in (a,b))
        rail_a,rail_b=(metrics[x]["rail_density"] for x in (a,b))
        if not road_a>0 or not road_b>0:
            # Source missing road density is *not* imputed to a genuine
            # road-free habitat: omit from this control-identified analysis.
            continue
        share=expected_rarefied_overlap(p["probs"],q["probs"])
        if share is None:continue
        x={"a":a,"b":b,
           "connected":conn,
           "both_connector_source_rail":bool(
               rail["cell_component_cover"][a]["n_source_components"]
               and rail["cell_component_cover"][b]["n_source_components"]),
           "both_actual_source_rail_length":rail_a>0 and rail_b>0,
           "distance_km":geographic_distance_km(p["bbox"],q["bbox"]),
           "log_record_effort_ratio":abs(log(p["n_total_source_records"]/q["n_total_source_records"])),
           "log_mean_record_count":(log(p["n_total_source_records"])+log(q["n_total_source_records"]))/2,
           "road_mean_log_density":(log(road_a+ROAD_PSEUDO_DENSITY)+log(road_b+ROAD_PSEUDO_DENSITY))/2,
           "road_absolute_log_ratio":abs(log((road_a+ROAD_PSEUDO_DENSITY)/(road_b+ROAD_PSEUDO_DENSITY))),
           "source_rarefied20_introduced_similarity":share,
           "source_unrarefied_introduced_jaccard":jaccard(p["counts"],q["counts"]),
        }
        rows.append(x)
    qc={
        "n_record_eligible_intro_plant_cells":len(ids),
        "n_cell_rail_presence_disagreements_connector_vs_clipped_length":sum(
            ((rail["cell_component_cover"][cid]["n_source_components"]>0)
             != (metrics[cid]["rail_density"]>0))
            for cid in ids
        ),
        "n_cell_pairs_pass_both_cell_road_density_positive":len(rows)
    }
    return rows,qc


def compare_matched_pairs(rows, require_both_rails=False, match_road=False):
    """Greedy matching fixed by source covariates, never response similarity."""
    source=rows if not require_both_rails else [r for r in rows if r["both_connector_source_rail"]]
    pos=sorted([r for r in source if r["connected"]],key=lambda x:(x["distance_km"],x["a"],x["b"]))
    neg=[r for r in source if not r["connected"]]
    used=set();pair_deltas=[]
    for x in pos:
        opts=[]
        for i,y in enumerate(neg):
            if i in used:continue
            distdiff=abs(x["distance_km"]-y["distance_km"])
            distmax=max(5.,.33*max(x["distance_km"],y["distance_km"]))
            if distdiff>distmax:continue
            meandiff=abs(x["log_mean_record_count"]-y["log_mean_record_count"])
            obsratio=abs(x["log_record_effort_ratio"]-y["log_record_effort_ratio"])
            if meandiff>1. or obsratio>1.:continue
            road_mean_diff=abs(x["road_mean_log_density"]-y["road_mean_log_density"])
            road_balance_diff=abs(x["road_absolute_log_ratio"]-y["road_absolute_log_ratio"])
            if match_road and (road_mean_diff>.5 or road_balance_diff>.7):
                continue
            score=(distdiff/distmax)+meandiff+obsratio
            if match_road:score+=road_mean_diff+road_balance_diff
            opts.append((score,i))
        if opts:
            _,idx=min(opts)
            used.add(idx)
            y=neg[idx]
            pair_deltas.append(x["source_rarefied20_introduced_similarity"]-
                               y["source_rarefied20_introduced_similarity"])
    adequate=len(pair_deltas)>=N_MATCH_MIN and len(pos)>=N_MATCH_MIN and len(neg)>=N_MATCH_MIN
    return {
        "n_connected_source_pairs":len(pos),
        "n_unconnected_source_pairs":len(neg),
        "n_covariate_matched_source_pair_contrasts":len(pair_deltas),
        "median_city_pair_similarity_contrast":round(median(pair_deltas),6) if adequate else None,
        "source_support_gate_pass":adequate,
        "treatment": "both cells have Overture connector source" if require_both_rails else "all source-eligible cells",
        "source_match_includes_road_intensity":match_road,
        "not_road_network_distance_or_train_causal_effect":True
    }


def city_summary(plant,rail,road,keys,period):
    rows,qc=build_city_pairs(plant,rail,road,keys,period)
    results={}
    for label,kwargs in (
        ("distance_effort_all",{}),
        ("road_intensity_all",{"match_road":True}),
        ("distance_effort_both_rail",{"require_both_rails":True}),
        ("road_intensity_both_rail",{"require_both_rails":True,"match_road":True})
    ):
        results[label]=compare_matched_pairs(rows,**kwargs)
    return {"source_coverage":qc,"comparisons":results}


def aggregate_blocks(regions,period,label):
    valid=[]
    for name,entry in regions.items():
        p=entry["periods"].get(period)
        if not p:continue
        c=p["comparisons"][label]
        if c["source_support_gate_pass"]:
            valid.append((name,c["median_city_pair_similarity_contrast"]))
    vals=[v for n,v in valid]
    return {
        "n_informative_city_blocks":len(vals),
        "positive_blocks":sum(v>0 for v in vals),
        "negative_blocks":sum(v<0 for v in vals),
        "median_city_block_source_contrast":round(median(vals),6) if vals else None,
        "country_block_minimum_for_evaluation":N_INDEPENDENT_CITY_MIN,
        "minimum_geographic_replication_met":len(vals)>=N_INDEPENDENT_CITY_MIN,
        "not_a_p_value_or_causal_effect":True
    }


def analyze(plant,rail,griis,roads):
    if not (plant.get("source_complete") and rail.get("complete")
            and roads.get("complete")):
        raise ValueError("Incomplete prerequisite geography, rail or road source")
    expected=set(COUNTRIES)
    if any(set(r.get("regions",{}))!=expected for r in (plant,rail,roads)):
        raise ValueError("Ten-country source frame mismatch")
    aliases=valid_griis(griis)
    out={"utc":datetime.now(timezone.utc).isoformat(),
         "study":"GRIIS introduced plant source composition vs rail connection conditional on road INTENSITY",
         "regions":{},"country_block":{},"rarefied_plant_source_records":ALIEN_RAREFACTION_RECORDS,
         "not_an_estimate_of_railway_propagule_transport":True}
    for name,country in COUNTRIES.items():
        if plant["regions"][name]["country"]!=country or rail["regions"][name]["country"]!=country or roads["regions"][name]["country"]!=country:
            raise ValueError("Source city-country mismatch")
        out["regions"][name]={"country":country,"periods":{}}
        keys=aliases.get(country)
        if keys is None:
            out["regions"][name]["source_exclusion"]="no adequate source GRIIS checklist, taxon status unknown"
            continue
        for period in PERIODS:
            result=city_summary(plant["regions"][name],rail["regions"][name],
                                roads["regions"][name],keys,period)
            out["regions"][name]["periods"][period]=result
            a=result["comparisons"]["distance_effort_all"]
            b=result["comparisons"]["road_intensity_all"]
            print(json.dumps({
                "city":name,"country":country,"period":period,
                "n_eligible_species_cells":result["source_coverage"]["n_record_eligible_intro_plant_cells"],
                "rail_presence_spatial_disagreements":result["source_coverage"]["n_cell_rail_presence_disagreements_connector_vs_clipped_length"],
                "base_matched":a["n_covariate_matched_source_pair_contrasts"],
                "road_density_matched":b["n_covariate_matched_source_pair_contrasts"],
                "base_delta":a["median_city_pair_similarity_contrast"],
                "road_density_delta":b["median_city_pair_similarity_contrast"],
            },ensure_ascii=False),flush=True)
    for period in PERIODS:
        out["country_block"][period]={
            l:aggregate_blocks(out["regions"],period,l)
            for l in ("distance_effort_all","road_intensity_all",
                      "distance_effort_both_rail","road_intensity_both_rail")
        }
    out["limits"]=[
        "Road length density is not network reachability or measured road traffic. This is a *partial* road confound screen.",
        "Plant taxon presence is indexed recording, not systematic repeated ecological census or true absence.",
        "GRIIS nonmembership means unknown, not native, and checklists have different completeness/vintages.",
        "No actual freight train activity, railside seed trapping, grass-cutting, wind treatment or adaptation genetics measured.",
        "Source rail network is 2026 while occurrence periods end 2017/2025; no historical treatment chronology.",
        "Cell pair observations share sites; only city blocks are independent summaries. No pseudo-independent p-values.",
        "One selected urban study window per country and sub-50km rectangular cells are not 25/50km equal-area ecological regional replicates.",
        "If road-density calipers remove support, omit estimates rather than interpreting null as zero effect."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gbif-source",required=True)
    p.add_argument("--rail-source",required=True)
    p.add_argument("--griis-source",required=True)
    p.add_argument("--transport-source",required=True)
    p.add_argument("--output",default="results/audit/alien_rail_road_partial_control_v27.json")
    a=p.parse_args()
    source=[json.loads(Path(x).read_text(encoding="utf-8")) for x in (
        a.gbif_source,a.rail_source,a.griis_source,a.transport_source)]
    result=analyze(*source)
    o=Path(a.output);o.parent.mkdir(parents=True,exist_ok=True)
    o.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"blocks":result["country_block"],
                      "not_railway_seed_transport":result["not_an_estimate_of_railway_propagule_transport"]},
                     ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
