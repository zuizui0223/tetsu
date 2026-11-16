#!/usr/bin/env python3
"""v2.10 five-source introduced plant × RAIL connectivity versus ROAD PATH.

Reuse exact observed GRIIS/GBIF 90 cell × period source, mapped physical
rail graph, measured road km density and new actual Overture physical
ROAD connector shortest-path proxies. Match *connected vs disconnected
railway* pairs on road geography and record-count effort, now including
road network *path length* and network circuity, not just length density.

This is exploratory and posthoc to prior v2.7 outcome; report all source
eligibility and a whole prespecified grid of matching calipers, NEVER
select the caliper which gives a wanted sign or p-value. One chosen city
per country is the replicate, not overlapping cell pairs.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from statistics import median

from analyze_griis_alien_rail_source_v25 import valid_griis
from analyze_griis_rail_with_road_density_v27 import build_city_pairs
from analyze_macro_rail_assemblages_v22 import COUNTRIES, PERIODS

SENSITIVITY_CALIPERS={
    "tight":(0.25,0.25),
    "main":(0.45,0.40),
    "loose":(0.70,0.70),
}
MIN_PER_CITY_MATCHES=3
MIN_INDEPENDENT_CITIES=5
ROAD_MAX_CIRCUITY=10.


def validate_road_graph(city_name,country,road_graph,rail_graph):
    if road_graph.get("country")!=country or rail_graph.get("country")!=country:
        raise ValueError("Road and rail country registry disagrees")
    rows=road_graph.get("cell_pair_road_paths")
    if not isinstance(rows,list) or len(rows)!=36:
        raise ValueError(f"Incomplete source road paths: {city_name}")
    r_pairs=rail_graph.get("cell_pair_connection_descriptors")
    if not isinstance(r_pairs,list) or len(r_pairs)!=36:
        raise ValueError(f"Incomplete source rail connectivity: {city_name}")
    rail_ids={tuple(sorted([r["a"],r["b"]])) for r in r_pairs}
    result={}
    for r in rows:
        a,b=r.get("a"),r.get("b")
        if not isinstance(a,str) or not isinstance(b,str) or a==b:
            raise ValueError("Invalid source city pair names")
        key=tuple(sorted((a,b)))
        if key in result or key not in rail_ids:
            raise ValueError("Road graph geographic cells disagree with rail")
        if not isinstance(r.get("within_clipped_road_graph_physically_connected"),bool):
            raise ValueError("Road connectivity not source-verified")
        length=r.get("optimistic_road_connector_path_km")
        circuity=r.get("optimistic_road_to_euclidean_ratio")
        good=(r.get("within_clipped_road_graph_physically_connected") is True
              and r.get("source_centroids_snapped_within_5km") is True
              and isinstance(length,(int,float)) and not isinstance(length,bool)
              and isinstance(circuity,(int,float)) and not isinstance(circuity,bool)
              and math.isfinite(length) and math.isfinite(circuity)
              and length>0 and 0.85<=circuity<=ROAD_MAX_CIRCUITY)
        result[key]={
            "network_length_km":float(length) if good else None,
            "circuity":float(circuity) if good else None,
            "road_network_source_usable":good,
            "mapped_road_connected":r["within_clipped_road_graph_physically_connected"]
        }
    if set(result)!=rail_ids:
        raise ValueError("Road graph incomplete pair support")
    return result


def add_road_paths(city,country,rail_city,road_city):
    path_lookup=validate_road_graph(city,country,road_city,rail_city)
    for rec in path_lookup.values():
        if rec["road_network_source_usable"] and not rec["mapped_road_connected"]:
            raise ValueError("Invalid positive road path on disconnected graph")
    return path_lookup


def source_pair_records(plant,rail,griis,transport,road_graph,city,period):
    aliases=valid_griis(griis).get(COUNTRIES[city])
    if aliases is None:
        return [],{"source_excluded":"missing/insufficient country GRIIS checklist"}
    plants=plant["regions"][city]
    railway=rail["regions"][city]
    roads=transport["regions"][city]
    routes=road_graph["regions"][city]
    if any(x["country"]!=COUNTRIES[city] for x in (plants,railway,roads,routes)):
        raise ValueError("Cross-source country identity mismatch")
    records,qc=build_city_pairs(plants,railway,roads,aliases,period)
    source=add_road_paths(city,COUNTRIES[city],railway,routes)
    joined=[]
    for r in records:
        key=tuple(sorted((r["a"],r["b"])))
        if key not in source:
            raise ValueError("Source-eligible botanical pair missing road graph")
        entry=source[key]
        if not entry["road_network_source_usable"]:
            continue
        joined.append({**r,
             "road_network_shortest_path_km":entry["network_length_km"],
             "road_network_circuity":entry["circuity"],
        })
    return joined,{
        **qc,
        "n_botanical_source_pairs_before_road_path":len(records),
        "n_botanical_source_pairs_with_road_route":len(joined),
        "n_source_pairs_road_route_not_usable":len(records)-len(joined),
        "n_all_36_city_road_pairs_with_route":sum(r["road_network_source_usable"] for r in source.values())
    }


def compare(rows,require_both_rail=False,calipers=(0.45,0.40)):
    """Greedy 1:1 without reuse, pre-outcome covariates only.

    Road route matching is conditional on both cells road-snappable,
    and old geography, plant source-record effort and road-density gates.
    It cannot correct unknown train operation, species sampling, climate,
    municipal freight access or omitted national rail paths.
    """
    if any(not isinstance(v,(int,float)) or v<=0 for v in calipers):
        raise ValueError("Invalid road-match calipers")
    subset=rows if not require_both_rail else [
        r for r in rows if r["both_connector_source_rail"]
    ]
    positive=sorted((r for r in subset if r["connected"]),
                    key=lambda r:(r["distance_km"],r["a"],r["b"]))
    negative=[r for r in subset if not r["connected"]]
    used=set()
    matches=[]
    for pos in positive:
        alternatives=[]
        for i,neg in enumerate(negative):
            if i in used:continue
            dd=abs(pos["distance_km"]-neg["distance_km"])
            max_distance_diff=max(5.,.33*max(pos["distance_km"],neg["distance_km"]))
            if dd>max_distance_diff:continue
            gdiff=abs(pos["log_mean_record_count"]-neg["log_mean_record_count"])
            ediff=abs(pos["log_record_effort_ratio"]-neg["log_record_effort_ratio"])
            if gdiff>1.0 or ediff>1.0:continue
            road_mean=abs(pos["road_mean_log_density"]-neg["road_mean_log_density"])
            road_balance=abs(pos["road_absolute_log_ratio"]-neg["road_absolute_log_ratio"])
            if road_mean>0.5 or road_balance>0.7:continue
            distance_log=abs(math.log(pos["road_network_shortest_path_km"]/
                                      neg["road_network_shortest_path_km"]))
            circuity_log=abs(math.log(pos["road_network_circuity"]/
                                      neg["road_network_circuity"]))
            if distance_log>calipers[0] or circuity_log>calipers[1]:
                continue
            score=dd/max_distance_diff+gdiff+ediff+road_mean+road_balance
            score+=distance_log/max(calipers[0],.01)+circuity_log/max(calipers[1],.01)
            alternatives.append((score,i))
        if alternatives:
            _,index=min(alternatives)
            used.add(index)
            ctrl=negative[index]
            matches.append(pos["source_rarefied20_introduced_similarity"]-
                           ctrl["source_rarefied20_introduced_similarity"])
    good=(len(positive)>=MIN_PER_CITY_MATCHES and
          len(negative)>=MIN_PER_CITY_MATCHES and
          len(matches)>=MIN_PER_CITY_MATCHES)
    return {
        "n_source_botanical_pairs_with_rail_connected":len(positive),
        "n_source_botanical_pairs_with_rail_unconnected":len(negative),
        "n_matched_exposed_control_pair_contrasts":len(matches),
        "city_level_median_similarity_delta":round(median(matches),6) if good else None,
        "city_comparison_supported":good,
        "uses_road_density_and_shortest_paths":True,
        "both_cells_need_mapped_rail":require_both_rail,
        "caliper_log_road_path":calipers[0],
        "caliper_log_road_circuity":calipers[1],
        "not_a_causal_train_seed_dispersal_estimate":True,
    }


def summarize_city_blocks(results,period,comparison_label,caliper_label):
    values=[]
    for name,v in results.items():
        per=v["periods"].get(period)
        if not per:continue
        x=per["comparisons"][comparison_label][caliper_label]
        if x["city_comparison_supported"]:
            values.append((name,x["city_level_median_similarity_delta"]))
    return {
        "n_informative_independent_city_blocks":len(values),
        "positive_blocks":sum(v>0 for _,v in values),
        "negative_blocks":sum(v<0 for _,v in values),
        "median_city_level_descriptive_delta":round(median(v for _,v in values),6) if values else None,
        "minimum_five_city_blocks":len(values)>=MIN_INDEPENDENT_CITIES,
        "source_city_deltas":{city:value for city,value in values},
        "not_causal_railway_effect_or_significance_test":True,
    }


def analyze(plant,rail,griis,transport,road_graph):
    if not(plant.get("source_complete") and rail.get("complete") and
           griis.get("complete") and transport.get("complete") and road_graph.get("complete")):
        raise ValueError("Incomplete fifth source evidence; do not fit railway result")
    for source in (plant,rail,transport,road_graph):
        if set(source.get("regions",{}))!=set(COUNTRIES):
            raise ValueError("Ten-city source population not complete")
    checked=valid_griis(griis)
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "study":"Posthoc five-source GRIIS-introduced observed composition vs mapped rail physical connector, conditional on physical ROAD routes",
        "source_scope":"one preselected city per country; NOT independently sampled national macro ecoregions",
        "country_introduced_source_keys":"nine GRIIS checklists; TW unassessed and KR key support inadequate",
        "regions":{},"block_summaries":{},
        "not_train_seed_transport_or_evolutionary_adaptation_effect":True,
        "source_results_exploratory_given_prior_v27_outcomes":True,
    }
    for city,country in COUNTRIES.items():
        out["regions"][city]={"country":country,"periods":{}}
        if checked.get(country) is None:
            out["regions"][city]["source_exclusion"]="No adequate GRIIS introduced taxon keys"
            continue
        for period in PERIODS:
            observations,qc=source_pair_records(plant,rail,griis,transport,road_graph,city,period)
            result={"recording_and_road_support":qc,"comparisons":{}}
            for tag,need in (("all",False),("both_rail",True)):
                result["comparisons"][tag]={
                    label:compare(observations,require_both_rail=need,calipers=caliper)
                    for label,caliper in SENSITIVITY_CALIPERS.items()
                }
            out["regions"][city]["periods"][period]=result
            x=result["comparisons"]["all"]["main"]
            print(json.dumps({
                "city":city,"period":period,
                "source_pairs_with_road_routes":qc["n_botanical_source_pairs_with_road_route"],
                "roadpath_matched":x["n_matched_exposed_control_pair_contrasts"],
                "median_diff":x["city_level_median_similarity_delta"],
                "strict_bothrail_median":result["comparisons"]["both_rail"]["main"]["city_level_median_similarity_delta"]
            },ensure_ascii=False),flush=True)
    for period in PERIODS:
        out["block_summaries"][period]={
            tag:{
                label:summarize_city_blocks(out["regions"],period,tag,label)
                for label in SENSITIVITY_CALIPERS
            } for tag in ("all","both_rail")
        }
    out["caveats"]=[
        "The actual ROAD graph paths are optimistic physical routes: no legal turn restrictions, motor-vehicle access, observed traffic, or freight routing.",
        "Current Overture roads and rails (2026) are not exposure histories for GBIF observation periods 2010–2017 and 2018–2025.",
        "Spatial graph routes are truncated by urban window boundary and use nearest connector to cell centroid rather than real botanical records.",
        "Only GRIIS-listed species are included; unlisted introduced/native status is unknown.",
        "Indexed occurrence rarefaction does not account for source event duplicates, observers, georeference errors, or habitat variation.",
        "Small city cells lack independent country×ecoregion block replication, including climate and port controls.",
        "Matched cell pairs share sites, no per-pair p-value, and comparison calipers were developed after earlier outcomes: all tested calipers must be reported.",
        "A zero or negative source association does not prove absence of actual train-mediated propagule dispersal."
    ]
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--gbif-source",required=True)
    a.add_argument("--rail-source",required=True)
    a.add_argument("--griis-source",required=True)
    a.add_argument("--transport-source",required=True)
    a.add_argument("--road-path-source",required=True)
    a.add_argument("--output",default="results/audit/griis_rail_road_paths_v210.json")
    p=a.parse_args()
    sources=[json.loads(Path(path).read_text(encoding="utf-8")) for path in (
        p.gbif_source,p.rail_source,p.griis_source,
        p.transport_source,p.road_path_source)]
    out=analyze(*sources)
    dest=Path(p.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"n_country_windows":len(out["regions"]),
        "block_summary":out["block_summaries"],
        "mechanism_unidentified":out["not_train_seed_transport_or_evolutionary_adaptation_effect"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
