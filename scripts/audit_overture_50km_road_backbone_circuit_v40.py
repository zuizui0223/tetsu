#!/usr/bin/env python3
"""v4.0 road NETWORK electrical control: MAJOR road backbone sensitivity.

Purpose: test if Circuit-style alternative physical paths are unique to
railways or also abundant in an independent ROAD network. True Overture
2026-10 connector junctions, source geometries and identical 50km
geographic cell anchors / 75km source halo are used.

To remain under GitHub Action RAM, use a fixed subset of motorable
backbone road classes motorway/trunk/primary/secondary (no tertiary,
residential, unclassified); this is a MAJOR-ROAD CIRCUIT SENSITIVITY,
NOT full ROAD-network Circuitscape. The full ROAD shortest-path source
used in v3.3 includes more classes and is NOT directly an equivalent
resistance denominator. This script recomputes its own backbone
Dijkstra on the exact graph used for circuit calculations.

Graph edge resistance is mapped road length in km-equivalent units,
not traffic flow, one-way, allowed freight vehicles or plant transport.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json,math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from statistics import median

from pyproj import Geod
from shapely import wkb

from audit_overture_50km_macro_rail_v30 import (
    BLOCKS,GEOD,macro_geometry,bbox_with_halo,open_network_connection
)
from audit_overture_50km_road_paths_v33 import (
    network_source as road_network_source, _dijkstra
)
from audit_overture_road_network_v28 import _valid_road_segment
from tetsu.circuit_graph import add_conductance,pairwise_effective_resistance

MAJOR_ROAD_CLASSES=("motorway","trunk","primary","secondary")
STRATEGIES=("nearest_in_cell","nearest_largest_component_in_cell")
MAX_SEGMENTS_PER_BLOCK=180000
MAX_REDUCED_SOLVER_NODES=35000
ROAD_RESISTANCE_MIN_KM=1e-6


def filtered_query(block,classes=MAJOR_ROAD_CLASSES,count_only=False):
    if not classes or any(x not in MAJOR_ROAD_CLASSES for x in classes):
        raise ValueError("Non-prespecified major-road source class")
    w,s,e,n=bbox_with_halo(block)
    selected=",".join("'"+s+"'" for s in classes)
    cols="COUNT(*)" if count_only else 'id, "class" AS road_class, connectors, ST_AsWKB(geometry) AS geometry_wkb'
    return f"""SELECT {cols} FROM read_parquet(
        's3://overturemaps-us-west-2/release/2026-10-21.0/theme=transportation/type=segment/*',
        hive_partitioning=1)
        WHERE subtype='road' AND "class" IN ({selected})
        AND bbox.xmin<={e} AND bbox.xmax>={w}
        AND bbox.ymin<={n} AND bbox.ymax>={s}""" + ("" if count_only else " ORDER BY id")


def major_road_circuit_graph(rows):
    graph={}
    ids=set()
    stats={"n_source_road_segments":0,"n_source_connector_links":0,
           "n_parallel_connector_pairs_in_source":0,
           "n_equal_at_connector_links":0}
    for raw in rows:
        sid,cls,geom,length,links=_valid_road_segment(raw)
        if cls not in MAJOR_ROAD_CLASSES:
            raise ValueError("Unexpected source road-class inclusion")
        if sid in ids:
            raise ValueError("Duplicated mapped road source ID")
        ids.add(sid)
        for (at,u),(bt,v) in zip(links,links[1:]):
            if u in graph and v in graph[u]:
                stats["n_parallel_connector_pairs_in_source"]+=1
            if at==bt:
                stats["n_equal_at_connector_links"]+=1
            resistance_km=max(ROAD_RESISTANCE_MIN_KM,length*(bt-at)/1000.)
            add_conductance(graph,u,v,1./resistance_km)
            stats["n_source_connector_links"]+=1
        stats["n_source_road_segments"]+=1
    if not graph:
        raise ValueError("No physical source road connectors to solve")
    stats["n_mapped_original_connector_nodes"]=len(graph)
    stats["n_distinct_undirected_connector_pairs"]=sum(len(x) for x in graph.values())//2
    return graph,stats


def analyze_block(block,rows,max_nodes=MAX_REDUCED_SOLVER_NODES):
    if block not in BLOCKS:
        raise ValueError("Region is not fixed in macro source frame")
    shortest,anchors,short_summary=road_network_source(block,rows)
    circuit,source=major_road_circuit_graph(rows)
    if set(shortest)!=set(circuit):
        raise ValueError("Road shortest and circuit source IDs differ")
    terminals=sorted({v["node"] for m in anchors.values() for v in m.values()})
    if len(terminals)<2:
        raise ValueError("Too few source road cell anchors")
    resistance,solver=pairwise_effective_resistance(
        circuit,terminals,max_reduced_nodes=max_nodes
    )
    out={
        "country":BLOCKS[block][0],
        "source_classes":list(MAJOR_ROAD_CLASSES),
        "road_backbone_not_full_network":True,
        "source_graph":source,"anchor_source":short_summary,"sparse_solver":solver,
        "strategies":{}
    }
    all_cells=[f"{block}_x{i}_y{j}" for j in range(4) for i in range(4)]
    for strategy in STRATEGIES:
        src=anchors[strategy]
        distance={}
        for a,t in src.items():
            targets=[v["node"] for name,v in src.items()
                     if name!=a and v["component"]==t["component"]]
            if targets:distance[a]=_dijkstra(shortest,t["node"],targets)
        pairs=[];qvals=[];ratio=[]
        for a,b in combinations(sorted(all_cells),2):
            aa=src.get(a);bb=src.get(b)
            item={
                "a":a,"b":b,
                "both_have_backbone_road_anchor":bool(aa and bb),
                "physically_connected_in_source_halo":False,
                "source_backbone_shortest_km":None,
                "source_backbone_R_eff_equiv_km":None,
                "source_backbone_Q_one_minus_Reff_over_shortest":None,
                "source_centroid_access_augmented_Q":None
            }
            if aa and bb and aa["component"]==bb["component"]:
                d=distance.get(a,{}).get(bb["node"])
                if aa["node"]==bb["node"]:d=0.
                k=tuple(sorted((aa["node"],bb["node"])))
                r=resistance.get(k) if k[0]!=k[1] else 0.
                if r is None or d is None:
                    raise ValueError("Resistor and shortest major-road source path connectivity disagree")
                min_km=d/1000.
                if r < -1e-8 or r > min_km+max(1e-7,min_km*1e-5):
                    raise ValueError("Road resistance violates Dijkstra path upper bound")
                item["physically_connected_in_source_halo"]=True
                item["source_backbone_shortest_km"]=round(min_km,6)
                item["source_backbone_R_eff_equiv_km"]=round(r,6)
                if min_km>1e-8:
                    q=max(0.,min(1.,1.-r/min_km))
                    item["source_backbone_Q_one_minus_Reff_over_shortest"]=round(q,7)
                    qvals.append(q);ratio.append(1.-q)
                access=aa["snap_m"]+bb["snap_m"]
                if d+access>0:
                    item["source_centroid_access_augmented_Q"]=round(
                        1.-(r*1000.+access)/(d+access),7)
            pairs.append(item)
        out["strategies"][strategy]={
            "n_50km_source_cells":16,
            "n_backbone_road_anchored_cells":len(src),
            "n_connected_source_pairs":sum(p["physically_connected_in_source_halo"] for p in pairs),
            "n_nonzero_distance_source_pairs":len(qvals),
            "median_source_Q_redundancy":round(median(qvals),6) if qvals else None,
            "median_source_R_over_shortest":round(median(ratio),6) if ratio else None,
            "n_pairs_with_5pct_redundancy":sum(x>=.05 for x in qvals),
            "ratio_result_is_backbone_only_not_full_road":True,
            "pair_source_diagnostics":pairs,
        }
    return out


def run(blocks):
    if not blocks or len(blocks)!=len(set(blocks)) or any(b not in BLOCKS for b in blocks):
        raise ValueError("Invalid outcome-independent country macro block list")
    report={
        "utc":datetime.now(timezone.utc).isoformat(),
        "study":"Fixed major-road 4-class Overture physical backbone circuit control",
        "source_classes":list(MAJOR_ROAD_CLASSES),
        "complete":False,"regions":{},"errors":{},
        "no_railway_or_road_plant_effect_measured":True,
    }
    db=open_network_connection()
    db.execute("INSTALL spatial; LOAD spatial")
    try:
        for b in blocks:
            try:
                n=int(db.execute(filtered_query(b,count_only=True)).fetchone()[0])
                if n<2:
                    raise ValueError("No mapped source major roads")
                if n>MAX_SEGMENTS_PER_BLOCK:
                    raise ValueError(f"Source size exceeds frozen solver budget: {n}>{MAX_SEGMENTS_PER_BLOCK}; use separate scalable solver, no fake circuit")
                data=db.execute(filtered_query(b)).fetchall()
                if len(data)!=n:
                    raise ValueError("Preflight major-road count disagrees with downloaded rows")
                source=analyze_block(b,data)
                report["regions"][b]=source
                part=source["strategies"]["nearest_in_cell"]
                print(json.dumps({
                    "block":b,"source_segments":n,
                    "n_road_backbone_anchors":part["n_backbone_road_anchored_cells"],
                    "source_circuit_pair_count":part["n_connected_source_pairs"],
                    "median_backbone_Q":part["median_source_Q_redundancy"],
                    "solver_max_reduced_nodes":source["sparse_solver"]["n_reduced_nodes_max"],
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                report["errors"][b]=f"{type(exc).__name__}: {str(exc)[:240]}"
                print(json.dumps({"source_error":b,"reason":report["errors"][b]},ensure_ascii=False),flush=True)
    finally:
        db.close()
    report["complete"]=not report["errors"] and len(report["regions"])==len(blocks)
    report["limitations"]=[
        "Major roads comprise only motorway/trunk/primary/secondary classes, not full mapped motorable road graph. No direct unqualified comparison to the full road shortest routes.",
        "Source road physical connectivity is NOT legal vehicle journey, freight traffic, shipment volume or propagule transport.",
        "Electrical resistance proportional to physical road km is a hypothetical resistor source, as with rail model.",
        "The 75km map halo may truncate longer alternative roads, raising effective resistance.",
        "A sparse solver input graph above frozen source/reduced node budgets is excluded, not silently approximated.",
        "No introduced plant response or actual active train operations have been used in this source-only control."
    ]
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--blocks",nargs="+",choices=list(BLOCKS),default=["us_midwest"])
    p.add_argument("--output",default="results/audit/road_backbone_circuit_v40.json")
    a=p.parse_args()
    r=run(a.blocks)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":r["complete"],"n_blocks":len(r["regions"]),"errors":r["errors"]},ensure_ascii=False),flush=True)
    if not r["complete"]:
        raise SystemExit("Real major road circuit control unavailable; no false control estimate")
if __name__=="__main__":
    main()
