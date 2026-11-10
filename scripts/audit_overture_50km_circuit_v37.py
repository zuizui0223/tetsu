#!/usr/bin/env python3
"""v3.7 test Circuitscape NETWORK theory against actual mapped RAIL paths.

Source rail 2026-10 Overture physical connector graph is reconstructed
from original source segment records; parallel routes are combined by
SUM OF CONDUCTANCES, not their shortest length. Normalized edge
resistance = geometric length in kilometers, not calibrated biological
seed dispersal resistance. This is a SciPy exact sparse Laplacian solve
of the network Circuitscape equation, *NOT* a run of Julia Circuitscape.

The PURPOSE is to see whether all-paths effective resistance differs
meaningfully from shortest-path distance in true pre-frozen 50km source
blocks. Two anchor strategies from v3.2, NO plant outcome input.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime,timezone
from itertools import combinations
import json
import math
from pathlib import Path
import sys
from statistics import median

from shapely import wkb

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.circuit_graph import (
    add_conductance,components,pairwise_effective_resistance,
    write_circuitscape_network
)
from tetsu.rail_graph import validate_segment
from audit_overture_50km_macro_rail_v30 import (
    BLOCKS,RAIL_CLASSES,GEOD,query_macro_rail,open_network_connection
)
from audit_overture_50km_rail_paths_v32 import source_network,path_m

STRATEGIES=("nearest_in_cell","nearest_largest_component_in_cell")
MIN_RESISTANCE_KM=1e-6  # v3.2 source epsilon 0.001m -> km
MIN_REDUNDANCY_MEANINGFUL=0.05
MIN_REDUNDANCY_TINY=0.01


def load_physical_multigraph_conductance(rows):
    """Return source network with true parallel tracks retained as conductance."""
    graph={}
    seen_ids=set()
    stats={"n_source_segments":0,"n_source_links":0,
           "n_parallel_connector_links":0,
           "n_same_at_connector_links":0,
           "source_resistance_per_km_default":1.}
    for sid,source_class,connectors,source_wkb in rows:
        if not isinstance(sid,str) or not sid or sid in seen_ids:
            raise ValueError("Missing/duplicate Overture mapped segment ID")
        seen_ids.add(sid)
        if source_class not in RAIL_CLASSES:
            raise ValueError("Nonapproved rail class in physical resistor source")
        shape=wkb.loads(bytes(source_wkb))
        if shape.geom_type!="LineString" or shape.is_empty or not shape.is_valid:
            raise ValueError("Invalid source rail line")
        length_m=abs(float(GEOD.geometry_length(shape)))
        record=validate_segment({
            "id":sid,"subtype":"rail","class":source_class,
            "length_m":length_m,"connectors":connectors
        },RAIL_CLASSES)
        if record is None:
            raise ValueError("Non-approved source rail segment was silently included")
        for (a,u),(b,v) in zip(record["connectors"],record["connectors"][1:]):
            if u in graph and v in graph[u]:
                stats["n_parallel_connector_links"]+=1
            if a==b:stats["n_same_at_connector_links"]+=1
            resistance_km=max(MIN_RESISTANCE_KM,(b-a)*length_m/1000)
            add_conductance(graph,u,v,1./resistance_km)
            stats["n_source_links"]+=1
        stats["n_source_segments"]+=1
    if not graph:raise ValueError("Zero true physical source rail connector links")
    comp=components(graph)
    n=len(graph)
    E=stats["n_source_links"]
    stats.update({
        "n_original_graph_nodes":n,
        "n_parallel_aggregated_undirected_connector_pairs":sum(len(a) for a in graph.values())//2,
        "n_components":len(comp),
        "multigraph_cycle_rank":E-n+len(comp),
        "source_parallel_linkages_not_dropped_for_circuit":True,
        "not_observed_train_service_or_freight_flow":True,
    })
    if stats["multigraph_cycle_rank"]<0:
        raise ValueError("Physical multigraph cycle rank impossible")
    return graph,stats


def diagnostic(block,rows,export_dir=None,max_reduced_nodes=35000):
    # v3.2 provides the EXISTING outcome-blind location anchors and
    # SOURCE-Dijkstra's shortest route. Same source, same geographic frame.
    min_graph,original_components,anchors,cell_src,graph_source=source_network(block,rows)
    conductance,source_stats=load_physical_multigraph_conductance(rows)
    if set(min_graph)!=set(conductance):
        raise ValueError("Dijkstra and all-paths circuit source connector sets differ")
    terminals=sorted({v["node"] for a in anchors.values() for v in a.values()})
    if len(terminals)<2:
        raise ValueError("Too few source-identified 50km macro terminal connectors")
    exact,solver=pairwise_effective_resistance(conductance,terminals,max_reduced_nodes)
    output={
        "country":BLOCKS[block][0],
        "block":block,
        "model_name":"Circuitscape-network-equivalent sparse Laplacian / source physical rail length resistances",
        "actual_Julia_Circuitscape_executed":False,
        "resistance_unit":"hypothetical resistance units = geodesic km of rail source edge; NOT measured seed transmission",
        "source_graph":source_stats,
        "original_v32_anchor_source":graph_source,
        "sparse_solver":solver,
        "strategies":{}
    }
    if export_dir is not None:
        # Local opt-in, may be restricted by source licensing. NOT CI/uploaded.
        output["local_Circuitscape_network_files"]=write_circuitscape_network(
            conductance,terminals,Path(export_dir)/block
        )
    for strategy in STRATEGIES:
        src=anchors[strategy]
        rows_out=[]
        dcache={}
        for a,source in src.items():
            targets=[v["node"] for b,v in src.items() if b!=a and
                     v["component"]==source["component"]]
            if targets:
                dcache[a]=path_m(min_graph,source["node"],targets)
        all_cells=sorted(cell_src)
        for a,b in combinations(all_cells,2):
            p=src.get(a);q=src.get(b)
            entry={
                "a":a,"b":b,
                "both_source_rail_centroid_anchors":bool(p and q),
                "physical_circuit_reachable":False,
                "R_effective_rail_equiv_km":None,
                "shortest_rail_connector_km":None,
                "R_effective_divided_by_Dijkstra":None,
                "independent_path_redundancy_1minusRoverD":None,
                "centroid_augmented_R_eff_divided_by_Dijkstra":None,
            }
            if p and q and p["component"]==q["component"]:
                v=exact.get(tuple(sorted([p["node"],q["node"]])))
                d=dcache.get(a,{}).get(q["node"])
                if p["node"]==q["node"]:
                    v=0.;d=0.
                if v is None or d is None:
                    raise ValueError("Circuit solver and shortest-rail path disagree on connectivity")
                # v3.2 Dijkstra returns METERS; circuit edge resistances
                # are normalized in KM-equivalent units. Always convert
                # before comparing! Otherwise a fabricated 1000x circuit
                # redundancy appears in all source regions.
                shortest_km=d/1000.
                if v< -1e-7 or d<0 or v>shortest_km+max(1e-7,shortest_km*1e-5):
                    raise ValueError(f"Effective resistance violates shortest-path bound: {v} > {shortest_km}")
                entry["physical_circuit_reachable"]=True
                entry["R_effective_rail_equiv_km"]=round(v,6)
                entry["shortest_rail_connector_km"]=round(shortest_km,6)
                if shortest_km>1e-8:
                    ratio=min(1.,max(0.,v/shortest_km))
                    entry["R_effective_divided_by_Dijkstra"]=round(ratio,8)
                    entry["independent_path_redundancy_1minusRoverD"]=round(1.-ratio,8)
                total_access=p["snap_m"]+q["snap_m"]
                entry["centroid_augmented_R_eff_divided_by_Dijkstra"]=round(
                    (v*1000+total_access)/(d+total_access),8
                ) if d+total_access>0 else None
            rows_out.append(entry)
        values=[v["independent_path_redundancy_1minusRoverD"] for v in rows_out
                if v["independent_path_redundancy_1minusRoverD"] is not None]
        ratios=[v["R_effective_divided_by_Dijkstra"] for v in rows_out
                if v["R_effective_divided_by_Dijkstra"] is not None]
        output["strategies"][strategy]={
            "n_50km_source_cells":16,
            "n_source_rail_anchored_cells":len(src),
            "n_source_cell_pairs":120,
            "n_circuit_reachable_source_pairs":sum(v["physical_circuit_reachable"] for v in rows_out),
            "n_nonzero_shortest_path_terminal_pairs":len(values),
            "n_pairs_with_at_least_1pct_independent_path_redundancy":sum(x>=MIN_REDUNDANCY_TINY for x in values),
            "n_pairs_with_at_least_5pct_independent_path_redundancy":sum(x>=MIN_REDUNDANCY_MEANINGFUL for x in values),
            "fraction_pairs_with_5pct_redundancy":round(sum(x>=MIN_REDUNDANCY_MEANINGFUL for x in values)/len(values),6) if values else None,
            "median_redundancy":round(median(values),6) if values else None,
            "median_R_effective_over_shortest":round(median(ratios),6) if ratios else None,
            "min_R_effective_over_shortest":round(min(ratios),6) if ratios else None,
            "n_R_effective_gt_shortest_bound_violations":0,
            "additive_information_over_shortest_path_unknown_without_external_validation":True,
            "pairwise_source_diagnostics":rows_out,
        }
    output["limits"]=[
        "All edge resistances assumed proportional to mapped rail length in km: there is no validated railway seed-propagule permeability model.",
        "Multiple track links sharing endpoint IDs have SUMMED conductance; shortest paths retain only minimum segment length.",
        "Connectivity and effective resistance are symmetric physical network models; train schedules, direction, transfer, freight carriage, herbicide, trackside seed establishment unknown.",
        "A sparse grounded-Laplacian exact numeric solve is circuit-theoretic and compatible with Circuitscape NETWORK math, but this run is NOT verified Circuitscape.jl output.",
        "Rarefied species, rail-vs-road ecological coefficient, climate/observer sampling, country ecoregions are not analyzed in this source-only step.",
        "75km source halo can exclude paths, making local graph effective resistance larger than real whole-network potential.",
        "Near-tree networks give effective resistance near shortest-path distance; before ecological fitting, quantify redundancy and reject constant explanatory variables."
    ]
    return output


def run(blocks,export_dir=None):
    if not blocks or any(x not in BLOCKS for x in blocks) or len(set(blocks))!=len(blocks):
        raise ValueError("Geographical source blocks must use fixed 50km country frame")
    out={"utc":datetime.now(timezone.utc).isoformat(),
         "source":"2026-10 Overture real physical railway mapped source connector network",
         "model":"Mathematically exact sparse resistor network for 50km source-referenced anchor pairs; Circuitscape-compatible but Julia not run",
         "regions":{},"errors":{},
         "not_train_dispersal_result":True}
    db=open_network_connection()
    db.execute("INSTALL spatial; LOAD spatial;")
    try:
        for block in blocks:
            try:
                rows=db.execute(query_macro_rail(block)).fetchall()
                result=diagnostic(block,rows,export_dir=export_dir)
                out["regions"][block]=result
                a=result["strategies"]["nearest_in_cell"]
                print(json.dumps({
                    "block":block,
                    "source_rail_segments":result["source_graph"]["n_source_segments"],
                    "source_cycle_rank":result["source_graph"]["multigraph_cycle_rank"],
                    "n_original_connectors":result["source_graph"]["n_original_graph_nodes"],
                    "n_solver_reduced_nodes_max":result["sparse_solver"]["n_reduced_nodes_max"],
                    "n_rail_reachable_pairs":a["n_circuit_reachable_source_pairs"],
                    "n_5pct_redundant_pairs":a["n_pairs_with_at_least_5pct_independent_path_redundancy"],
                    "median_1_minus_RoverD":a["median_redundancy"],
                    "median_RoverD":a["median_R_effective_over_shortest"],
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                out["errors"][block]=f"{type(exc).__name__}: {str(exc)[:330]}"
                print(json.dumps({"source_error":block,"detail":out["errors"][block]},ensure_ascii=False),flush=True)
    finally:
        db.close()
    out["complete"]=not out["errors"] and len(out["regions"])==len(blocks)
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--blocks",nargs="+",choices=list(BLOCKS),default=["us_midwest","fr_burgundy"])
    a.add_argument("--output",default="results/audit/50km_rail_circuit_resistance_v37.json")
    a.add_argument("--local-circuitscape-export-dir",default=None)
    args=a.parse_args()
    res=run(args.blocks,args.local_circuitscape_export_dir)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":res["complete"],"n_blocks":len(res["regions"]),"errors":res["errors"]},ensure_ascii=False),flush=True)
    if not res["complete"]:
        raise SystemExit("Circuit resistance source quality gate failed; never infer seeded rail movement")


if __name__=="__main__":
    main()
