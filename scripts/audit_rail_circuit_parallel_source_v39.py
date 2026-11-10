#!/usr/bin/env python3
"""v3.9 sceptical Circuitscape railway redundancy source mapping audit.

Overture 2026 physical mapped rail segments, pre-frozen JP/US/FR
50km equal-area blocks. Recompute pairwise effective resistance for
two *source data* assumptions, both on identical terminals/graphs:

A) SUM all parallel source link conductances (v3.7 baseline)
B) SHORTEST_ONLY retain one source edge between identical connector IDs
   (conservative possible-map-duplicate sensitivity; different
   physical routes/junctions are still retained).

The second model is not a claim that real parallel tracks are duplicates.
For connected anchors, R_sum <= R_shortest_only <= Dijkstra_shortest.

No flora, species, freight, observer data or causal effect is used.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json,math
from pathlib import Path
from statistics import median

from tetsu.circuit_graph import pairwise_effective_resistance
from audit_overture_50km_macro_rail_v30 import (
    BLOCKS,query_macro_rail,open_network_connection
)
from audit_overture_50km_rail_paths_v32 import source_network,path_m
from audit_overture_50km_circuit_v37 import (
    STRATEGIES,load_physical_multigraph_conductance
)

POLICIES=("sum","shortest_only")
MIN_BIOLOGICAL_CLAIM=False


def audit_block(block,rows,max_reduced_nodes=35000):
    shortest,comp,anchors,cell_summary,meta=source_network(block,rows)
    all_terminals=sorted({x["node"] for strat in STRATEGIES
                          for x in anchors[strat].values()})
    if len(all_terminals)<2:
        raise ValueError("Insufficient source anchor nodes")
    circuits={};stats={};graph_meta={}
    for policy in POLICIES:
        g,source=load_physical_multigraph_conductance(rows,parallel_policy=policy)
        if set(shortest)!=set(g):
            raise ValueError("Changing parallel policy incorrectly removed source vertices")
        values,details=pairwise_effective_resistance(
            g,all_terminals,max_reduced_nodes=max_reduced_nodes)
        circuits[policy]=values;stats[policy]=details;graph_meta[policy]=source
    results={}
    for strategy in STRATEGIES:
        chosen=anchors[strategy]
        distances={}
        for a,source in chosen.items():
            lookup=[p["node"] for b,p in chosen.items()
                    if a!=b and p["component"]==source["component"]]
            if lookup:distances[a]=path_m(shortest,source["node"],lookup)
        delta=[];sum_redundancy=[];one_redundancy=[];rows_out=[]
        for a,b in combinations(sorted(cell_summary),2):
            ia=chosen.get(a);ib=chosen.get(b)
            row={
                "a":a,"b":b,
                "network_source_connected":bool(ia and ib and ia["component"]==ib["component"]),
                "shortest_rail_km":None,
                "R_sum_over_shortest":None,
                "R_single_link_over_shortest":None,
                "Q_all_parallel_links":None,
                "Q_single_link_per_connector_pair":None,
                "parallel_source_edge_Q_increment":None
            }
            if row["network_source_connected"]:
                d=(distances.get(a,{}) or {}).get(ib["node"])
                if ia["node"]==ib["node"]:d=0.
                key=tuple(sorted([ia["node"],ib["node"]]))
                r_sum=circuits["sum"].get(key) if key[0]!=key[1] else 0.
                r_one=circuits["shortest_only"].get(key) if key[0]!=key[1] else 0.
                if d is None or r_sum is None or r_one is None:
                    raise ValueError("Connectivity does not match source graph paths")
                distance=d/1000.
                if not (distance>=0 and math.isfinite(distance)
                        and -1e-6<=r_sum<=r_one+1e-6
                        and r_one<=distance+max(1e-7,distance*1e-5)):
                    raise ValueError("Source physical resistor ordering violated")
                row["shortest_rail_km"]=round(distance,6)
                if distance>1e-8:
                    ratio_sum=min(1.,max(0.,r_sum/distance))
                    ratio_one=min(1.,max(0.,r_one/distance))
                    q_sum=1.-ratio_sum;q_one=1.-ratio_one
                    row.update({
                        "R_sum_over_shortest":round(ratio_sum,7),
                        "R_single_link_over_shortest":round(ratio_one,7),
                        "Q_all_parallel_links":round(q_sum,7),
                        "Q_single_link_per_connector_pair":round(q_one,7),
                        "parallel_source_edge_Q_increment":round(q_sum-q_one,7)
                    })
                    sum_redundancy.append(q_sum)
                    one_redundancy.append(q_one)
                    delta.append(q_sum-q_one)
            rows_out.append(row)
        if not sum_redundancy:raise ValueError("No circuit-supported source cell pairs")
        results[strategy]={
            "n_source_50km_cells":16,
            "n_source_focal_rail_cells":len(chosen),
            "n_valid_nonzero_connector_paths":len(delta),
            "median_Q_with_parallel_sum":round(median(sum_redundancy),6),
            "median_Q_with_parallel_collapsed":round(median(one_redundancy),6),
            "median_Q_attributable_to_parallel_source_links":round(median(delta),6),
            "n_pairs_Q_collapsed_at_least_5pct":sum(x>=.05 for x in one_redundancy),
            "n_pairs_Q_parallel_sum_at_least_5pct":sum(x>=.05 for x in sum_redundancy),
            "source_parallel_labels_sensitivity_not_independent_train_paths":True,
            "pair_diagnostics":rows_out
        }
    return {
        "country":BLOCKS[block][0],
        "source_graph_policies":graph_meta,
        "sparse_solver_summaries":stats,
        "strategies":results,
        "no_plant_response_analyzed":True,
        "no_rail_operation_or_cargo_flow_inferred":True,
    }


def run(blocks):
    if not blocks or len(set(blocks))!=len(blocks) or any(b not in BLOCKS for b in blocks):
        raise ValueError("Source-frozen block frame required")
    db=open_network_connection()
    db.execute("INSTALL spatial; LOAD spatial")
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "question":"Does Circuitscape-equivalent rail redundancy persist without parallel source connector labels?",
        "regions":{},"errors":{},
        "not_railway_seed_dispersal_result":True,
    }
    try:
        for b in blocks:
            try:
                rows=db.execute(query_macro_rail(b)).fetchall()
                v=audit_block(b,rows)
                out["regions"][b]=v
                chosen=v["strategies"]["nearest_in_cell"]
                links=v["source_graph_policies"]["sum"]
                print(json.dumps({
                    "block":b,
                    "source_parallel_connector_links":links["n_parallel_connector_links"],
                    "n_reachable_cell_pairs":chosen["n_valid_nonzero_connector_paths"],
                    "Q_sum_median":chosen["median_Q_with_parallel_sum"],
                    "Q_dedup_median":chosen["median_Q_with_parallel_collapsed"],
                    "Q_parallel_increment_median":chosen["median_Q_attributable_to_parallel_source_links"],
                    "n_redundant_if_dedup":chosen["n_pairs_Q_collapsed_at_least_5pct"],
                },ensure_ascii=False),flush=True)
            except Exception as e:
                out["errors"][b]=f"{type(e).__name__}: {str(e)[:260]}"
                print(json.dumps({"source_error":b,"error":out["errors"][b]},ensure_ascii=False),flush=True)
    finally:
        db.close()
    out["complete"]=len(out["regions"])==len(blocks) and not out["errors"]
    out["interpretation_limits"]=[
        "Shortest-only parallel-policy sensitivity tests mapping/duplicated label risk, not confirmed physical one-track operation.",
        "All source resistor parameters remain hypothetical track length, not actual seed movement or train service.",
        "Real distinct alternative routes with different source junction IDs remain present under either policy.",
        "Only within 75km halo mapped graph, not national schedule or freight routes.",
        "There are NO vegetation/ecological model results in this geometry audit.",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--blocks",nargs="+",choices=list(BLOCKS),
                   default=["jp_tohoku","us_midwest","fr_burgundy"])
    p.add_argument("--output",default="results/audit/rail_circuit_parallel_source_sensitivity_v39.json")
    x=p.parse_args()
    out=run(x.blocks)
    dest=Path(x.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"n_blocks":len(out["regions"]),
                      "source_errors":out["errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Parallel rail source audit incomplete; no circuit inference")
if __name__=="__main__":
    main()
