#!/usr/bin/env python3
"""v3.8: test circuit-theory redundancy against observed alien plant similarity.

The independent Overture source circuit network v3.7 is already computed
BEFORE any biological response is used by this script. Source-verified
country-GRIIS speciesKey × 50km cell × time, Overture rail Dijkstra,
Overture road shortest paths and rail effective resistance share the
SAME predeclared geographical pair IDs (US Midwest and FR Burgundy).

This is a descriptive *posthoc* incremental diagnostic, NOT an
independently preregistered ecological hypothesis, no p-value,
no national independent replication (only 2 macro blocks).
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json
import math
from pathlib import Path

from analyze_macro_50km_rail_road_aliens_v34 import (
    MACRO_BLOCKS, STRATEGIES, partial_rank_r, _eligible_both_period_cells,
    _routes, _quartiles
)
from analyze_macro_rail_assemblages_v22 import expected_rarefied_overlap

PERIODS=("2010_2017","2018_2025")
MIN_PAIR_ROWS=20


def circuit_pairs(circuit,block,strategy):
    region=circuit["regions"][block]
    if region.get("actual_Julia_Circuitscape_executed") is not False:
        raise ValueError("Source must declare whether Julia Circuitscape was actually run")
    entries=region["strategies"][strategy]["pairwise_source_diagnostics"]
    if len(entries)!=120:
        raise ValueError("Incomplete pairwise circuit resistor network source")
    out={}
    for p in entries:
        a,b=p["a"],p["b"]
        k=tuple(sorted((a,b)))
        if k in out or a==b:
            raise ValueError("Duplicate circuit physical source pair")
        out[k]=p
    return out


def analyze_city(block,flora,rail,road,circuit,allowed_blocks=MACRO_BLOCKS):
    eligible=_eligible_both_period_cells(block,flora)
    if block not in allowed_blocks or len(eligible)<5:
        raise ValueError("Insufficient country-GRIIS prequalified source cells")
    result={
        "n_two_period_eligible_GRIIS_50km_cells":len(eligible),
        "periods":{},
    }
    for period in PERIODS:
        result["periods"][period]={}
        for strategy in STRATEGIES:
            shortest=_routes(rail,block,strategy,"pair_routes")
            roads=_routes(road,block,strategy,"pair_paths")
            circuits=circuit_pairs(circuit,block,strategy)
            if set(shortest)!=set(roads) or set(roads)!=set(circuits):
                raise ValueError("Source pair universe for road/rail/circuit mismatches")
            rows=[]
            for a,b in combinations(sorted(eligible),2):
                key=(a,b)
                r=shortest[key];d=roads[key];c=circuits[key]
                q=c.get("independent_path_redundancy_1minusRoverD")
                circuit_ratio=c.get("centroid_augmented_R_eff_divided_by_Dijkstra")
                railkm=r.get("cell_centroid_snap_augmented_rail_path_km")
                roadkm=d.get("centroid_snap_augmented_road_path_km")
                geokm=r.get("source_cell_centroid_geodesic_km")
                if not all(isinstance(v,(float,int)) and not isinstance(v,bool)
                           and math.isfinite(v) for v in (q,circuit_ratio,railkm,roadkm,geokm)):
                    continue
                if not (0<=q<1.00001 and 0<circuit_ratio<=1.00001 and
                        railkm>0 and roadkm>0 and geokm>0):
                    raise ValueError("Invalid physical circuit source or route denominator")
                if r.get("rail_to_geodesic_ratio") is None or d.get("road_to_geodesic_ratio") is None:
                    continue
                if r["rail_to_geodesic_ratio"]<.95 or d["road_to_geodesic_ratio"]<.95:
                    continue
                e1=eligible[a][period];e2=eligible[b][period]
                obs1=set(e1["counts"]);obs2=set(e2["counts"])
                rawJ=len(obs1&obs2)/len(obs1|obs2)
                rare=expected_rarefied_overlap(e1["rarified_probs"],e2["rarified_probs"])
                if rare is None:
                    continue
                n1=e1["all_plant_indexed_records"];n2=e2["all_plant_indexed_records"]
                rows.append({
                    "x_log_short_rail_over_road":math.log(railkm/roadkm),
                    "x_log_circuit_rail_over_road":math.log((railkm*circuit_ratio)/roadkm),
                    "x_redundancy_Q":q,
                    "log_geo":math.log(geokm),
                    "log_road":math.log(roadkm),
                    "log_short_rail":math.log(railkm),
                    "log_count_ratio":abs(math.log(n1/n2)),
                    "log_count_geomean":(math.log(n1)+math.log(n2))/2,
                    "rarefied20_source_similarity":rare,
                    "raw_observed_key_jaccard":rawJ,
                })
            summary={
                "n_dependent_cell_pair_rows":len(rows),
                "n_two_period_qualified_cell_sources":len(eligible),
                "Q_alternative_path_fraction_quartiles":_quartiles([r["x_redundancy_Q"] for r in rows]),
                "n_Q_at_least_5pct":sum(r["x_redundancy_Q"]>=.05 for r in rows),
                "effect_of_Q_as_additional_predictor_p_value":None,
                "responses":{}
            }
            for response in ("rarefied20_source_similarity","raw_observed_key_jaccard"):
                Y=[r[response] for r in rows]
                base=[[r["log_geo"],r["log_road"],r["log_count_ratio"],
                       r["log_count_geomean"]] for r in rows]
                covQ=[b+[r["log_short_rail"]] for b,r in zip(base,rows)]
                summary["responses"][response]={
                    "partial_rank_shortest_rail_vs_road":partial_rank_r(
                        [r["x_log_short_rail_over_road"] for r in rows],Y,base),
                    "partial_rank_Circuitscape_equivalent_effective_rail_vs_road":partial_rank_r(
                        [r["x_log_circuit_rail_over_road"] for r in rows],Y,base),
                    "partial_rank_extra_path_redundancy_Q_at_fixed_shortest_rail":partial_rank_r(
                        [r["x_redundancy_Q"] for r in rows],Y,covQ),
                    "source_railway_propagule_flux_unidentified":True,
                }
            result["periods"][period][strategy]=summary
            print(json.dumps({
                "block":block,"period":period,"anchor":strategy,
                "n_source_pairs":len(rows),
                "Q_quartiles":summary["Q_alternative_path_fraction_quartiles"],
                "rho_shortest":summary["responses"]["rarefied20_source_similarity"]["partial_rank_shortest_rail_vs_road"],
                "rho_effective":summary["responses"]["rarefied20_source_similarity"]["partial_rank_Circuitscape_equivalent_effective_rail_vs_road"],
                "rho_Q_given_shortest":summary["responses"]["rarefied20_source_similarity"]["partial_rank_extra_path_redundancy_Q_at_fixed_shortest_rail"]
            },ensure_ascii=False),flush=True)
    return result


def run(flora,rail,road,circuit):
    if any(not x.get("complete") for x in (flora,rail,road,circuit)):
        raise ValueError("Missing/incomplete real source; cannot infer biological railway association")
    if set(road.get("regions",{}))!=set(MACRO_BLOCKS) or not set(MACRO_BLOCKS).issubset(flora.get("regions",{})):
        raise ValueError("Independent macro region sample has changed")
    if not set(MACRO_BLOCKS).issubset(rail.get("regions",{})) or not set(MACRO_BLOCKS).issubset(circuit.get("regions",{})):
        raise ValueError("Required frozen source network region missing")
    result={
        "utc":datetime.now(timezone.utc).isoformat(),
        "question":"Does Circuitscape-network extra route redundancy add noncausal association with observed introduced plant species beyond shortest railway path?",
        "model_definition":"Source circuit resistance = length_km physical resistor network; not verified seed carrying capacity",
        "only_two_independent_country_macro_blocks":True,
        "exploratory_posthoc_model_not_preregistered":True,
        "actual_Julia_Circuitscape_was_not_run":True,
        "regions":{},
    }
    for block in MACRO_BLOCKS:
        result["regions"][block]=analyze_city(block,flora,rail,road,circuit)
    result["interpretation_limits"]=[
        "Pairwise plant similarities are not independent; this descriptive partial-rank analysis deliberately reports NO P-VALUE.",
        "Resistor parameters are hypothetical railway lengths; not calibrated train-seed vector use, track activity or freight traffic.",
        "Sources are 2026-10 physical map graphs overlaid with 2010-25 historical GBIF observation records, not true before-after rail exposure.",
        "Country GRIIS status is an incomplete regional checklist. Nonmembers remain UNKNOWN, not native.",
        "Only two preselected independent macro blocks US Midwest and France Burgundy; may not generalize to 8-country recorded flora.",
        "GBIF speciesKey facet presence is recording opportunity, not verified ecological absence or systematic plot incidence.",
        "Resistor effective resistance and shortest path are topologically related; higher partial correlation does not establish prediction improvement.",
        "Current measurements are not Julia Circuitscape output; exact sparse electrical network solver and Circuitscape network input adapter are separately verified."
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--griis-source",required=True)
    p.add_argument("--rail-shortest-source",required=True)
    p.add_argument("--road-shortest-source",required=True)
    p.add_argument("--circuit-source",required=True)
    p.add_argument("--output",default="results/audit/50km_alien_circuit_incremental_v38.json")
    a=p.parse_args()
    z=[json.loads(Path(x).read_text(encoding="utf-8")) for x in (
        a.griis_source,a.rail_shortest_source,a.road_shortest_source,a.circuit_source
    )]
    result=run(*z)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,"n_blocks":len(result["regions"]),
                      "no_independent_species_network_effect_claim":True},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
