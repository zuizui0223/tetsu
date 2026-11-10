#!/usr/bin/env python3
"""v4.2: source-observed GRIIS plants vs independently measured RAIL and ROAD circuits.

Exploratory posthoc two-block check, NOT a replicated biological effect:
 - rail circuit redundancy from ALL mapped RAIL physical connectors
 - road circuit redundancy from ONLY four major road classes
 - same 50km cell pair IDs / 75km halo / independent GBIF+GRIIS source
 - physical ROAD shortest path includes more classes than the ROAD circuit
 - source 2026 routes do not identify historical freight train transport

All pairwise cells are dependent within the 2 macro regions. No p-values.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json,math
from pathlib import Path

from analyze_macro_50km_rail_road_aliens_v34 import (
    MACRO_BLOCKS,STRATEGIES,partial_rank_r,
    _eligible_both_period_cells,_routes,_quartiles
)
from analyze_macro_alien_circuit_increment_v38 import circuit_pairs
from analyze_macro_rail_assemblages_v22 import expected_rarefied_overlap

PERIODS=("2010_2017","2018_2025")
MIN_COUNT=20


def source_pairs(source,block,strategy):
    data=source["regions"][block]
    if data.get("road_backbone_not_full_network") is not True:
        raise ValueError("ROAD circuit source class/fidelity not declared")
    pairs=data["strategies"][strategy].get("pair_source_diagnostics")
    if not isinstance(pairs,list) or len(pairs)!=120:
        raise ValueError("Incomplete 16-cell ROAD circuit physical source")
    out={}
    for row in pairs:
        k=tuple(sorted((row["a"],row["b"])))
        if k in out or k[0]==k[1]:
            raise ValueError("Duplicated ROAD source cell pair")
        out[k]=row
    return out


def correlate_block(block,flora,rail_short,road_short,rail_circuit,road_circuits):
    elig=_eligible_both_period_cells(block,flora)
    if len(elig)<8:
        raise ValueError("Insufficient independent eligible source cells")
    out={"country":flora["regions"][block]["country"],
         "n_independent_geographic_plant_cells":len(elig),
         "periods":{}}
    for strategy in STRATEGIES:
        rshort=_routes(rail_short,block,strategy,"pair_routes")
        dshort=_routes(road_short,block,strategy,"pair_paths")
        rcirc=circuit_pairs(rail_circuit,block,strategy)
        dcirc=source_pairs(road_circuits,block,strategy)
        if not(set(rshort)==set(dshort)==set(rcirc)==set(dcirc)):
            raise ValueError("Rail/road physical source cell pair frame differs")
        for period in PERIODS:
            rows=[]
            for a,b in combinations(sorted(elig),2):
                k=(a,b);rr=rshort[k];dd=dshort[k];r=rcirc[k];d=dcirc[k]
                qrail=r.get("independent_path_redundancy_1minusRoverD")
                qroad=d.get("source_backbone_Q_one_minus_Reff_over_shortest")
                roadkm=dd.get("centroid_snap_augmented_road_path_km")
                railkm=rr.get("cell_centroid_snap_augmented_rail_path_km")
                geo=rr.get("source_cell_centroid_geodesic_km")
                if not all(isinstance(x,(int,float)) and not isinstance(x,bool)
                           and math.isfinite(x) for x in (qrail,qroad,roadkm,railkm,geo)):
                    continue
                if not (.0<=qrail<=1.000001 and .0<=qroad<=1.000001
                        and roadkm>0 and railkm>0 and geo>0):
                    raise ValueError("Unphysical resistor source or path distance")
                f1=elig[a][period];f2=elig[b][period]
                sim=expected_rarefied_overlap(f1["rarified_probs"],f2["rarified_probs"])
                if sim is None:continue
                n1,n2=(f1["all_plant_indexed_records"],f2["all_plant_indexed_records"])
                if n1<=0 or n2<=0:
                    raise ValueError("Invalid GBIF source record denominators")
                rows.append({
                    "Qrail":qrail,"Qroad":qroad,
                    "loggeo":math.log(geo),"lograil":math.log(railkm),
                    "logroad":math.log(roadkm),
                    "count_ratio":abs(math.log(n1/n2)),
                    "count_geomean":(math.log(n1)+math.log(n2))/2,
                    "plant_similarity":sim
                })
            result={
                "n_dependent_source_pairs_with_both_circuits":len(rows),
                "source_qualified_50km_cells":len(elig),
                "Qrail_quartiles":_quartiles([x["Qrail"] for x in rows]),
                "Qroad_major_backbone_quartiles":_quartiles([x["Qroad"] for x in rows]),
                "n_pairs_road_backbone_Q_at_least_5pct":sum(x["Qroad"]>=.05 for x in rows),
                "partial_rank_association":None,
                "NO_pairwise_p_values":True,
            }
            if len(rows)>=MIN_COUNT:
                base=[[x["loggeo"],x["lograil"],x["logroad"],
                       x["count_ratio"],x["count_geomean"]]
                      for x in rows]
                yy=[x["plant_similarity"] for x in rows]
                result["partial_rank_association"]={
                    "rail_Q_given_distance_roadroute_and_source_counts":
                        partial_rank_r([x["Qrail"] for x in rows],yy,base),
                    "road_backbone_Q_given_distance_roadroute_and_source_counts":
                        partial_rank_r([x["Qroad"] for x in rows],yy,base),
                    "rail_Q_given_road_backbone_Q_as_extra_covariate":
                        partial_rank_r([x["Qrail"] for x in rows],yy,
                                       [z+[x["Qroad"]] for z,x in zip(base,rows)]),
                    "road_Q_given_rail_Q_as_extra_covariate":
                        partial_rank_r([x["Qroad"] for x in rows],yy,
                                       [z+[x["Qrail"]] for z,x in zip(base,rows)])
                }
            out["periods"].setdefault(period,{})[strategy]=result
            rho=result["partial_rank_association"] or {}
            print(json.dumps({
                "block":block,"period":period,"anchor":strategy,
                "source_pairs":len(rows),
                "rail_Q_median":result["Qrail_quartiles"][1] if result["Qrail_quartiles"] else None,
                "road_major_Q_median":result["Qroad_major_backbone_quartiles"][1] if result["Qroad_major_backbone_quartiles"] else None,
                "rho_rail_Q_given_road_Q":rho.get("rail_Q_given_road_backbone_Q_as_extra_covariate"),
                "rho_road_Q_given_rail_Q":rho.get("road_Q_given_rail_Q_as_extra_covariate")
            },ensure_ascii=False),flush=True)
    return out


def run(flora,rail_short,road_short,rail_circuit,road_us,road_fr):
    if not all(x.get("complete") for x in (
            flora,rail_short,road_short,rail_circuit,road_us,road_fr)):
        raise ValueError("At least one indispensable real network or flora source absent")
    if set(road_short.get("regions",{}))!=set(MACRO_BLOCKS):
        raise ValueError("The frozen two-block road control frame changed")
    combined={"regions":{}}
    if set(road_us.get("regions",{}))!={"us_midwest"} or set(road_fr.get("regions",{}))!={"fr_burgundy"}:
        raise ValueError("Independent US and French road-circuit input regions changed")
    combined["regions"].update(road_us["regions"])
    combined["regions"].update(road_fr["regions"])
    output={
        "utc":datetime.now(timezone.utc).isoformat(),
        "question":"Does rail circuit path redundancy have additional exploratory source-flora association after controlling MAJOR road backbone circuit redundancy?",
        "only_two_independent_geographic_blocks":True,
        "rail_source_class_vs_road_backbone_class_not_equivalent":True,
        "not_actual_train_seed_transport":True,
        "regions":{}
    }
    for block in MACRO_BLOCKS:
        if not all(block in x.get("regions",{}) for x in (flora,rail_short,road_short,rail_circuit)):
            raise ValueError("Missing source region in rail/flora data")
        output["regions"][block]=correlate_block(block,flora,rail_short,
                         road_short,rail_circuit,combined)
    output["caveats"]=[
        "Five methods combined; pilot is posthoc, not preregistered and has only US Midwest and France Burgundy independent macro blocks.",
        "Road circuit includes only primary/secondary/trunk/motorway physical network while road-shortest path includes more road classes.",
        "Electrical source resistors are arbitrary physical road and rail kilometers and do not model actual freight, trains, station transfers or seeds.",
        "Both railway and road Q are dependent on network geometry/connector representation. Comparisons are structural, not equivalent biological conductance.",
        "No actual ecological absences, original GBIF DOI event downloads, historical route maps, ecoregion, climate, urban/port or observer survey controls.",
        "Cell pair observations are overlapping, so NO pairwise p-values or confidence intervals; an apparent signed partial rank is not generalizable biology."
    ]
    return output


def main():
    ap=argparse.ArgumentParser()
    for name in ("plants","rail","road","rail-circuit","road-circuit-us","road-circuit-fr"):
        ap.add_argument("--"+name,required=True)
    ap.add_argument("--output",default="results/audit/rail_circuit_vs_major_road_circuit_plant_v42.json")
    a=ap.parse_args()
    src=[json.loads(Path(getattr(a,k)).read_text(encoding="utf-8")) for k in (
        "plants","rail","road","rail_circuit","road_circuit_us","road_circuit_fr"
    )]
    result=run(*src)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,"n_macro_blocks":len(result["regions"]),
        "source_only_not_causal":result["not_actual_train_seed_transport"]}),flush=True)

if __name__=="__main__":
    main()
