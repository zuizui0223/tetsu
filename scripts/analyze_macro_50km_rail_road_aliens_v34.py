#!/usr/bin/env python3
"""v3.4 REAL observational 50km macro plant assemblages vs RAIL/ROAD routes.

Join three independently successful, preselected source runs:
  GRIIS-introduced GBIF speciesKey × 50km cell × period (v3.1)
  Overture physical RAIL shortest-path source (v3.2)
  Overture physical ROAD shortest-path source (v3.3)

Exploratory partial rank association, no p-values or causal effect:
do rail-competitive network-distance cell pairs have more similar
observed introduced plant taxon compositions, after rank adjustment
for straight-line distance, road route distance and GBIF indexed
recording opportunities? Two independent blocks only (US/FR).
Cells/pairs within one block are dependent; never treat N pairs as
independent macro replications. All two a priori anchor definitions and
two response metrics are reported, not selected by positive effect.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json
from pathlib import Path
import math
import numpy as np

from audit_gbif_equalarea_macroblocks_v29 import PERIODS,BLOCKS
from analyze_macro_rail_assemblages_v22 import incidence_probabilities,expected_rarefied_overlap

MACRO_BLOCKS=("us_midwest","fr_burgundy")
STRATEGIES=("nearest_in_cell","nearest_largest_component_in_cell")
RECORD_RAREFACTION=20
MIN_CELL_TOTAL=20
MIN_CELL_TAXA=5
MIN_VALID_CELLS=8
MIN_VALID_PAIRS=20
MAX_COUNT_DRIFT_RATE=0.25


def _rank(x):
    """Average-tie ranks, with no SciPy requirement."""
    x=np.asarray(x,dtype=float)
    order=np.argsort(x,kind="stable")
    r=np.empty(len(x),dtype=float)
    k=0
    while k<len(x):
        j=k+1
        while j<len(x) and x[order[j]]==x[order[k]]:
            j+=1
        r[order[k:j]]=(k+j+1)/2.
        k=j
    return r


def partial_rank_r(x,y,covariates):
    """Residualized Spearman-like partial correlation. NO pairwise p-value."""
    xx=np.asarray(x,dtype=float); yy=np.asarray(y,dtype=float)
    z=np.asarray(covariates,dtype=float)
    if xx.ndim!=1 or yy.ndim!=1 or z.ndim!=2 or len(xx)!=len(yy) or len(xx)!=len(z):
        raise ValueError("Invalid source comparison input dimensions")
    if len(xx)<MIN_VALID_PAIRS or not np.all(np.isfinite(xx)) or not np.all(np.isfinite(yy)) or not np.all(np.isfinite(z)):
        return None
    if len(np.unique(xx))<5 or len(np.unique(yy))<5:
        return None
    A=np.column_stack([np.ones(len(xx)),*[_rank(z[:,k]) for k in range(z.shape[1])]])
    if np.linalg.matrix_rank(A)<A.shape[1]:
        return None
    rx=_rank(xx);ry=_rank(yy)
    rx=rx-A@np.linalg.lstsq(A,rx,rcond=None)[0]
    ry=ry-A@np.linalg.lstsq(A,ry,rcond=None)[0]
    d=float(np.linalg.norm(rx)*np.linalg.norm(ry))
    if d<1e-9:return None
    return round(float(np.dot(rx,ry)/d),6)


def _eligible_both_period_cells(block,griis):
    city=griis["regions"][block]
    cells={}
    for cid,row in city["cells"].items():
        if row.get("source_pre_outcome_both_period_recording_eligible") is not True:
            continue
        per={}
        okay=True
        for a,b in PERIODS:
            key=f"{a}_{b}"
            c=row["periods"].get(key)
            if not isinstance(c,dict) or not c.get("full_facet_source_pagination_completed") or not c.get("alien_observed_minimum_5_species_and_20_records"):
                okay=False;break
            counts=c.get("observed_alien_specieskey_counts")
            if not isinstance(counts,dict) or len(counts)<MIN_CELL_TAXA or any(
                not isinstance(k,str) or isinstance(v,bool) or not isinstance(v,int) or v<=0 for k,v in counts.items()
            ):
                raise ValueError("Invalid country introduced species key source")
            n=sum(counts.values())
            if n<MIN_CELL_TOTAL or n!=c.get("n_indexed_GRIIS_introduced_plant_records"):
                raise ValueError("Alien indexed record denominator differs from species facets")
            current=c.get("indexed_all_plantae_source_records")
            prior=c.get("prior_2026_GBIF_indexed_source_record_count")
            if not isinstance(current,int) or not isinstance(prior,int) or prior<=0 or abs(current-prior)/prior>MAX_COUNT_DRIFT_RATE:
                okay=False;break
            per[key]={"counts":counts,
                      "rarified_probs":incidence_probabilities(counts,k=RECORD_RAREFACTION),
                      "all_plant_indexed_records":current}
        if okay:cells[cid]=per
    return cells


def _routes(artifact,block,strategy,key):
    city=artifact["regions"][block]
    data=city["routes"][strategy]
    rows=data[key]
    if not isinstance(rows,list) or len(rows)!=120:
        raise ValueError("Expected complete 16 choose 2 route source")
    pairs={}
    for r in rows:
        a,b=r.get("a"),r.get("b")
        if not isinstance(a,str) or not isinstance(b,str) or a==b:
            raise ValueError("Bad geographic source cell pair")
        pair=tuple(sorted((a,b)))
        if pair in pairs:raise ValueError("Duplicate physical route pair")
        pairs[pair]=r
    return pairs


def _quartiles(values):
    if not values:return None
    return [round(float(x),5) for x in np.quantile(values,[.25,.5,.75])]


def analyze_one_block(block,griis,rail,road):
    if block not in BLOCKS or griis["regions"][block]["country"]!=BLOCKS[block][0]:
        raise ValueError("Macro source block country mismatch")
    cells=_eligible_both_period_cells(block,griis)
    result={
        "country":BLOCKS[block][0],
        "n_source_qualified_two_period_plant_cells":len(cells),
        "periods":{}
    }
    for period in (f"{p[0]}_{p[1]}" for p in PERIODS):
        strategies={}
        for strat in STRATEGIES:
            ar=_routes(rail,block,strat,"pair_routes")
            rd=_routes(road,block,strat,"pair_paths")
            if set(ar)!=set(rd):raise ValueError("Rail and road source pair universe differs")
            recs=[]
            for a,b in combinations(sorted(cells),2):
                k=(a,b)
                x=ar[k];y=rd[k]
                railratio=x.get("rail_to_geodesic_ratio")
                roadratio=y.get("road_to_geodesic_ratio")
                if not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) and .95<=v<1e3
                           for v in (railratio,roadratio)):
                    continue
                railpath=x.get("cell_centroid_snap_augmented_rail_path_km")
                roadpath=y.get("centroid_snap_augmented_road_path_km")
                if not all(isinstance(v,(int,float)) and v>0 and math.isfinite(v) for v in (railpath,roadpath)):
                    continue
                pa=cells[a][period];pb=cells[b][period]
                keys=set(pa["counts"])|set(pb["counts"])
                if not keys:continue
                J=len(set(pa["counts"])&set(pb["counts"]))/len(keys)
                rare=expected_rarefied_overlap(pa["rarified_probs"],pb["rarified_probs"])
                if rare is None:continue
                n1,n2=pa["all_plant_indexed_records"],pb["all_plant_indexed_records"]
                if n1<=0 or n2<=0:raise ValueError("Missing source record denominator")
                recs.append({
                    "a":a,"b":b,
                    "log_rail_vs_road_route":math.log(railpath/roadpath),
                    "geodesic_km":x["source_cell_centroid_geodesic_km"],
                    "road_network_km":roadpath,
                    "log_source_plant_count_ratio":abs(math.log(n1/n2)),
                    "log_source_plant_count_geomean":(math.log(n1)+math.log(n2))/2.,
                    "source_observed_alien_jaccard":J,
                    "source_expected_rarefied20_intersection_union_ratio":rare,
                })
            summary={
                "n_plant_qualified_cells_with_complete_macro_routes":len(set(
                    p["a"] for p in recs)|set(p["b"] for p in recs)),
                "n_plant_pair_routes_with_both_network_lengths":len(recs),
                "network_rail_over_road_path_ratio_quantiles":_quartiles([
                    math.exp(r["log_rail_vs_road_route"]) for r in recs]),
                "observed_source_alien_jaccard_quantiles":_quartiles([
                    r["source_observed_alien_jaccard"] for r in recs]),
                "n_rail_more_direct_than_road_pairs":sum(
                    r["log_rail_vs_road_route"]<0 for r in recs),
                "responses":{},
            }
            for key in ("source_expected_rarefied20_intersection_union_ratio",
                        "source_observed_alien_jaccard"):
                X=[r["log_rail_vs_road_route"] for r in recs]
                Y=[r[key] for r in recs]
                Z=[[math.log(r["geodesic_km"]),math.log(r["road_network_km"]),
                    r["log_source_plant_count_ratio"],r["log_source_plant_count_geomean"]]
                   for r in recs]
                rho=partial_rank_r(X,Y,Z)
                summary["responses"][key]={
                    "descriptive_partial_rank_correlation":rho,
                    "n_paired_observations_NOT_independent":len(recs),
                    "n_unique_cells_source_qualified":len(cells),
                    "sign_relevant_to_rail_advantage":"negative correlation means relatively shorter rail journeys associated with higher observed taxon similarity",
                    "no_p_value_due_to_shared_cells_and_two_macro_blocks":True,
                }
            strategies[strat]=summary
        result["periods"][period]=strategies
    return result


def run(griis,rail,road):
    if not griis.get("complete") or not rail.get("complete") or not road.get("complete"):
        raise ValueError("Incomplete required independent sources")
    if set(road.get("regions",{}))!=set(MACRO_BLOCKS):
        raise ValueError("Road source region frame differs from fixed US FR pilot")
    if not set(MACRO_BLOCKS).issubset(rail.get("regions",{})) or not set(MACRO_BLOCKS).issubset(griis.get("regions",{})):
        raise ValueError("Rail and GRIIS source regions missing")
    result={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "status":"FIRST source-aligned 50km EQUAL-AREA two-block observational network and plant community comparison",
        "n_independent_geographic_blocks":2,
        "independent_geographic_replication_sufficient":False,
        "regions":{},
        "NO_causal_propagule_or_evolution_inference":True,
    }
    for block in MACRO_BLOCKS:
        v=analyze_one_block(block,griis,rail,road)
        result["regions"][block]=v
        for period in ("2010_2017","2018_2025"):
            for strat in STRATEGIES:
                row=v["periods"][period][strat]
                print(json.dumps({
                    "block":block,"period":period,"anchor_strategy":strat,
                    "n_qualified_two_period_cells":v["n_source_qualified_two_period_plant_cells"],
                    "n_valid_pairs":row["n_plant_pair_routes_with_both_network_lengths"],
                    "rail_vs_road_route_ratio_quartiles":row["network_rail_over_road_path_ratio_quantiles"],
                    "partial_rank_rarefied20":row["responses"]["source_expected_rarefied20_intersection_union_ratio"]["descriptive_partial_rank_correlation"],
                    "partial_rank_observed_Jaccard":row["responses"]["source_observed_alien_jaccard"]["descriptive_partial_rank_correlation"],
                },ensure_ascii=False),flush=True)
    result["limitations"]=[
        "Exactly two selected wide 50km source blocks (US Midwest and France Burgundy), so do NOT present a globally significant ecological relationship.",
        "The 2026 physical Overture networks are not historical train/freight routes for 2010–2017 or 2018–2025 GBIF observation periods.",
        "The 20-record rarefied GBIF facets and Jaccard are based on positive indexed source records, not comprehensive field assemblages; absence is not proven.",
        "Per-pair partial rank correlations control only physical Euclidean separation, road route distance and simple GBIF record-count proxies; climate/land/ecoregion/ports/observer sources absent.",
        "Cell centroid to rail/road nearest connector is a rough access proxy with possible off-cell botanical records.",
        "Road paths omit one-way/vehicle/freight permissions, traffic and graph paths beyond halo; rail paths omit active service.",
        "Plant species-key catalog is GRIIS country introduced lists, which vary in completeness; nonmembers retain unknown status.",
        "Pairwise observations within a 16-cell block are NOT independent; outputs deliberately have no pair-level p-value or confidence interval.",
        "Both predeclared anchor strategies and two abundance-normalization response definitions are reported regardless of effect direction.",
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--griis-macro-source",required=True)
    p.add_argument("--rail-path-source",required=True)
    p.add_argument("--road-path-source",required=True)
    p.add_argument("--output",default="results/audit/macro_50km_rail_road_aliens_v34.json")
    x=p.parse_args()
    r=run(*[json.loads(Path(f).read_text(encoding="utf-8")) for f in (
        x.griis_macro_source,x.rail_path_source,x.road_path_source)])
    dest=Path(x.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"completed":True,"two_block_source_comparison_only":True,
                      "n_regions":len(r["regions"])},ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
