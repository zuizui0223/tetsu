#!/usr/bin/env python3
"""Cross-taxon x rail-status double contrast on common support, adjusting jointly
for population and local geography. Same data as the two preceding analyses.
"""
import argparse,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from analyze_world_railstatus_population_local import (
    POP_SHA,POP,ACTIVE,INACTIVE,FULL,join_population,local_support,
    CATS,LOG,GRID_WIDTHS,LAND_THRESHOLDS,read,make_lookup,enrich_cells,inject_status,
    read_source_builtup,ORIGINAL_NE_SHA,WORLDCLIM_SHA,EXPECTED_SHA,MD5,sha_file,
    source_continents,worldclim_tif,climate_scale,source_climate_lookup,
    original_cell_rows,verify_original_sources,inject_built,block_weights,ho)
from analyze_world_railstatus_cross_taxa import GROUPS,key,convert,verify_replay,restrict

def flexible(r):
    o,b=r["log_observations"],r["log_observers"]
    return dict(r,effort_obs_sq=o*o,effort_obs_cu=o**3,
        effort_observers_sq=b*b,effort_observers_cu=b**3,effort_interaction=o*b)

EXTRA=("effort_obs_sq","effort_obs_cu","effort_observers_sq","effort_observers_cu","effort_interaction")
BASE=FULL+(POP,)
TARGETS={"active":{ACTIVE:1.},"inactive":{INACTIVE:1.},"active_minus_inactive":{ACTIVE:1.,INACTIVE:-1.}}

def fit(rows,fields,local):
    blocks=defaultdict(list);strata=defaultdict(list)
    for i,r in enumerate(rows):
        b=(r["continent"],r["block"]);blocks[b].append(i)
        strata[b if local else (r["continent"],)].append(i)
    qa={"n_cells":len(rows),"n_blocks":len(blocks),"n_strata":len(strata),
        "continents":dict(Counter(r["continent"] for r in rows))}
    if len(blocks)<20 or len(rows)-len(strata)-len(fields)<30:return dict(qa,status="insufficient support"),None
    w=block_weights(rows)
    X=np.array([[r[f] for f in fields]+[r["log_species"]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite source")
    for ix in strata.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None];A=X[:,:-1];y=X[:,-1]
    norms=np.sqrt((A*A).sum(axis=0));keep=norms>1e-10
    names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE}<=set(names):return dict(qa,status="no exposure variation"),None
    Z=A[:,keep]/norms[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return dict(qa,status="rank deficient"),None
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@y);e=y-Z@b
    coef=b/norms[keep]
    scores=np.array([Z[ix].T@e[ix] for ix in blocks.values()])@inv.T/norms[keep]
    estimates={};influence={}
    for t,terms in TARGETS.items():
        c=np.array([terms.get(f,0.) for f in names])
        estimates[t]=float(c@coef)
        influence[t]={k:float(v) for k,v in zip(blocks,scores@c)}
    return dict(qa,status="estimated",condition_number=float(s[0]/s[-1]),dropped_constant_fields=[f for f in fields if f not in names],effects={t:interval(estimates[t],influence[t]) for t in TARGETS}), (estimates,influence)

def interval(beta,score):
    n=len(score)
    if n<2:return None
    se=math.sqrt(n/(n-1)*sum(v*v for v in score.values()))
    return {"estimate":beta,"block_cluster_SE":se,"interval95_approx":[beta-1.96*se,beta+1.96*se],"n_blocks":n}

def paired(a,b):
    if a is None or b is None:return None
    ea,ia=a;eb,ib=b
    out={}
    for t in TARGETS:
        if set(ia[t])!=set(ib[t]):raise ValueError("Paired blocks differ")
        out[t]=interval(ea[t]-eb[t],{k:ia[t][k]-ib[t][k] for k in ia[t]})
    return out

def cv(rows,fields):
    base=tuple(f for f in fields if f!=ACTIVE)
    z=ho(rows,base,base+(ACTIVE,),"conditional active")
    d,m=z["mean_dMSE"],z["base_mean_MSE"]
    return dict(z,relative_gain=d/m if d is not None and m>0 else None)

def panel(rows_by_group,other):
    keys={r["key"] for r in rows_by_group["plants"]}&{r["key"] for r in rows_by_group[other]}
    rows={g:restrict(rows_by_group[g],keys) for g in ("plants",other)}
    # Local support is determined once from common geographic cells.
    lk={r["key"] for r in local_support(rows["plants"])}
    local={g:restrict(rows[g],lk) for g in rows}
    results={}
    for label,fields in (("linear_effort",BASE),("flexible_effort",BASE+EXTRA)):
        fits={};scores={}
        for g in rows:
            fits[g]={};scores[g]={}
            for name,rr,is_local in (("continent_FE",rows[g],False),("local_FE",local[g],True)):
                fit_result,score=fit(rr,fields,is_local);fits[g][name]=fit_result;scores[g][name]=score
            fits[g]["continent_holdout"]=cv(rows[g],fields) if rows[g] else None
        results[label]={"fits":fits,"plant_minus_other":{name:paired(scores["plants"][name],scores[other][name]) for name in ("continent_FE","local_FE")}}
    return {"other":other,"n_common_cells":len(keys),"n_local_cells":len(lk),"results":results}

def main():
    p=argparse.ArgumentParser()
    for n in ("taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    source,plant,rail,status,tree,road,pop=map(read,(a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    verify_replay(source,plant);checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Frozen climate/map mismatch")
    if pop["version"]!="v391-original-JRC-GHS-POP-2020-1km-global-counts" or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen population mismatch")
    if not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong rail status")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=climate_scale(raster)
        for width in GRID_WIDTHS:
            grid=source["grids"][str(width)]
            keys={key(z) for g in GROUPS for z in grid[g]}
            clim,_=source_climate_lookup(land,raster,scale,keys,width)
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
            pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(width)])
            geo={}
            for g in GROUPS:
                rows,_,_=original_cell_rows([convert(z) for z in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)],clim)
                rows,_=inject_built(rows,built[str(width)]);rows,_=inject_status(rows,status,width);geo[g]=rows
            for threshold in LAND_THRESHOLDS:
                rows={};qa={}
                for g in GROUPS:
                    m,lqa=enrich_cells(geo[g],rd,tr,width,threshold);m,pqa=join_population(m,pp,rd)
                    rows[g]=[flexible(r) for r in m];qa[g]={"land":lqa,"population":pqa,"n":len(m)}
                z={"width_m":width,"min_land_fraction":threshold,"support":qa,
                    "comparisons":[panel(rows,g) for g in GROUPS[1:]]}
                panels.append(z);print(json.dumps({"joint_cross_taxa_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-cross-taxa-population-local-v1","source_replay":checks,"source_qa":source["qa"],
        "primary":"plants versus birds, 50% land, pairwise matched local FE with linear effort; active minus inactive slope difference across taxa; 100/200 km both reported; insects/fungi secondary, flexible effort and 80% land sensitivity",
        "panels":panels,"limits":["same 2025 data, posthoc combined sensitivity not independent validation","recorded richness, not true richness or causal transport","approximate block-cluster intervals; shared blocks paired but spatial dependence between blocks remains","population is residential not mobility","fungal predictive folds may be insufficient","effort adjustment does not equal observer matching"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
