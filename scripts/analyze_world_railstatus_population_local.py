#!/usr/bin/env python3
"""Sensitivity of the existing whole-landscape plant rail-status signal to
residential population and within-1000km geographic contrasts. No new raw scan.
Not a causal active-vs-abandoned comparison: status and tagging are not random.
"""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_world_rail_status_richness_railstatus import (
    CATS, LOG, ROAD, BASE, LAND_THRESHOLDS, GRID_WIDTHS, CONTINENTS, ho,
    inject_status, read, make_lookup, enrich_cells, read_source_builtup,
    ORIGINAL_NE_SHA, WORLDCLIM_SHA, EXPECTED_SHA, MD5, sha_file,
    source_continents, worldclim_tif, climate_scale, source_climate_lookup,
    original_cell_rows, verify_original_sources)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import inject_built
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
POP_SHA="846c03076da12d9d998e6b9946a0774ca134ea76b0731d11483058a7f9185bbf"
POP="log_population_density"
ACTIVE=LOG["active_normal"];INACTIVE=LOG["inactive_normal"]
FULL=BASE+ROAD+tuple(LOG[c] for c in CATS)
GRID="source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"

def join_population(rows,pop,roads):
    matched=[];qa=Counter()
    for r in rows:
        z=pop.get(r["key"])
        if z is None:qa["no_population_source"]+=1;continue
        land=float(roads[r["key"]]["source_GRIP4_total_land_area_m2"])
        n=float(z["source_GHS_POP2020_people"]);a=float(z["source_GHS_POP2020_approx_valid_pixel_area_m2"])
        if not all(map(math.isfinite,(land,n,a))) or land<=0 or min(n,a)<0:
            raise ValueError("Invalid population / area source")
        if a<.9*land:qa["population_area_less_than_90pct_GRIP_land"]+=1;continue
        matched.append(dict(r,log_population_density=math.log1p(n/(land/1e6))))
    return matched,dict(qa)

def groupkey(r,local):
    return (r["continent"],r["block"]) if local else (r["continent"],)

def local_support(rows):
    groups=defaultdict(list)
    for r in rows:groups[groupkey(r,True)].append(r)
    return [r for rr in groups.values() if len(rr)>=2 for r in rr]

def fit(rows,fields,local,target="log_species"):
    groups=defaultdict(list);blocks=defaultdict(list)
    for i,r in enumerate(rows):
        groups[groupkey(r,local)].append(i);blocks[(r["continent"],r["block"])].append(i)
    out={"n_cells":len(rows),"n_FE_strata":len(groups),"n_blocks":len(blocks),
         "n_continents":len({r["continent"] for r in rows}),
         "n_blocks_varying_active":sum(len({rows[i][ACTIVE] for i in v})>1 for v in blocks.values()),
         "n_blocks_varying_inactive":sum(len({rows[i][INACTIVE] for i in v})>1 for v in blocks.values())}
    if len(rows)-len(groups)-len(fields)<30 or len(blocks)<20:
        return dict(out,status="insufficient support")
    w=block_weights(rows)
    X=np.array([[r[k] for k in fields]+[r[target]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite source")
    for ii in groups.values():
        X[ii]-=np.average(X[ii],axis=0,weights=w[ii])
    X*=np.sqrt(w)[:,None];A=X[:,:-1];y=X[:,-1]
    scale=np.sqrt((A*A).sum(axis=0));keep=scale>1e-10
    names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE}<=set(names):return dict(out,status="no residual exposure variation")
    Z=A[:,keep]/scale[keep]
    singular=np.linalg.svd(Z,compute_uv=False)
    if singular[-1]<1e-8*singular[0]:
        return dict(out,status="rank deficient",condition_number=float(singular[0]/max(singular[-1],1e-20)))
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@y)
    residual=y-Z@b
    # Standardized design is for conditioning only; return unstandardized slopes.
    beta=b/scale[keep]
    score=np.array([Z[ii].T@residual[ii] for ii in blocks.values()])
    U=score@inv.T/scale[keep]
    cov=U.T@U*len(blocks)/(len(blocks)-1)
    out.update(status="estimated",condition_number=float(singular[0]/singular[-1]),dropped_constant_fields=[f for f in fields if f not in names])
    def lin(co):
        c=np.array([co.get(k,0.) for k in names],float)
        v=float(c@beta);se=math.sqrt(max(0.,float(c@cov@c)))
        return {"estimate":v,"block_cluster_SE":se,"interval95_approx":[v-1.96*se,v+1.96*se]}
    out["effects"]={c:lin({LOG[c]:1.}) if LOG[c] in names else None for c in CATS}
    out["active_minus_inactive"]=lin({ACTIVE:1.,INACTIVE:-1.})
    return out

def compare_cv(rows):
    out={}
    for mode,base in (("without_population",FULL),("with_population",FULL+(POP,))):
        out[mode]={}
        for c in CATS:
            x=LOG[c];reduced=tuple(f for f in base if f!=x)
            r=ho(rows,reduced,reduced+(x,),c)
            d=r["mean_dMSE"];m=r["base_mean_MSE"]
            r["relative_gain"]=d/m if d is not None and m>0 else None
            out[mode][c]=r
    return out

def analyze(rows):
    local=local_support(rows)
    both=local_support([r for r in rows if r["active_normal_km"]>0 and r["inactive_normal_km"]>0])
    high=local_support([r for r in rows if r["observations"]>=200 and r["observers"]>=20])
    return {"n_pop_supported":len(rows),"n_both_status_positive":sum(r["active_normal_km"]>0 and r["inactive_normal_km"]>0 for r in rows),
        "same_support_continent_holdout":compare_cv(rows),
        "continent_FE_before_population":fit(rows,FULL,False),
        "continent_FE_after_population":fit(rows,FULL+(POP,),False),
        "local_FE_before_population":fit(local,FULL,True),
        "local_FE_after_population":fit(local,FULL+(POP,),True),
        "local_FE_high_effort_after_population":fit(high,FULL+(POP,),True),
        "local_FE_both_positive_after_population":fit(both,FULL+(POP,),True),
        "local_FE_leave_one_continent_out":{c:fit([r for r in local if r["continent"]!=c],FULL+(POP,),True)
                                              for c in CONTINENTS if any(r["continent"]==c for r in local)}}

def main():
    p=argparse.ArgumentParser()
    for n in ("plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    plant,rail,status,tree,road,pop=map(read,(a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":
        raise ValueError("Frozen road/tree source mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Frozen climate/map mismatch")
    if pop["version"]!="v391-original-JRC-GHS-POP-2020-1km-global-counts" or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:
        raise ValueError("Frozen population source mismatch")
    if not status["version"].startswith("railstatus-v1"):raise ValueError("Rail status source mismatch")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=climate_scale(raster)
        for width in GRID_WIDTHS:
            original=plant["source_grids"][str(width)][GRID]
            keys={(x["source_equal_area_grid_ix"],x["source_equal_area_grid_iy"]) for x in original}
            clim,_=source_climate_lookup(land,raster,scale,keys,width)
            rows,_,_=original_cell_rows(original,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)],clim)
            rows,_=inject_built(rows,built[str(width)]);rows,worst=inject_status(rows,status,width)
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
            pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(width)])
            for threshold in LAND_THRESHOLDS:
                base,qa=enrich_cells(rows,rd,tr,width,threshold);joined,pqa=join_population(base,pp,rd)
                z={"width_m":width,"min_land_fraction":threshold,"n_before_population":len(base),
                    "population_exclusions":pqa,"land_exclusions":qa,"max_active_km_replay_error":worst,**analyze(joined)}
                panels.append(z)
                print(json.dumps({"population_local_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-population-local-v1","source_replay":checks,"population_SHA256":POP_SHA,"panels":panels,
        "primary":"50% land, 100/200 km, joint status slopes within same 1000km geographic blocks plus GHS population; both-positive/high-effort and 80% support sensitivity",
        "limits":["2020 modeled residential population, not travel or actual operations","local block FE does not establish causal exchangeability",
        "both-positive cells restrict support and may induce selection","incomplete abandoned/status tags","conditional on observed effort and sampled cells",
        "whole-landscape recorded plant richness, not trackside flora; intervals approximate and blocks may be spatially dependent"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
