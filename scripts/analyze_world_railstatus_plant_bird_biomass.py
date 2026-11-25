#!/usr/bin/env python3
"""Paired plant and bird recorded richness on identical AGB-supported cells."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_world_railstatus_cross_taxa_population_local import (
    POP_SHA,POP,ACTIVE,INACTIVE,FULL,join_population,local_support,GRID_WIDTHS,
    LAND_THRESHOLDS,read,make_lookup,enrich_cells,inject_status,read_source_builtup,
    ORIGINAL_NE_SHA,WORLDCLIM_SHA,EXPECTED_SHA,MD5,sha_file,source_continents,
    worldclim_tif,climate_scale,source_climate_lookup,original_cell_rows,
    verify_original_sources,inject_built,block_weights,flexible,EXTRA,verify_replay,key,convert)
from aggregate_world_ESA_CCI_biomass import DOI,URL
BIO="log_woody_AGB_center50";IA="active_x_agb";II="inactive_x_agb"
BASE=FULL+(POP,"bird_log_observations","bird_log_observers")
FLEX=BASE+EXTRA+tuple("bird_"+x for x in EXTRA)
OUT=("plant_log_richness","bird_log_richness","plant_minus_bird")
TERMS={"active":{ACTIVE:1.},"inactive":{INACTIVE:1.},"active_minus_inactive":{ACTIVE:1.,INACTIVE:-1.},
       "woody_AGB":{BIO:1.},"active_x_AGB":{IA:1.},"inactive_x_AGB":{II:1.},
       "active_minus_inactive_x_AGB":{IA:1.,II:-1.},"plant_richness":{OUT[0]:1.}}

def fit(rows,fields,targets):
    g=defaultdict(list)
    for i,r in enumerate(rows):g[r["continent"],r["block"]].append(i)
    qa={"n_cells":len(rows),"n_blocks":len(g),"continents":dict(Counter(r["continent"] for r in rows))}
    if len(rows)<100 or len(g)<20 or len(rows)-len(g)-len(fields)<30:return dict(qa,status="insufficient support")
    w=block_weights(rows);X=np.array([[r[k] for k in fields+targets] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite joint design")
    for ix in g.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None];A=X[:,:len(fields)];Y=X[:,len(fields):]
    scale=np.sqrt(np.sum(A*A,axis=0));keep=scale>1e-10;names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE}<=set(names):return dict(qa,status="no status variation")
    Z=A[:,keep]/scale[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return dict(qa,status="rank deficient")
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@Y);e=Y-Z@b;coef=b/scale[keep,None]
    score=np.array([(inv@(Z[ix].T@e[ix]))/scale[keep,None] for ix in g.values()])
    effects={}
    for label,terms in TERMS.items():
        if not set(terms)<=set(names):continue
        c=np.array([terms.get(k,0.) for k in names]);v=c@coef
        if targets==OUT and abs(v[0]-v[1]-v[2])>1e-8:raise ValueError("Paired contrast identity failed")
        u=np.einsum("i,bio->bo",c,score);se=np.sqrt(len(g)/(len(g)-1)*(u*u).sum(axis=0))
        effects[label]={k:{"estimate":float(a),"interval95_approx":[float(a-1.96*b),float(a+1.96*b)]} for k,a,b in zip(targets,v,se)}
    return dict(qa,status="estimated",condition_number=float(s[0]/s[-1]),effects=effects)

def join(rows,agb20,agb24,width,roads):
    p={r["key"]:r for r in rows["plants"]};b={r["key"]:r for r in rows["birds"]}
    keys=set(p)&set(b);out=[];qa=Counter()
    for k in sorted(keys):
        x,y=p[k],b[k]
        if any(x[f]!=y[f] for f in FULL+(POP,) if f not in ("log_observations","log_observers")):raise ValueError("Geographic fields differ")
        a,c=agb20.get(k),agb24.get(k);land=roads[k]["source_GRIP4_total_land_area_m2"]
        if a is None or c is None:qa["missing_AGB"]+=1;continue
        if any(not 0<z["valid_area_m2"]<=width*width*1.000001 or not math.isfinite(z["agb_Mg_ha"]) or z["agb_Mg_ha"]<0 for z in (a,c)):raise ValueError("Invalid AGB support")
        if min(a["valid_area_m2"],c["valid_area_m2"])<.9*land:qa["AGB_area_below_90pct_GRIP_land"]+=1;continue
        z=dict(x,plant_log_richness=x["log_species"],bird_log_richness=y["log_species"],plant_minus_bird=x["log_species"]-y["log_species"],
               bird_observations=y["observations"],bird_observers=y["observers"],agb2020=a["agb_Mg_ha"],agb2024=c["agb_Mg_ha"])
        for f in ("log_observations","log_observers",*EXTRA):z["bird_"+f]=y[f]
        out.append(z)
    return out,{"n_joint_before_AGB":len(keys),"excluded":dict(qa)}

def with_epoch(rows,year):
    out=[]
    for r in rows:
        v=math.log1p(r["agb"+str(year)])-math.log1p(50.)
        out.append(dict(r,**{BIO:v,IA:r[ACTIVE]*v,II:r[INACTIVE]*v}))
    return out

def models(rows,year,flex=False):
    r=with_epoch(local_support(rows),year);fields=FLEX if flex else BASE
    return {"n":len(r),"without_biomass":fit(r,fields,OUT),"additive_biomass":fit(r,fields+(BIO,),OUT),
            "biomass_status_interaction":fit(r,fields+(BIO,IA,II),OUT),"biomass_as_outcome":fit(r,fields,(BIO,)),
            "bird_plant_richness_association":fit(r,fields+(BIO,OUT[0]),(OUT[1],))}
def panel(rows):
    return {"n_biomass_supported":len(rows),"primary_2024":models(rows,2024),"flexible_effort_2024":models(rows,2024,True),
            "same_cells_2020":models(rows,2020),
            "high_effort_both_taxa":models([r for r in rows if min(r["observations"],r["bird_observations"])>=200 and min(r["observers"],r["bird_observers"])>=20],2024),
            "both_status_positive":models([r for r in rows if r["active_normal_km"]>0 and r["inactive_normal_km"]>0],2024),
            "leave_continent_out":{c:models([r for r in rows if r["continent"]!=c],2024) for c in sorted({r["continent"] for r in rows})}}

def main():
    p=argparse.ArgumentParser()
    for n in ("agb","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();agb,source,plant,rail,status,tree,road,pop=map(read,(a.agb,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    verify_replay(source,plant);checks=verify_original_sources(plant,rail)
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or agb["DOI"]!=DOI or agb["source_url"]!=URL or set(agb["grids"])!={"2020","2024"}:raise ValueError("Wrong original AGB artifact")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    if not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong rail status")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=climate_scale(raster)
        for width in GRID_WIDTHS:
            grid=source["grids"][str(width)]
            keys={key(z) for g in ("plants","birds") for z in grid[g]}
            clim,_=source_climate_lookup(land,raster,scale,keys,width)
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
            pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(width)])
            ag20=make_lookup(agb["grids"]["2020"][str(width)]);ag24=make_lookup(agb["grids"]["2024"][str(width)])
            geo={}
            for g in ("plants","birds"):
                rows,_,_=original_cell_rows([convert(z) for z in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)],clim)
                rows,_=inject_built(rows,built[str(width)]);rows,_=inject_status(rows,status,width);geo[g]=rows
            for threshold in LAND_THRESHOLDS:
                rows={};qa={}
                for g in ("plants","birds"):
                    m,lqa=enrich_cells(geo[g],rd,tr,width,threshold);m,pqa=join_population(m,pp,rd)
                    rows[g]=[flexible(r) for r in m];qa[g]={"land":lqa,"population":pqa,"n":len(m)}
                joined,bqa=join(rows,ag20,ag24,width,rd)
                z={"width_m":width,"min_land_fraction":threshold,"support":qa,"biomass_join":bqa,**panel(joined)}
                panels.append(z);print(json.dumps({"biomass_cross_taxa_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-plant-bird-woody-biomass-association-v1","source_replay":checks,
            "biomass_source":{k:v for k,v in agb.items() if k!="grids"},"panels":panels,
            "primary":"100km land80, 2024 woody AGB, same cells with both taxa effort and local block FE; direct plant-minus-bird difference in active-minus-inactive x centered log(1+AGB) interaction",
            "limits":["woody AGB, not all plant biomass or productivity","recorded richness, not animal abundance/biomass","mapped 2026 status, not measured operation","conditional observational association, no causal mediation","joint geographic and observation selection","approximate block intervals"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
