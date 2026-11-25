#!/usr/bin/env python3
"""Actual crossing vs noncrossing neighboring plant landscapes with rails in
both endpoints. Equal distances, exposure controlled, endpoint FE sensitivity.
"""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr
from analyze_world_railstatus_population_local import (
    POP,POP_SHA,join_population,inject_status,CATS,LOG,BASE,ROAD,GRID_WIDTHS,
    LAND_THRESHOLDS,read,make_lookup,enrich_cells,read_source_builtup,
    ORIGINAL_NE_SHA,WORLDCLIM_SHA,EXPECTED_SHA,MD5,sha_file,source_continents,
    worldclim_tif,climate_scale,source_climate_lookup,original_cell_rows,
    verify_original_sources,inject_built,block_weights)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import continent_holdout
GRID="source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"
PAIRS="source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"
ACTIVE="cross_active_normal";INACTIVE="cross_inactive_normal"
CROSS=tuple("cross_"+c for c in CATS)
CELL_FIELDS=tuple(dict.fromkeys(BASE+ROAD+(POP,"log_species")+tuple(LOG[c] for c in CATS)))
# Endpoint FE absorbs all additive mean terms, but not dissimilarity terms.
DIFF=tuple("diff_"+k for k in CELL_FIELDS)
MEAN=tuple("mean_"+k for k in CELL_FIELDS)
PRED=MEAN+DIFF+("horizontal",)+CROSS

def build_pairs(pairs,rows,edges):
    by={r["key"]:r for r in rows};out=[]
    for p in pairs:
        a=tuple(p["source_cell_A_original_grid_ix_iy"]);b=tuple(p["source_cell_B_original_grid_ix_iy"])
        if b[0]-a[0]+b[1]-a[1]!=1 or sum(abs(x-y) for x,y in zip(a,b))!=1:raise ValueError("Nonadjacent pair")
        x=by.get(a);y=by.get(b)
        if x is None or y is None or x["continent"]!=y["continent"] or x["block"]!=y["block"]:continue
        if min(x["active_normal_km"],y["active_normal_km"])<=0:continue
        if x["species"]!=p["source_grid_A_original_plant_species_n"] or y["species"]!=p["source_grid_B_original_plant_species_n"]:raise ValueError("Original species count mismatch")
        shared=p["source_shared_plants_species_n"];union=p["source_union_plants_species_n"]
        if union!=x["species"]+y["species"]-shared or union<=0:raise ValueError("Jaccard components invalid")
        j=float(p["source_biological_Jaccard_similarity_between_ENTIRE_landscape_cells"])
        if not 0<=j<=1 or abs(j-shared/union)>1e-7:raise ValueError("Jaccard mismatch")
        r={"a":a,"b":b,"continent":x["continent"],"block":x["block"],"jaccard":j,
           "simpson":shared/min(x["species"],y["species"]),"horizontal":float(b[0]!=a[0]),
           "both_inactive":min(x["inactive_normal_km"],y["inactive_normal_km"])>0}
        for k in CELL_FIELDS:r["mean_"+k]=(x[k]+y[k])/2;r["diff_"+k]=abs(x[k]-y[k])
        for c in CATS:
            n=edges[c].get((a,b),0);r["cross_"+c]=float(n>0);r["log_cross_"+c]=math.log1p(n)
            if n>0 and min(x[c+"_km"],y[c+"_km"])<=0:raise ValueError("Geometry crossing inconsistent with prior rail km")
        out.append(r)
    return out

def fit(rows,target,endpoint,fields=PRED):
    blocks=defaultdict(list)
    for i,r in enumerate(rows):blocks[(r["continent"],r["block"])].append(i)
    qa={"n_pairs":len(rows),"n_blocks":len(blocks),"n_endpoints":len({r[k] for r in rows for k in ("a","b")}),
        "n_cross_active":sum(r[ACTIVE]>0 for r in rows),"n_cross_inactive":sum(r[INACTIVE]>0 for r in rows)}
    if len(blocks)<20 or len(rows)<100:return dict(qa,status="insufficient support")
    w=block_weights(rows);X=np.array([[r[k] for k in fields]+[r[target]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite pair")
    if endpoint:
        nodes={k:i for i,k in enumerate(sorted({r[q] for r in rows for q in ("a","b")}))}
        ii=np.repeat(np.arange(len(rows)),2);jj=np.array([nodes[r[k]] for r in rows for k in ("a","b")])
        D=csr_matrix((np.repeat(np.sqrt(w),2),(ii,jj)),shape=(len(rows),len(nodes)))
        X*=np.sqrt(w)[:,None]
        for j in range(X.shape[1]):
            sol=lsmr(D,X[:,j],atol=1e-12,btol=1e-12,maxiter=4000)
            if sol[1] not in (0,1,2,4,5):raise ValueError("Endpoint FE projection failed")
            X[:,j]-=D@sol[0]
        if np.max(np.abs(D.T@X))>1e-7:raise ValueError("Endpoint FE residual check")
    else:
        for ix in blocks.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
        X*=np.sqrt(w)[:,None]
    A=X[:,:-1];y=X[:,-1];norm=np.sqrt((A*A).sum(axis=0))
    # Relative to pre-projection raw feature norm; avoid rescaling numerical FE noise.
    orig=np.sqrt(np.sum(w[:,None]*np.array([[r[k] for k in fields] for r in rows],float)**2,axis=0))
    keep=norm>np.maximum(1e-9,orig*1e-8);names=[f for f,k in zip(fields,keep) if k]
    if ACTIVE not in names:return dict(qa,status="no residual active edge variation")
    Z=A[:,keep]/norm[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return dict(qa,status="rank deficient")
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@y);e=y-Z@b;coef=b/norm[keep]
    scores=np.array([Z[ix].T@e[ix] for ix in blocks.values()])@inv.T/norm[keep]
    idx=names.index(ACTIVE)
    residual_blocks=sum(np.linalg.norm(Z[ix,idx])>1e-8 for ix in blocks.values())
    qa["n_blocks_residual_active"]=residual_blocks
    if residual_blocks<20:return dict(qa,status="insufficient residual crossing support")
    def effect(terms):
        if not set(terms)<=set(names):return None
        c=np.array([terms.get(k,0.) for k in names]);v=float(c@coef)
        se=math.sqrt(len(blocks)/(len(blocks)-1)*float(np.sum((scores@c)**2)))
        return {"estimate":v,"interval95_approx":[v-1.96*se,v+1.96*se]}
    return dict(qa,status="estimated",condition_number=float(s[0]/s[-1]),
                effects={c:effect({c:1.}) for c in CROSS},active_minus_inactive=effect({ACTIVE:1.,INACTIVE:-1.}))

def evaluate(rows):
    z={}
    for target in ("jaccard","simpson"):
        z[target]={"block_FE":fit(rows,target,False),"endpoint_FE":fit(rows,target,True,DIFF+("horizontal",)+CROSS)}
        if rows:
            base=tuple(k for k in PRED if k!=ACTIVE)
            cv=continent_holdout(rows,base,base+(ACTIVE,),target,True,50,300,target)
            z[target]["continent_holdout"]=cv
    return z

def main():
    p=argparse.ArgumentParser()
    for n in ("edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    edge,plant,rail,status,tree,road,pop=map(read,(a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail)
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"]:raise ValueError("Source edge replay invalid")
    if any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Status lengths not replayed")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=climate_scale(raster)
        for width in GRID_WIDTHS:
            g=plant["source_grids"][str(width)];original=g[GRID]
            keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for r in original}
            clim,_=source_climate_lookup(land,raster,scale,keys,width)
            rows,_,_=original_cell_rows(original,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)],clim)
            rows,_=inject_built(rows,built[str(width)]);rows,_=inject_status(rows,status,width)
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
            pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(width)])
            el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(width)][c]} for c in CATS}
            for threshold in LAND_THRESHOLDS:
                m,_=enrich_cells(rows,rd,tr,width,threshold);m,_=join_population(m,pp,rd)
                pairs=build_pairs(g[PAIRS],m,el)
                out={"width_m":width,"min_land_fraction":threshold,"main":evaluate(pairs),
                     "both_status_at_both_endpoints":evaluate([r for r in pairs if r["both_inactive"]])}
                panels.append(out);print(json.dumps({"direct_boundary_panel":out},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-direct-boundary-beta-v1","replay":checks,"panels":panels,
            "primary":"100km, land>=50%, both endpoints active positive, same 1000km block; actual crossing binary vs no crossing, Jaccard; endpoint FE sensitivity; 200km, land>=80%, Simpson similarity secondary",
            "limits":["crossing is physical geometry, not service/transport","whole-landscape recorded similarity, not railway flora","same-block restriction avoids cross-block shared endpoints; blocks not guaranteed independent","endpoint FE cannot absorb dyad-specific confounding","crossing absent from map not certified absent; endpoint-only crossing may be missed","secondary endpoints and scales not independent replications"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
