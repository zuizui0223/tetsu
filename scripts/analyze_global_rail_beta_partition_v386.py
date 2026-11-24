#!/usr/bin/env python3
"""v386 source Jaccard dissimilarity partition, same v382/v385 samples."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_global_rail_road_tree_richness_v382 import (
    enrich_cells,make_lookup,pair_enrich,PM,LAND_THRESHOLDS)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    verify_original_sources,source_climate_lookup,original_cell_rows,continent_holdout,
    GRID_WIDTHS,CONTINENTS,MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,block_weights)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import read_source_builtup,inject_built,ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
from validate_global_rail_local_matched_placebo_v385 import EXPECTED
TARGETS=("beta_jac","beta_jtu","beta_jne")
RAIL=("log_mean_rail_km","absdiff_log_rail_km")
FIELDS=PM["both"]
EXPECTED_MEAN_RAIL_JACCARD={
 (100000,.5):-.003406061772347097,(100000,.8):-.003004747423075928,
 (200000,.5):-.0009290443529109554,(200000,.8):-.0008299615657538045}
def partition(a,s1,s2,union):
    if any(not isinstance(x,int) or isinstance(x,bool) for x in (a,s1,s2,union)):
        raise ValueError("Original counts must be exact integers")
    b,c=s1-a,s2-a
    if min(a,b,c)<0 or min(s1,s2)<=0 or union!=a+b+c or union<=0:
        raise ValueError("Inconsistent original pair species counts")
    jac=(b+c)/union
    m=min(b,c)
    # With nonempty endpoint floras, a+2m is always positive.
    jtu=2*m/(a+2*m)
    jne=jac-jtu
    if not -1e-12<=jne<=1 or not 0<=jtu<=1:
        raise ValueError("Invalid Jaccard partition")
    return jac,jtu,max(0.,jne)

def pair_partition(neighbors,cells):
    rows=pair_enrich(neighbors,cells)
    bycell={r["key"]:r for r in cells}
    source=[]
    seen=set()
    for p in neighbors:
        a=tuple(p["source_cell_A_original_grid_ix_iy"]);b=tuple(p["source_cell_B_original_grid_ix_iy"])
        if a not in bycell or b not in bycell or bycell[a]["continent"]!=bycell[b]["continent"]:continue
        if (a,b) in seen:raise ValueError("Duplicate source pair")
        seen.add((a,b));source.append((p,a,b))
    if len(rows)!=len(source):raise ValueError("Pair membership differs from exact v382")
    for row,(p,a,b) in zip(rows,source):
        vals=partition(p["source_shared_plants_species_n"],p["source_grid_A_original_plant_species_n"],
                       p["source_grid_B_original_plant_species_n"],p["source_union_plants_species_n"])
        if abs(1-row["jaccard"]-vals[0])>5.0001e-9:
            raise ValueError("Original rounded Jaccard not reproduced by source counts")
        row.update(dict(zip(TARGETS,vals)))
        row["endpoint_blocks"]=tuple(sorted({(row["continent"],bycell[k]["block"]) for k in (a,b)}))
    return rows

def direction(rows):
    if len(rows)<100:return None
    X=np.asarray([[r[k] for k in FIELDS+TARGETS] for r in rows],float)
    w=block_weights(rows)
    by=defaultdict(list)
    for i,r in enumerate(rows):by[r["continent"]].append(i)
    for ids in by.values():X[ids]-=np.average(X[ids],axis=0,weights=w[ids])
    X*=np.sqrt(w)[:,None]
    out={}
    for var in RAIL:
        i=FIELDS.index(var)
        z=X[:,i];Y=X[:,-3:];other=np.delete(X[:,:-3],i,axis=1)
        rz=z-other@np.linalg.lstsq(other,z,rcond=1e-10)[0]
        ry=Y-other@np.linalg.lstsq(other,Y,rcond=1e-10)[0]
        denom=float(rz@rz)
        if denom<=1e-10:return None
        slopes=rz@ry/denom
        if abs(slopes[0]-slopes[1]-slopes[2])>1e-10:
            raise ValueError("Component regression slopes do not add to total")
        out[var]=dict(zip(TARGETS,map(float,slopes)))
    return out

def delete_blocks(rows):
    # Delete ALL pairs touching either endpoint of the deleted spatial block;
    # weights and continent means refitted, unlike simply deleting pair centers.
    blocks=sorted({b for r in rows for b in r["endpoint_blocks"]})
    vals={v:{t:[] for t in TARGETS} for v in RAIL}
    for b in blocks:
        coef=direction([r for r in rows if b not in r["endpoint_blocks"]])
        if coef is None:continue
        for v in RAIL:
            for t in TARGETS:vals[v][t].append(coef[v][t])
    return {"n_original_endpoint_blocks":len(blocks),"ranges":{
        v:{t:{"n_valid":len(s),"min":min(s) if s else None,"max":max(s) if s else None,
               "positive":sum(x>0 for x in s),"negative":sum(x<0 for x in s)}
           for t,s in targets.items()} for v,targets in vals.items()}}

def panel(cells,neighbors,w,g):
    expected_n,expected_p,_,expected_j=EXPECTED[w,g]
    rows=pair_partition(neighbors,cells)
    if len(cells)!=expected_n or len(rows)!=expected_p:raise ValueError("v382 sample mismatch")
    prior=continent_holdout(rows,PM["base_road"],PM["both"],"jaccard",True,
        MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,"v386 v382 replay")
    if abs(prior["mean_holdout_continent_delta_MSE"]-expected_j)>1e-9:
        raise ValueError("v382 original Jaccard holdout mismatch")
    coef=direction(rows)
    if coef is None or abs(coef[RAIL[0]]["beta_jac"]+EXPECTED_MEAN_RAIL_JACCARD[w,g])>1e-8:
        raise ValueError("v385 signed total Jaccard coefficient not reproduced")
    cv={t:continent_holdout(rows,PM["base_road"],PM["both"],t,True,
        MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,t) for t in TARGETS}
    out={"width_m":w,"min_land_fraction":g,"n_cells":len(cells),"n_pairs":len(rows),
         "original_v382_Jaccard_replayed":True,
         "blockweighted_mean_component":{t:float(np.average([r[t] for r in rows],weights=block_weights(rows))) for t in TARGETS},
         "conditional_descriptive_coefficients":coef,
         "endpoint_block_deletion":delete_blocks(rows),
         "by_continent_direction":{c:direction([r for r in rows if r["continent"]==c]) for c in CONTINENTS},
         "heldout_component_models":cv}
    print(json.dumps({"v386_panel":{k:v for k,v in out.items() if k!="heldout_component_models"},
                      "component_heldout_delta":{t:v["mean_holdout_continent_delta_MSE"] for t,v in cv.items()}},allow_nan=False),flush=True)
    return out
def main():
    p=argparse.ArgumentParser()
    for k in ("plant","rail","built","tree","road","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text())
    plant,rail,tree,road=map(read,(a.plant,a.rail,a.tree,a.road))
    checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong original tree source")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or not road["version"].startswith("v381b"):raise ValueError("Wrong original road source")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong map/climate")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne)
    raster,mem,_=worldclim_tif(a.climate);out=[]
    try:
        scale,_=climate_scale(raster)
        for w in GRID_WIDTHS:
            g=plant["source_grids"][str(w)]
            cells=g["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            neighbors=g["source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"]
            keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for r in cells}
            clim,_=source_climate_lookup(land,raster,scale,keys,w)
            rows,_,_=original_cell_rows(cells,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
            rows,_=inject_built(rows,built[str(w)])
            if len(rows)!={100000:3455,200000:1442}[w]:raise ValueError("Original v376c panel mismatch")
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
            for gate in LAND_THRESHOLDS:
                m,_=enrich_cells(rows,rd,tr,w,gate);out.append(panel(m,neighbors,w,gate))
    finally:raster.close();mem.close()
    result={"version":"v386-exact-world-landscape-Jaccard-turnover-nestedness","source_checks":checks,"panels":out,
        "limits":["posthoc description of original v382/v385 samples, not independent field replication",
         "beta_jne is nestedness-resultant dissimilarity, not raw richness difference or a demonstrated process",
         "observed presence-only landscape checklists and unequal detection, not actual track vegetation",
         "1000km endpoint-block deletion range is sensitivity, not confidence interval or p value",
         "road and tree historical proxies; remaining access/population/spatial confounding"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
