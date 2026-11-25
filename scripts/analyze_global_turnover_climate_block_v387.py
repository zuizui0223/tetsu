#!/usr/bin/env python3
"""v387: stricter regional environmental and spatial control of v386 turnover."""
import argparse,json,math,zipfile
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from pyproj import Transformer,Geod
from rasterio.io import MemoryFile
from analyze_global_rail_beta_partition_v386 import pair_partition,direction,EXPECTED,RAIL,TARGETS,FIELDS
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,PM,LAND_THRESHOLDS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    verify_original_sources,source_climate_lookup,original_cell_rows,continent_holdout,block_weights,
    GRID_WIDTHS,CONTINENTS,MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import read_source_builtup,inject_built,ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file,read_temp
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
BIO=(4,12,15)
BIOFIELDS=("bio4_raw","log1p_bio12","bio15_raw")
EXTRA=tuple(x+k for k in BIOFIELDS for x in ("mean_","absdiff_"))+("log_geodesic_km",)
FULL=FIELDS+EXTRA
EXPECTED_TURNOVER={(100000,.5):.00579548466154781,(100000,.8):.0056617046631512544,
                   (200000,.5):.002259228697217621,(200000,.8):.0005373706045001456}
def center(key,width,tf):
    return tf.transform(X_START+(key[0]+.5)*width,Y_START+(key[1]+.5)*width)
def load_extra_climate(archive,keys,width):
    tf=Transformer.from_crs(6933,4326,always_xy=True)
    coords={k:center(k,width,tf) for k in keys}
    out={k:{} for k in keys};qa={}
    with zipfile.ZipFile(archive) as z:
        for b,f in zip(BIO,BIOFIELDS):
            members=[n for n in z.namelist() if n.endswith("bio_"+str(b)+".tif")]
            if len(members)!=1:raise ValueError("Original WorldClim BIO member not unique")
            with MemoryFile(z.read(members[0])) as mf,mf.open() as src:
                if src.crs is None or src.crs.to_epsg()!=4326 or src.count!=1:
                    raise ValueError("Original climate raster CRS mismatch")
                good=[]
                for k,(lon,lat) in coords.items():
                    v=read_temp(src,lat,lon)
                    if v is None or v<0:continue
                    out[k][f]=math.log1p(v) if b==12 else v
                    good.append(v)
                if not good:raise ValueError("No original extra climate support")
                qa[f]={"source_member":members[0],"n_valid":len(good),"raw_min":min(good),"raw_max":max(good)}
    complete={k:v for k,v in out.items() if len(v)==len(BIOFIELDS)}
    return complete,coords,qa
def add_pair_covariates(neighbors,cells,climate,coords):
    rows=pair_partition(neighbors,cells)
    bycell={r["key"]:r for r in cells}
    src=[]
    for p in neighbors:
        a=tuple(p["source_cell_A_original_grid_ix_iy"]);b=tuple(p["source_cell_B_original_grid_ix_iy"])
        if a in bycell and b in bycell and bycell[a]["continent"]==bycell[b]["continent"]:src.append((a,b))
    if len(src)!=len(rows):raise ValueError("Pair ordering/source support mismatch")
    geo=Geod(ellps="WGS84");out=[]
    for r,(a,b) in zip(rows,src):
        if a not in climate or b not in climate:continue
        q=dict(r)
        for k in BIOFIELDS:
            x,y=climate[a][k],climate[b][k]
            q["mean_"+k]=(x+y)/2
            q["absdiff_"+k]=abs(x-y)
        dist=geo.inv(*coords[a],*coords[b])[2]/1000
        if not math.isfinite(dist) or dist<=0:raise ValueError("Invalid source adjacent-cell geodesic separation")
        q["log_geodesic_km"]=math.log(dist)
        # Both pair-center and endpoint grouping retained for spatial sensitivity.
        q["both_endpoints_same_block"]=len(q["endpoint_blocks"])==1
        out.append(q)
    return out
def coefficient(rows,fields,group_mode):
    # Use identical response, weights and eligible pair support within each comparison.
    if len(rows)<100:return None
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        key=r["continent"] if group_mode=="continent" else (r["continent"],r["block"])
        groups[key].append(i)
    if len(groups)<2 and group_mode!="continent":return None
    if group_mode=="block" and sum(len(v)>1 for v in groups.values())<10:return None
    w=block_weights(rows)
    X=np.array([[r[k] for k in fields+TARGETS] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite original regression matrix")
    for ix in groups.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None]
    result={}
    for v in RAIL:
        j=fields.index(v);z=X[:,j];Y=X[:,-3:]
        other=np.delete(X[:,:-3],j,axis=1)
        scale=np.sqrt((other**2).sum(axis=0));ok=scale>1e-12
        other=other[:,ok]/scale[ok]
        rz=z-other@np.linalg.lstsq(other,z,rcond=1e-10)[0]
        ry=Y-other@np.linalg.lstsq(other,Y,rcond=1e-10)[0]
        var=float(rz@rz)
        if var<=1e-10:return None
        beta=rz@ry/var
        if abs(beta[0]-beta[1]-beta[2])>1e-10:raise ValueError("Component slopes not additive")
        result[v]={"coefficients":dict(zip(TARGETS,map(float,beta))),"residual_variance":var}
    return {"n_pairs":len(rows),"n_groups":len(groups),"n_pair_blocks":len({(r["continent"],r["block"]) for r in rows}),"direction":result}
def delete_endpoint_blocks(rows):
    blocks=sorted({b for r in rows for b in r["endpoint_blocks"]})
    v=[]
    for b in blocks:
        fit=coefficient([r for r in rows if b not in r["endpoint_blocks"]],FULL,"block")
        if fit is not None:v.append(fit["direction"][RAIL[0]]["coefficients"]["beta_jtu"])
    return {"n_endpoint_blocks":len(blocks),"n_valid":len(v),"min":min(v) if v else None,
            "max":max(v) if v else None,"positive":sum(x>0 for x in v),"negative":sum(x<0 for x in v)}
def panel(cells,neighbors,extra,coords,w,gate):
    old=pair_partition(neighbors,cells)
    n,p,_,_=EXPECTED[w,gate]
    if len(cells)!=n or len(old)!=p:raise ValueError("Original v386 pair denominator mismatch")
    d=direction(old)
    if d is None or abs(d[RAIL[0]]["beta_jtu"]-EXPECTED_TURNOVER[w,gate])>1e-9:
        raise ValueError("Original v386 turnover coefficient not reproduced")
    rows=add_pair_covariates(neighbors,cells,extra,coords)
    out={"width_m":w,"min_land_fraction":gate,"original_pairs":p,"complete_climate_pairs":len(rows),
         "complete_by_continent":dict(Counter(r["continent"] for r in rows)),"models":[]}
    # Primary full matched pair sample; sensitivity restricts both endpoints
    # to same block, so block FE is not an arbitrary pair midpoint grouping.
    for mode,r in (("all_complete_pairs",rows),("same_endpoint_block",[x for x in rows if x["both_endpoints_same_block"]])):
        fits={"original_covariates_continent_FE":coefficient(r,FIELDS,"continent"),
              "extended_climate_continent_FE":coefficient(r,FULL,"continent"),
              "extended_climate_block_FE":coefficient(r,FULL,"block")}
        fit=fits["extended_climate_block_FE"]
        per={c:coefficient([x for x in r if x["continent"]==c],FULL,"block") for c in CONTINENTS}
        cv=continent_holdout(r,PM["base_road"]+EXTRA,PM["both"]+EXTRA,"beta_jtu",True,
            MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,"v387 turnover extended climate")
        out["models"].append({"support":mode,"n_pairs":len(r),"fits":fits,
             "block_FE_by_continent":per,"block_FE_endpoint_deletion":delete_endpoint_blocks(r) if fit else None,
             "turnover_heldout":cv})
    print(json.dumps({"v387_panel":out},allow_nan=False),flush=True)
    return out
def main():
    p=argparse.ArgumentParser()
    for k in ("plant","rail","built","tree","road","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text())
    plant,rail,tree,road=map(read,(a.plant,a.rail,a.tree,a.road))
    checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree source")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or not road["version"].startswith("v381b"):raise ValueError("Wrong road source")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong original climate/map")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne)
    raster,mem,_=worldclim_tif(a.climate);out=[];qa={}
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
            if len(rows)!={100000:3455,200000:1442}[w]:raise ValueError("Original source denominator mismatch")
            extra,coords,q=load_extra_climate(a.climate,{r["key"] for r in rows},w);qa[str(w)]=q
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
            for gate in LAND_THRESHOLDS:
                m,_=enrich_cells(rows,rd,tr,w,gate);out.append(panel(m,neighbors,extra,coords,w,gate))
    finally:raster.close();mem.close()
    result={"version":"v387-original-v386-turnover-climate-and-local-block-retest",
            "source_checks":checks,"additional_source_climate_QA":qa,"panels":out,
            "limitations":["posthoc v386 same observational outcome, not independent field replication",
            "BIO4 and BIO15 use original raw units, BIO12 log1p raw; cell center macroclimate, not local heterogeneity",
            "block FE removes between-block differences, not all local geography, habitat, survey bias or human population",
            "source road/tree temporal mismatch; no direct mechanism, causality or actual railway flora",
            "endpoint-block deletion is sensitivity not CI or p-value"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
