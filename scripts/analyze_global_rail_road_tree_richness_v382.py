#!/usr/bin/env python3
"""v382: independent road/land/tree controls on original v369 worldwide plant responses."""
import argparse, json, math
from pathlib import Path
from collections import Counter
import numpy as np
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    verify_original_sources,source_climate_lookup,original_cell_rows,pair_rows,continent_holdout,
    fit_predict,metrics,block_weights,RICH_BASE,BETA_BASE,GRID_WIDTHS,CONTINENTS,
    MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import (
    read_source_builtup,inject_built,RICH_BUILT,ORIGINAL_NE_SHA)
from generate_global_equal_area_climate_land_grid_v351 import (
    WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file)
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5

LAND_THRESHOLDS=(.5,.8)
MIN_ROAD_LAND_SUPPORT=.95
MIN_TREE_GRID_SUPPORT=.95
def read(path):return json.loads(Path(path).read_text())
def make_lookup(rows):
    d={}
    for x in rows:
        k=(x["source_equal_area_grid_ix"],x["source_equal_area_grid_iy"])
        if k in d:raise ValueError("Duplicate original source grid row")
        d[k]=x
    return d

def enrich_cells(cells,roads,trees,w,min_land):
    out=[];qa=Counter()
    for r in cells:
        k=r["key"];rd=roads.get(k);tr=trees.get(k)
        if rd is None or tr is None:
            qa["missing_source_road_or_tree"]+=1;continue
        land=float(rd["source_GRIP4_land_fraction_of_cell"])
        valid=float(rd["source_GRIP4_road_coverage_of_source_land"])
        tv=float(tr["source_valid_area_fraction_of_grid_approx"])
        if not 0<land<=1.05 or not 0<=valid<=1.0000001 or not 0<tv<=1.05:
            raise ValueError("Source coverage range invalid")
        if land<min_land:qa["insufficient_land_area"]+=1;continue
        if valid<MIN_ROAD_LAND_SUPPORT:qa["insufficient_known_road_land"]+=1;continue
        if tv<MIN_TREE_GRID_SUPPORT:qa["insufficient_tree_raster_area"]+=1;continue
        density=rd["source_GRIP4_observed_road_density_m_per_km2_land"]
        tree=float(tr["source_tree_consensus_area_m2"])/(w*w)
        if density is None or not math.isfinite(density) or density<0 or not 0<=tree<=1.05:
            raise ValueError("Eligible original predictor invalid")
        r=dict(r,land_fraction=land,log_road_density=math.log1p(density),
                 tree_consensus_fullcell=tree)
        r["rail_tree"]=r["log_rail_km"]*tree
        r["road_tree"]=r["log_road_density"]*tree
        out.append(r)
    return out,dict(qa)
BASE=RICH_BASE+(RICH_BUILT,"land_fraction","tree_consensus_fullcell")
MODELS={
 "base":BASE,
 "base_rail":BASE+("log_rail_km",),
 "base_road":BASE+("log_road_density",),
 "both":BASE+("log_road_density","log_rail_km"),
 "road_interaction":BASE+("log_road_density","log_rail_km","road_tree"),
 "both_interactions":BASE+("log_road_density","log_rail_km","road_tree","rail_tree")}
CONTRASTS=(("base","base_rail"),("base","base_road"),("base_road","both"),
           ("base_rail","both"),("both","road_interaction"),("road_interaction","both_interactions"))
def heldout(rows,target,models,contrasts,clip,min_test,min_train):
    out={}
    for a,b in contrasts:
        out[a+"__to__"+b]=continent_holdout(rows,models[a],models[b],target,clip,
            min_test,min_train,target+" "+a+" vs "+b)
    return out
def within_continent_direction(rows,target,models,interaction_name):
    # Direction, NOT causal coefficients or nominal independent-cell p-values.
    # Fit with continent FE via weighted demeaning, on complete cases.
    fields=models
    if len(rows)<100:return {"n":len(rows),"estimate":None}
    w=block_weights(rows)
    X=np.array([[r[k] for k in fields]+[r[target]] for r in rows],float)
    for c in CONTINENTS:
        ix=np.array([i for i,r in enumerate(rows) if r["continent"]==c],int)
        if ix.size:X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None]
    k=fields.index(interaction_name)
    z=X[:,k];Y=X[:,-1];other=np.delete(X[:,:-1],k,axis=1)
    def res(v):return v-other@np.linalg.lstsq(other,v,rcond=1e-10)[0]
    rz,ry=res(z),res(Y)
    denom=float(rz@rz)
    if denom<=1e-12:return {"n":len(rows),"estimate":None,"residual_variance":denom}
    return {"n":len(rows),"estimate":float(rz@ry/denom),"residual_variance":denom,
        "n_1000km_blocks":len({(r["continent"],r["block"]) for r in rows})}
def pair_enrich(neighbors,cells):
    m={r["key"]:r for r in cells};pairs,_=pair_rows(neighbors,m,{})
    keys=[]
    for p in neighbors:
        a=tuple(p["source_cell_A_original_grid_ix_iy"]);b=tuple(p["source_cell_B_original_grid_ix_iy"])
        if a in m and b in m and m[a]["continent"]==m[b]["continent"]:keys.append((a,b))
    if len(keys)!=len(pairs):raise ValueError("Original adjacent pair correspondence mismatch")
    for row,(a,b) in zip(pairs,keys):
        l,r=m[a],m[b]
        for v in (RICH_BUILT,"land_fraction","tree_consensus_fullcell","log_road_density"):
            row["mean_"+v]=(l[v]+r[v])/2
            row["absdiff_"+v]=abs(l[v]-r[v])
        row["rail_tree"]=row["log_mean_rail_km"]*row["mean_tree_consensus_fullcell"]
        row["road_tree"]=row["mean_log_road_density"]*row["mean_tree_consensus_fullcell"]
    return pairs
PB=BETA_BASE+tuple(k+v for v in (RICH_BUILT,"land_fraction","tree_consensus_fullcell") for k in ("mean_","absdiff_"))
RR=("mean_log_road_density","absdiff_log_road_density")
RL=("log_mean_rail_km","absdiff_log_rail_km")
PM={"base":PB,"base_rail":PB+RL,"base_road":PB+RR,"both":PB+RR+RL,
    "road_interaction":PB+RR+RL+("road_tree",),"both_interactions":PB+RR+RL+("road_tree","rail_tree")}
def main():
    p=argparse.ArgumentParser()
    for name in ("plant","rail","built","tree","road","ne","climate","out"):p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();plant,rail,tree,road=map(read,(a.plant,a.rail,a.tree,a.road))
    checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:
        raise ValueError("Tree data not source verified v380b")
    if road["source_zip_MD5"]!=MD5 or not road["version"].startswith("v381b"):
        raise ValueError("Road data not land-weighted v381b original source")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:
        raise ValueError("Original NE / WorldClim source mismatch")
    built,_=read_source_builtup(a.built)
    land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    out=[]
    try:
        scale,_=climate_scale(raster)
        for width in GRID_WIDTHS:
            g=plant["source_grids"][str(width)]
            original=g["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            neighbor=g["source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"]
            keys={(x["source_equal_area_grid_ix"],x["source_equal_area_grid_iy"]) for x in original}
            clim,_=source_climate_lookup(land,raster,scale,keys,width)
            rows,_,_=original_cell_rows(original,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)],clim)
            expected={100000:3455,200000:1442}[width]
            if len(rows)!=expected:raise ValueError("Original v376c cell denominator not reproduced")
            rows,_=inject_built(rows,built[str(width)])
            if len(rows)!=expected:raise ValueError("Original GHSL full panel not reproduced")
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
            for minland in LAND_THRESHOLDS:
                matched,qa=enrich_cells(rows,rd,tr,width,minland)
                pairs=pair_enrich(neighbor,matched)
                rich=heldout(matched,"log_species",MODELS,CONTRASTS,False,MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2)
                beta=heldout(pairs,"jaccard",PM,CONTRASTS,True,MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS)
                direction={}
                for target,table,fields in (("log_species",matched,MODELS),("jaccard",pairs,PM)):
                    direction[target]={x:within_continent_direction(table,target,fields["both_interactions"],x) for x in ("rail_tree","road_tree")}
                    direction[target]["per_continent"]={c:{x:within_continent_direction([r for r in table if r["continent"]==c],target,fields["both_interactions"],x) for x in ("rail_tree","road_tree")} for c in CONTINENTS}
                panel={"width_m":width,"min_land_fraction":minland,"original_eligible":expected,
                       "n_matched":len(matched),"n_pairs":len(pairs),"excluded":qa,
                       "matched_continents":dict(Counter(r["continent"] for r in matched)),
                       "richness_prediction":rich,"jaccard_prediction":beta,"interaction_directions":direction}
                out.append(panel)
                print(json.dumps({"panel":{k:v for k,v in panel.items() if k not in ("richness_prediction","jaccard_prediction")},
                   "richness_delta":{k:v["mean_holdout_continent_delta_MSE"] for k,v in rich.items()},
                   "jaccard_delta":{k:v["mean_holdout_continent_delta_MSE"] for k,v in beta.items()}},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"v382-original-world-plant-richness-and-Jaccard-tree-road-falsification",
        "source_replay":checks,"panels":out,"scope":"whole landscape observed plant richness/Jaccard, NOT GloNAF fraction or actual rail flora",
        "limitations":["historical GRIP4 roads and tree consensus vs 2025 plants and 2026 rails",
          "tree share of whole grid area, not land-only forest cover; land fraction covariate plus inland sensitivity",
          "road inventory incomplete and human accessibility confounding remains",
          "not causal dispersal, competition or forest invasion resistance",
          "descriptive direction estimates, no independent-cell p-values"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
