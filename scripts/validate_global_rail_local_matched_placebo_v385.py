#!/usr/bin/env python3
"""v385 outcome-blind matched local rail exposure swap, source v382 exact replay."""
import argparse, json, math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from analyze_global_rail_road_tree_richness_v382 import (
    enrich_cells,make_lookup,pair_enrich,MODELS,PM,LAND_THRESHOLDS,within_continent_direction)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    verify_original_sources,source_climate_lookup,original_cell_rows,continent_holdout,
    GRID_WIDTHS,CONTINENTS,MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import read_source_builtup,inject_built,RICH_BUILT,ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
from validate_world_railkm_original_plant_placebo_v375 import tech_summary
FEATURES=("bio1","abs_lat","log_road_density",RICH_BUILT,"tree_consensus_fullcell","land_fraction","log_observations","log_observers")
DRAWS=99
SEED=38520261124
CALIPER=.75
EXPECTED={
 (100000,.5):(3350,4559,.0009221421712878227,.00011921337426956458),
 (100000,.8):(2959,3981,.000273406949724687,.0000806332565667083),
 (200000,.5):(1396,1894,.0021098589224810096,-.00005313450370118157),
 (200000,.8):(1171,1558,-.00004651742157152083,-.00009058569159828046)}
def match_pairs(rows):
    """Deterministic greedy matching in continent × 1000km block.
    Outcomes, rail km, plant species richness never enter matching.
    Per-continent SD scales, all-feature caliper .75 SD; no replacement.
    """
    cont=defaultdict(list);blocks=defaultdict(list)
    for i,r in enumerate(rows):
        cont[r["continent"]].append(i)
        blocks[(r["continent"],r["block"])].append(i)
    X=np.array([[r[k] for k in FEATURES] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Source nonfinite matching covariate")
    Z=X.copy()
    for ids in cont.values():
        scale=X[ids].std(axis=0);scale[scale<1e-12]=1
        Z[ids]=(X[ids]-X[ids].mean(axis=0))/scale
    out=[];distances=[]
    for ids in blocks.values():
        candidate=[]
        for k,i in enumerate(ids):
            for j in ids[k+1:]:
                d=np.abs(Z[i]-Z[j])
                if d.max()<=CALIPER:
                    candidate.append((float(d@d),min(i,j),max(i,j)))
        used=set()
        for d,i,j in sorted(candidate):
            if i in used or j in used:continue
            used.update((i,j));out.append((i,j));distances.append(d**.5)
    return out,{"matched_cells":2*len(out),"total_cells":len(rows),
        "matched_cells_by_continent":dict(Counter(rows[i]["continent"] for pair in out for i in pair)),
        "matched_blocks":len({(rows[i]["continent"],rows[i]["block"]) for i,j in out}),
        "mean_scaled_match_distance":float(np.mean(distances)) if distances else None,
        "max_scaled_match_distance":max(distances) if distances else None}
def swapped(rows,matches,rng):
    # Every swap uses original source values, never updates the real rows.
    out=[dict(r) for r in rows]
    changed=0
    for i,j in matches:
        if rng.integers(0,2)==1:
            for a,b in ((i,j),(j,i)):
                r=out[a];k=rows[b]["rail_km"]
                if r["rail_km"]!=k:changed+=1
                r["rail_km"]=k;r["log_rail_km"]=math.log1p(k)
                r["rail_tree"]=r["log_rail_km"]*r["tree_consensus_fullcell"]
    if sorted(r["rail_km"] for r in out)!=sorted(r["rail_km"] for r in rows):
        raise ValueError("Rail exposure multiset changed")
    return out,changed
def fit(rows,neighbors):
    pairs=pair_enrich(neighbors,rows)
    r=continent_holdout(rows,MODELS["base_road"],MODELS["both"],"log_species",False,
                        MIN_CONTINENT_RICH_CELLS,2*MIN_CONTINENT_RICH_CELLS,"v385 richness")
    b=continent_holdout(pairs,PM["base_road"],PM["both"],"jaccard",True,
                        MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,"v385 Jaccard")
    return r,b,pairs
def delta(m):return m["mean_holdout_continent_delta_MSE"]
def analyze_panel(rows,neighbors,w,g):
    n,p,r,b=EXPECTED[w,g]
    rr,bb,pairs=fit(rows,neighbors)
    if len(rows)!=n or len(pairs)!=p or abs(delta(rr)-r)>1e-9 or abs(delta(bb)-b)>1e-9:
        raise ValueError("Exact v382 road-adjusted real-map replay failed")
    matches,qa=match_pairs(rows)
    qa["pairs_with_distinct_rail_exposure"]=sum(rows[i]["rail_km"]!=rows[j]["rail_km"] for i,j in matches)
    qa["movable_distinct_rail_cells"]=2*qa["pairs_with_distinct_rail_exposure"]
    qa["movable_fraction"]=qa["movable_distinct_rail_cells"]/len(rows)
    # On the actual map, estimate sign of mean and difference rail
    # exposure separately; positive predictive gain alone is unsigned.
    directions={"richness_rail":within_continent_direction(rows,"log_species",MODELS["both"],"log_rail_km"),
        "Jaccard_mean_rail":within_continent_direction(pairs,"jaccard",PM["both"],"log_mean_rail_km"),
        "Jaccard_rail_difference":within_continent_direction(pairs,"jaccard",PM["both"],"absdiff_log_rail_km")}
    out={"width_m":w,"min_land_fraction":g,"matching":qa,"real_direction":directions,
         "real_road_adjusted_delta_MSE":{"richness":r,"jaccard":b},
         "original_eligible_holdouts":{"richness":rr["n_original_continents_eligible_as_heldout_folds"],
                                        "jaccard":bb["n_original_continents_eligible_as_heldout_folds"]}}
    if qa["movable_fraction"]<.20 or qa["matched_blocks"]<20:
        out["placebo_status"]="insufficient movable exposure or local blocks; no technical rank"
        return out
    rng=np.random.default_rng(SEED+w+int(100*g))
    vals={"richness":[],"jaccard":[]};changes=[]
    for i in range(DRAWS):
        fake,c=swapped(rows,matches,rng);changes.append(c)
        fr,fb,fp=fit(fake,neighbors)
        if len(fp)!=p:raise ValueError("Placebo altered original biological sample")
        if (fr["n_original_continents_eligible_as_heldout_folds"]!=rr["n_original_continents_eligible_as_heldout_folds"]
            or fb["n_original_continents_eligible_as_heldout_folds"]!=bb["n_original_continents_eligible_as_heldout_folds"]):
            raise ValueError("Placebo altered heldout support")
        vals["richness"].append(delta(fr));vals["jaccard"].append(delta(fb))
    out["placebo_status"]="technical local matched swaps, not spatial exchangeability proof"
    out["changed_cells_mean"]=float(np.mean(changes))
    out["changed_cells_min_max"]=[min(changes),max(changes)]
    out["conditional_placebo"]={k:tech_summary(v,vals[k]) for k,v in (("richness",r),("jaccard",b))}
    return out
def main():
    p=argparse.ArgumentParser()
    for k in ("plant","rail","built","tree","road","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text())
    plant,rail,tree,road=map(read,(a.plant,a.rail,a.tree,a.road))
    checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree source")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or not road["version"].startswith("v381b"):raise ValueError("Wrong road source")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong land/climate source")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne)
    raster,mem,_=worldclim_tif(a.climate);panels=[]
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
            if len(rows)!={100000:3455,200000:1442}[w]:raise ValueError("Original v376c source panel mismatch")
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
            for gate in LAND_THRESHOLDS:
                m,_=enrich_cells(rows,rd,tr,w,gate)
                result=analyze_panel(m,neighbors,w,gate);panels.append(result)
                print(json.dumps({"v385_panel":result},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    out={"version":"v385-local-environment-matched-rail-exposure-swap","source_checks":checks,
         "matching_features":FEATURES,"max_each_feature_continental_sd_caliper":CALIPER,"draws":DRAWS,"seed":SEED,
         "panels":panels,"limits":["rail exposure swapped within 1000km blocks and covariate caliper, not original geometry or rail route relocation",
         "no proof of exchangeability, residual spatial confounding possible",
         "same posthoc v382 outcomes, not independent validation","not causal field inference or actual rail habitat","no individual GPS or IDs"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
