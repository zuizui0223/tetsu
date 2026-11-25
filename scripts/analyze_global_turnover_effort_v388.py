#!/usr/bin/env python3
"""v388: observation-effort sensitivity of observed landscape turnover."""
import argparse,json,math
from pathlib import Path
from collections import Counter
from analyze_global_rail_beta_partition_v386 import pair_partition,direction,delete_blocks,EXPECTED,RAIL
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,PM,LAND_THRESHOLDS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
 verify_original_sources,source_climate_lookup,original_cell_rows,continent_holdout,
 GRID_WIDTHS,CONTINENTS,MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import (
 read_source_builtup,inject_built,ORIGINAL_NE_SHA)
from generate_global_equal_area_climate_land_grid_v351 import (
 WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file)
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
from analyze_global_turnover_climate_block_v387 import EXPECTED_TURNOVER
def close_ratio(a,b,max_ratio):
    if min(a,b)<=0:return False
    return max(a,b)/min(a,b)<=max_ratio

def support(a,b,mode):
    # Source endpoint effort and richness only, not observed Jaccard,
    # shared species, rail exposure or model residuals.
    if mode=="all":return True
    if mode in ("effort_balance","balance_richness","high_effort_balance"):
        if not (close_ratio(a["observations"],b["observations"],2.) and
                close_ratio(a["observers"],b["observers"],2.)):
            return False
    if mode=="balance_richness":
        return close_ratio(a["species"],b["species"],2.)
    if mode=="high_effort_balance":
        return min(a["observations"],b["observations"])>=200 and min(a["observers"],b["observers"])>=20
    if mode=="effort_balance":return True
    raise ValueError("Unknown frozen support")

def endpoint_rows(neighbors,cells):
    m={r["key"]:r for r in cells}
    rows=pair_partition(neighbors,cells)
    endpoints=[]
    for p in neighbors:
        a=tuple(p["source_cell_A_original_grid_ix_iy"]);b=tuple(p["source_cell_B_original_grid_ix_iy"])
        if a in m and b in m and m[a]["continent"]==m[b]["continent"]:
            endpoints.append((m[a],m[b]))
    if len(rows)!=len(endpoints):raise ValueError("Original pair support mismatch")
    return list(zip(rows,endpoints))
def panel(cells,neighbors,w,gate):
    ep=endpoint_rows(neighbors,cells)
    n,p,_,_=EXPECTED[w,gate]
    if len(cells)!=n or len(ep)!=p:raise ValueError("Original v386 panel changed")
    base=direction([r for r,_ in ep])
    if base is None or abs(base[RAIL[0]]["beta_jtu"]-EXPECTED_TURNOVER[w,gate])>1e-9:
        raise ValueError("Original v386 signed turnover coefficient not reproduced")
    out=[]
    for mode in ("all","effort_balance","balance_richness","high_effort_balance"):
        rows=[r for r,(a,b) in ep if support(a,b,mode)]
        d=direction(rows)
        blocks={b for r in rows for b in r["endpoint_blocks"]}
        if len(blocks)<20:d=None
        # Spatial support criterion independent of observed sign.
        by=Counter(r["continent"] for r in rows)
        cv=continent_holdout(rows,PM["base_road"],PM["both"],"beta_jtu",True,
                MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,"v388 turnover") if d is not None else None
        rec={"support":mode,"n_pairs":len(rows),"n_endpoint_blocks":len(blocks),
             "pairs_by_continent":dict(by),"direction":d,
             "direction_by_continent":{c:direction([r for r in rows if r["continent"]==c]) for c in CONTINENTS},
             "endpoint_block_deletion":delete_blocks(rows) if d is not None and mode!="all" else None,
             "turnover_heldout":cv}
        out.append(rec)
    result={"width_m":w,"min_land_fraction":gate,"original_pairs":p,"support_models":out}
    print(json.dumps({"v388_panel":result},allow_nan=False),flush=True)
    return result
def main():
    p=argparse.ArgumentParser()
    for k in ("plant","rail","built","tree","road","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text())
    plant,rail,tree,road=map(read,(a.plant,a.rail,a.tree,a.road))
    checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong source tree")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or not road["version"].startswith("v381b"):raise ValueError("Wrong source road")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong map/climate")
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
            if len(rows)!={100000:3455,200000:1442}[w]:raise ValueError("Original v376c source mismatch")
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
            for gate in LAND_THRESHOLDS:
                matched,_=enrich_cells(rows,rd,tr,w,gate)
                panels.append(panel(matched,neighbors,w,gate))
    finally:raster.close();mem.close()
    out={"version":"v388-original-v386-effort-balanced-turnover-sensitivity",
         "source_checks":checks,"panels":panels,
         "limits":["posthoc same original samples; conditional subset comparisons, not independent evidence",
         "equal record and observer counts do not imply equal detection or survey design",
         "richness balance may condition on an ecological response, sensitivity only",
         "no raw records, coverage estimators, observer overlap, individual histories or occupancy model",
         "block deletion ranges are not confidence intervals; no causal transport inference"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
