#!/usr/bin/env python3
"""v370 TRUE complete global rail-km density per 100km and 200km equal-area
grid, from exact original full 16 SHA-pinned Overture 2026 rail Parquets.

This is PHYSICAL RAILWAY INFRASTRUCTURE exposure (km in original grid),
not active train frequency, trackbed floral community, native species
range, seed movement or railway-caused vegetation homogenization.

All source GIS geometry/ID stays transient in original GitHub Actions
and the output is 100km/200km rail-km aggregates ONLY.
"""
from __future__ import annotations
import argparse,json,math,time
from collections import Counter,defaultdict
from pathlib import Path

from build_global_Overture_physical_connector_components_v361 import (
    validate_files,SOURCE_RELEASE
)
from sample_entire_global_GERS_railkm_WGS84_v368b import (
    stream_source_rails,candidate_segment,seg_m,
    EXPECTED_CANDIDATE_SEGMENTS,EXPECTED_CANDIDATE_KM
)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    X_START,X_END,Y_START,Y_END,CELL_SIZE_M
)

LIMIT_ORIGINAL_CROSSED_CELLS_PER_LINE=500
MAX_ALLOWED_MAPPED_GERS_RAIL_GRID_KM_MISMATCH=2.
SOURCE_EXPECTED_GRID_N=(X_END-X_START)//CELL_SIZE_M,(Y_END-Y_START)//CELL_SIZE_M

def tile_index(x,y):
    if not math.isfinite(x) or not math.isfinite(y):
        return None
    if X_START<=x<X_END and Y_START<=y<Y_END:
        return (int((x-X_START)//CELL_SIZE_M),
                int((y-Y_START)//CELL_SIZE_M))
    return None

def clip_source_line_to_tiles(source_geom,original_m,forward,inverse,geod):
    from shapely.geometry import box
    from shapely.ops import transform
    coords=list(source_geom.coords)
    if len(coords)<2:raise ValueError("GERS rail line not two original WGS vertices")
    lons=[v[0] for v in coords]
    if max(lons)-min(lons)>180:
        # Antimeridian discontinuity in EPSG:6933 map projection: never
        # draw a fake source railway across 35,000km of world grid.
        return {},original_m,"source_dateline_crossing"
    proj=transform(forward.transform,source_geom)
    bx=proj.bounds
    a=tile_index(bx[0],bx[1]);b=tile_index(bx[2],bx[3])
    # grid tile at exact upper bound potentially =grid outer boundary;
    # outside global source map is treated as unassignable not silently zero.
    if a is None or b is None:
        return {},original_m,"source_outside_100km_world_grid"
    if a==b:
        return {a:original_m},0.0,"source_one_tile"
    nx=b[0]-a[0]+1
    ny=b[1]-a[1]+1
    if nx*ny>LIMIT_ORIGINAL_CROSSED_CELLS_PER_LINE:
        raise ValueError("Original long GERS line crosses too many global 100km tiles")
    result={}
    for i in range(a[0],b[0]+1):
        for j in range(a[1],b[1]+1):
            west=X_START+i*CELL_SIZE_M
            south=Y_START+j*CELL_SIZE_M
            cut=proj.intersection(box(west,south,west+CELL_SIZE_M,south+CELL_SIZE_M))
            if cut.is_empty:continue
            def segments(g):
                if g.geom_type=="LineString" and len(g.coords)>=2 and g.length>0:
                    yield g
                elif hasattr(g,"geoms"):
                    for v in g.geoms:yield from segments(v)
            m=0.
            for part in segments(cut):
                m+=abs(float(geod.geometry_length(transform(inverse.transform,part))))
            if m>0:result[i,j]=m
    delta=original_m-sum(result.values())
    if abs(delta)>max(.5,.008*original_m):
        # Accept only small finite WGS84 projection/chord approximation
        # differences, fail on source dateline/map geometry artifacts.
        raise ValueError("Original GERS railway within-grid geodesic length unaccounted beyond 0.8%")
    return result,0.0,"source_multitile_clip"

def classify_world_source(paths):
    from pyproj import Transformer,Geod
    trans=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    inv=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    geod=Geod(ellps="WGS84")
    cells=defaultdict(lambda:{"source_rail_length_metres":0.,
        "original_source_GERS_segmentation_instances":0,
        "source_class_segment_intersection_counts":Counter()})
    qa=Counter();source_total_m=0.;excluded_source_m=0.
    started=time.monotonic()
    for row in stream_source_rails(paths):
        qa["source_all_original_GERS_rail_segmentation_rows"]+=1
        if not candidate_segment(row):continue
        geom,L=seg_m(row,geod)
        source_total_m+=L
        qa["source_normal_active_candidate_rail_original_source_segments"]+=1
        assignments,omitted,status=clip_source_line_to_tiles(geom,L,trans,inv,geod)
        qa[status]+=1
        excluded_source_m+=omitted
        for key,length in assignments.items():
            entry=cells[key]
            entry["source_rail_length_metres"]+=length
            entry["original_source_GERS_segmentation_instances"]+=1
            entry["source_class_segment_intersection_counts"][row["class"]]+=1
        if qa["source_normal_active_candidate_rail_original_source_segments"]%200000==0:
            print(json.dumps({"original_global_source_rail_segments_processed":
                qa["source_normal_active_candidate_rail_original_source_segments"],
                "source_total_original_rail_km":round(source_total_m/1000,3),
                "source_elapsed_seconds":round(time.monotonic()-started,1)}),flush=True)
    if qa["source_all_original_GERS_rail_segmentation_rows"]!=2054855:
        raise ValueError("Full original 128 Parquet all-rail GERS source count changed")
    if qa["source_normal_active_candidate_rail_original_source_segments"]!=EXPECTED_CANDIDATE_SEGMENTS:
        raise ValueError("v361 original world normal rail candidate segment count drift")
    if abs(source_total_m/1000-EXPECTED_CANDIDATE_KM)>.04:
        raise ValueError("Original 2026 full normal rail geography geodesic km differ from v361 census")
    accounted=sum(v["source_rail_length_metres"] for v in cells.values())
    d=source_total_m-accounted-excluded_source_m
    if abs(d)/1000>MAX_ALLOWED_MAPPED_GERS_RAIL_GRID_KM_MISMATCH:
        raise ValueError("Global rail-km GIS cell assignment lost >2km of genuine source infrastructure")
    return cells,qa,source_total_m,excluded_source_m,d

def rollup(grid,cell_m):
    # 100km x/y origin is multiple of 200km; exact 2x2 tile collapse
    grouped=defaultdict(lambda:{"rail_m":0.,"seg_n":0})
    factor=cell_m//CELL_SIZE_M
    for (ix,iy),cell in grid.items():
        p=grouped[ix//factor,iy//factor]
        p["rail_m"]+=cell["source_rail_length_metres"]
        p["seg_n"]+=cell["original_source_GERS_segmentation_instances"]
    out=[]
    for (ix,iy),cell in sorted(grouped.items()):
        out.append({"source_equal_area_grid_ix":ix,
                    "source_equal_area_grid_iy":iy,
                    "source_grid_width_m":cell_m,
                    "source_real_original_Overture_normal_rail_km_clipped_into_cell":round(cell["rail_m"]/1000.,5),
                    "source_original_GERS_segmentation_cell_overlap_n_NOT_train_routes":cell["seg_n"]})
    return out

def run(rail_dir,census):
    source=json.loads(Path(census).read_text(encoding="utf8"))
    paths=validate_files(rail_dir,source)
    grid,qa,L,missing,delta=classify_world_source(paths)
    aggregations={"100000":rollup(grid,100000),"200000":rollup(grid,200000)}
    totals={k:round(sum(x["source_real_original_Overture_normal_rail_km_clipped_into_cell"]
        for x in val),3) for k,val in aggregations.items()}
    return {"version":"v370-original-all-128-world-rail-2026-100km-200km-geodesic-km-exposure",
        "source_Oveture_release":SOURCE_RELEASE,
        "source_original_fullworld_rail_128_files_SHA_verified":True,
        "source_fullworld_candidate_normal_rail_km":round(L/1000,3),
        "original_source_rail_spatially_unassigned_dateline_or_map_border_km":round(missing/1000,5),
        "original_projection_clipping_residual_metres":round(delta,3),
        "source_true_global_rail_geographic_grid_QA":dict(qa),
        "source_global_gridded_rail_exposure_km_by_cell_width":aggregations,
        "source_100km_200km_grid_rail_km_total_within_mapped_cells":totals,
        "rail_source_measure_is_2026_geographic_centerline_km_NOT_operated_train_service_or_ballast_plants":True,
        "taxon_community_not_joined_or_calculated_here":True,
        "no_original_GERS_source_railway_coordinates_or_GERS_UUIDs_republished":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--rail-dir",required=True)
    p.add_argument("--original-v357-census",required=True)
    p.add_argument("--output-summary",required=True)
    p.add_argument("--output-grid",required=True)
    a=p.parse_args()
    obj=run(a.rail_dir,a.original_v357_census)
    out=Path(a.output_grid);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    compact={k:v for k,v in obj.items() if k!="source_global_gridded_rail_exposure_km_by_cell_width"}
    ss=Path(a.output_summary);ss.parent.mkdir(parents=True,exist_ok=True)
    ss.write_text(json.dumps(compact,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(compact,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
