#!/usr/bin/env python3
"""v374: original GLOBAL railway line-length probability sample climate validation.

2048 true fixed 2026 Overture conventional-rail WGS84 points drawn with
replacement by original full-world GERS geodesic length v368b, exact source
summary/seed checks. WorldClim 2.1 BIO1, NE continents as in v351.
Compare DIRECT climate *at actual sampled rail point* with 100/200km grid
centres (v369/v370 original EPSG:6933 origin), and the independent v373
rail-km-by-100/200km LAND-centre classification.

This does NOT observe the botanical presence, seed carriage or microclimate.
Only aggregate results leave CI. No original point coordinates or GERS IDs.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median
from shapely.geometry import Point
from pyproj import Transformer

from generate_global_equal_area_climate_land_grid_v351 import (
    WORLDCLIM_SHA,CONTINENTS,source_continents,worldclim_tif,
    read_temp,climate_scale,sha_file,temp_band
)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
from sample_entire_global_GERS_railkm_WGS84_v368b import (
    SAMPLE_K,SEED,EXPECTED_CANDIDATE_KM,EXPECTED_CANDIDATE_SEGMENTS
)

SEED_SHA="0f78310226db6e4448fcf56467b4ac0bc38c185afa25569adee304aba2634d2b"
EXPECTED_V373_KM=2237043.481
WIDTHS=(100_000,200_000)
X_STOP=17_400_000
Y_STOP=7_400_000
BANDS=("cool","temperate","warm")
CATEGORIES=CONTINENTS+("unassigned_land_map_units",)

def valid_source_sample(obj):
    meta=obj.get("source_global_rail_sample_meta",{})
    points=obj.get("original_public_rail_point_draws_no_botanist_coordinates")
    if (meta.get("sample_seed_SHA256")!=SEED_SHA
        or SEED_SHA!=hashlib.sha256(SEED.encode()).hexdigest()
        or meta.get("n_independent_uniform_world_rail_km_draws_with_replacement")!=SAMPLE_K
        or not meta.get("original_128_global_source_parquet_complete_verified")
        or abs(meta.get("sampling_frame_original_fullworld_candidate_normal_rail_metres",0)/
               1000-EXPECTED_CANDIDATE_KM)>.05):
        raise ValueError("Not pinned original global independent 2048 rail-km PPS sample")
    if not isinstance(points,list) or len(points)!=SAMPLE_K:
        raise ValueError("Original source 2048 rail-km sampled points missing")
    idx=[r["original_independent_rail_km_point_draw_index"] for r in points]
    if sorted(idx)!=list(range(SAMPLE_K)):
        raise ValueError("Original source point draw IDs missing/duplicate")
    nc=Counter()
    for row in points:
        lon=row.get("source_original_2026_WGS84_track_point_lon")
        lat=row.get("source_original_2026_WGS84_track_point_lat")
        if (not isinstance(lon,(float,int)) or not isinstance(lat,(float,int))
            or not math.isfinite(lon) or not math.isfinite(lat)
            or not -180<=lon<=180 or not -90<=lat<=90):
            raise ValueError("Source point WGS84 lat/lon invalid")
        if row.get("source_site_geocode_is_synthetic_rail_km_point_NOT_plant_or_observer_GPS") is not True:
            raise ValueError("Expected independent railway centerline points, not citizen observations")
        nc[row.get("source_original_rail_class")]+=1
        if row.get("source_original_rail_class") not in ("standard_gauge","narrow_gauge","unknown"):
            raise ValueError("Atypical source rail class in conventional rail frame")
    if nc!=Counter(meta.get("source_sampled_original_rail_class_counts",{})):
        raise ValueError("Rail class source tally and whole 2048 geometry point source disagree")
    return points,meta

def grid_from_true_point(lon,lat,width,forward,inverse):
    x,y=forward.transform(lon,lat)
    if not all(math.isfinite(v) for v in (x,y)):return None
    if not (X_START<=x<X_STOP and Y_START<=y<Y_STOP):return None
    ix=math.floor((x-X_START)/width)
    iy=math.floor((y-Y_START)/width)
    center_x=X_START+(ix+.5)*width
    center_y=Y_START+(iy+.5)*width
    lon0,lat0=inverse.transform(center_x,center_y)
    if not all(math.isfinite(v) for v in (lon0,lat0)):return None
    if not -180<=lon0<=180 or not -90<=lat0<=90:return None
    return (int(ix),int(iy)),(float(lon0),float(lat0))

def land_source_continent(lat,lon,land):
    point=Point(lon,lat)
    hits=[name for name in CONTINENTS if land[name].covers(point)]
    if len(hits)!=1:return "unassigned_land_map_units"
    return hits[0]

def safe_band(v):
    return temp_band(v) if v is not None and math.isfinite(v) else None

def weighted_counts(rows,key,condition=lambda _:True):
    counts=Counter(r[key] for r in rows if r[key] is not None and condition(r))
    total=sum(counts.values())
    return {"n_valid":total,"band_counts":{band:counts.get(band,0) for band in BANDS},
            "band_fractions_valid":{band:(counts.get(band,0)/total if total else None)
                 for band in BANDS}}

def num_summary(nums):
    a=sorted(x for x in nums if x is not None and math.isfinite(x))
    if not a:return {"n":0,"mean":None,"median":None,"p10":None,"p90":None,
                     "min":None,"max":None}
    n=len(a)
    return {"n":n,"mean":mean(a),"median":median(a),
            "p10":a[int(.1*(n-1))],"p90":a[int(.9*(n-1))],
            "min":a[0],"max":a[-1]}

def source_v373(v373):
    d=v373.get("results",[])
    if not isinstance(d,list) or {x.get("source_grid_width_m") for x in d}!={100_000,200_000}:
        raise ValueError("Not authentic v373 global 100km/200km climate-band results")
    valid={}
    for row in d:
        q=row.get("source_full_global_climate_context_QA",{})
        if abs(q.get("original_full_global_rail_geodesic_km",0)-EXPECTED_V373_KM)>.1:
            raise ValueError("Original v373 rail km census mismatch")
        c=row.get("source_each_of_six_continents_climate_band",{})
        if set(c)!=set(CONTINENTS):raise ValueError("v373 original six continents missing")
        valid[row["source_grid_width_m"]]=row
    return valid

def analyze(original_points,original_land,worldclim,verified_v373):
    points_obj=json.loads(Path(original_points).read_text(encoding="utf8"))
    points,metadata=valid_source_sample(points_obj)
    prior=source_v373(json.loads(Path(verified_v373).read_text(encoding="utf8")))
    if sha_file(worldclim)!=WORLDCLIM_SHA:
        raise ValueError("Original pinned WorldClim BIO1 source SHA mismatch")
    land,landqa=source_continents(original_land)
    raster,mem,sourceqa=worldclim_tif(worldclim)
    forward=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    inverse=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    outrows=[]
    try:
        scale,scaleqa=climate_scale(raster)
        for r in points:
            lon=float(r["source_original_2026_WGS84_track_point_lon"])
            lat=float(r["source_original_2026_WGS84_track_point_lat"])
            cont=land_source_continent(lat,lon,land)
            t=read_temp(raster,lat,lon)
            t=(t/scale if t is not None else None)
            if t is not None and not -70<=t<=50:
                raise ValueError("Real railway location WorldClim physical BIO1 invalid")
            item={"continent":cont,"direct_BIO1_C":t,
                  "direct_temperature_band":safe_band(t),
                  "bridge_tunnel":bool(r.get("source_overpass_or_tunnel_candidate")),
                  "rail_class":r["source_original_rail_class"]}
            for width in WIDTHS:
                entry=grid_from_true_point(lon,lat,width,forward,inverse)
                if entry is None:
                    item[f"center_{width}_BIO1_C"]=None
                    item[f"center_{width}_band"]=None
                    item[f"center_{width}_continent"]=None
                    continue
                (ix,iy),(xc,yc)=entry
                item[f"center_{width}_continent"]=land_source_continent(yc,xc,land)
                tm=read_temp(raster,yc,xc)
                tm=(tm/scale if tm is not None else None)
                item[f"center_{width}_BIO1_C"]=tm
                item[f"center_{width}_band"]=safe_band(tm)
            outrows.append(item)
    finally:
        raster.close();mem.close()
    if len(outrows)!=SAMPLE_K:
        raise AssertionError("Lost original railway length-PPS locations")
    legitimate=[r for r in outrows if r["continent"] in CONTINENTS and
                r["direct_temperature_band"] is not None]
    if len(legitimate)<int(.85*SAMPLE_K):
        raise ValueError("Global source point-level land/climate coverage below 85%")
    continent_results={}
    for cont in CONTINENTS:
        selected=[r for r in outrows if r["continent"]==cont]
        valid=[r for r in selected if r["direct_temperature_band"] is not None]
        continent_results[cont]={
             "n_original_rail_km_PPS_points_total":len(selected),
             "n_direct_WorldClim_BIO1_valid":len(valid),
             "temperature_C":num_summary([x["direct_BIO1_C"] for x in valid]),
             "direct_source_temperature_band":weighted_counts(valid,"direct_temperature_band"),
             "n_original_bridge_tunnel_PPS_points":sum(x["bridge_tunnel"] for x in selected),
             "surface_only_band":weighted_counts(valid,"direct_temperature_band",
                                                  lambda x:not x["bridge_tunnel"])
        }
    global_temps=weighted_counts(legitimate,"direct_temperature_band")
    resolutions={}
    for width in WIDTHS:
        prefix=f"center_{width}"
        same=[r for r in legitimate if r[f"{prefix}_continent"] in CONTINENTS
                  and r[f"{prefix}_band"] is not None]
        same_cont=[r for r in same if r["continent"]==r[f"{prefix}_continent"]]
        diff=[r["direct_BIO1_C"]-r[f"{prefix}_BIO1_C"] for r in same_cont]
        discord=sum(r["direct_temperature_band"]!=r[f"{prefix}_band"] for r in same_cont)
        rows_v373=prior[width]["source_each_of_six_continents_climate_band"]
        known_km=sum(d["source_railkm_in_original_climate_band"]
                 for c in rows_v373.values() for d in c.values())
        km_by_band={band:sum(c[band]["source_railkm_in_original_climate_band"]
                             for c in rows_v373.values()) for band in BANDS}
        expected={band:km_by_band[band]/known_km for band in BANDS}
        resolutions[str(width)]={
            "original_rail_point_samples_with_both_railpoint_and_valid_grid_center":len(same),
            "n_railpoint_samples_where_center_is_same_continent":len(same_cont),
            "n_direct_source_rail_points_with_invalid_or_ocean_grid_center":len(legitimate)-len(same),
            "n_direct_vs_center_macroclimate_band_changed":discord,
            "fraction_direct_vs_center_macroclimate_band_changed":discord/len(same_cont) if same_cont else None,
            "direct_railpoint_minus_center_BIO1_C":num_summary(diff),
            "v373_original_fullrailkm_weighted_100km_or_200km_classified_km":known_km,
            "v373_original_fullrailkm_BIO1_classification_fractions":expected,
            "direct_source_2048_PPS_band_fractions_valid":global_temps["band_fractions_valid"],
            "PPS_minus_full_census_cellcenter_band_fractions":{
                band:global_temps["band_fractions_valid"][band]-expected[band] for band in BANDS},
            "actual_railpoint_vs_center_location_climate_not_railway_microclimate":True}
    return {
        "version":"v374-genuine-world-2048-independent-length-PPS-railpoint-WGS84-physical-climate",
        "source_original_Overture_release":metadata["original_source_Overture_release"],
        "source_true_GERS_global_rail_km_population":metadata["sampling_frame_original_fullworld_candidate_normal_rail_metres"]/1000,
        "source_length_PPS_sample_seed_sha256":SEED_SHA,
        "n_global_original_independent_railkm_samples_with_replacement":SAMPLE_K,
        "n_valid_land_and_WorldClim_railkm_samples":len(legitimate),
        "n_naturalearth_unassigned_original_railkm_samples":sum(x["continent"]=="unassigned_land_map_units" for x in outrows),
        "n_known_continent_but_no_BIO1_samples":sum(x["continent"] in CONTINENTS and x["direct_BIO1_C"] is None for x in outrows),
        "sample_source_rail_classes":dict(Counter(r["rail_class"] for r in outrows)),
        "sample_original_bridge_tunnel_candidates":sum(r["bridge_tunnel"] for r in outrows),
        "source_original_Natural_Earth_map_quality":landqa,
        "source_original_WorldClim_BIO1_quality":sourceqa,
        "source_original_WorldClim_BIO1_temperature_scale":scaleqa,
        "source_global_direct_railkm_PPS_BIO1_C":num_summary([x["direct_BIO1_C"] for x in legitimate]),
        "source_global_direct_railkm_PPS_temperature_bands":global_temps,
        "source_six_original_continent_PPS_sampled_temperature":continent_results,
        "source_true_railpoint_vs_coarse_v373_grid_cell_center":resolutions,
        "scientific_limits":[
            "These are physical GLOBAL railway centreline positions not vegetation, habitat quality or live train service",
            "Overture rail line multiple tracks/overlaps could be counted as multiple geometry length, standard-gauge includes unverified live operation",
            "WorldClim2.1 1970-2000 BIO1 is MACRO climate and never measured ballast surface or species thermal preference",
            "Random rail-km points are sampled with replacement; 2048 rail locations are not 2048 independent ecological communities",
            "Nearest grid centre climate may miss coast, island, true mountainous microclimate and river gorge",
            "Unclassified original 2048 rail-km points remain unclassified, no favorable replacement sampling",
            "Source land map 110m and 100/200km gridding use genuinely different spatial support, difference is resolution audit",
            "Source 2025 iNaturalist original photo-linked whole-landscape plant taxon data not used in this physical climate comparison",
            "No geographic original points, exact railway UUID, botanical observer IDs or individual rare plant GPS published"]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-2048-railkm-points",required=True)
    p.add_argument("--natural-earth",required=True)
    p.add_argument("--worldclim-zip",required=True)
    p.add_argument("--verified-v373-global-grid-climate",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=analyze(a.original_2048_railkm_points,a.natural_earth,a.worldclim_zip,
                  a.verified_v373_global_grid_climate)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf8")
    compact={k:v for k,v in result.items()
         if k in ("version","n_global_original_independent_railkm_samples_with_replacement",
             "n_valid_land_and_WorldClim_railkm_samples",
             "n_naturalearth_unassigned_original_railkm_samples",
             "source_global_direct_railkm_PPS_BIO1_C",
             "source_global_direct_railkm_PPS_temperature_bands",
             "source_six_original_continent_PPS_sampled_temperature",
             "source_true_railpoint_vs_coarse_v373_grid_cell_center")}
    print(json.dumps(compact,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
