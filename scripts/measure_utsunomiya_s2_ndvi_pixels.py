#!/usr/bin/env python3
"""REAL Sentinel-2 L2A red/NIR/SCL pixel pilot at the 2023 Utsunomiya LRT.

Actual pixel values are read as tiny HTTP-range requests from public COGs.
A small outcome-blind set of fixed transects along independently verified
2023 N05 passenger line geometry is compared to offsets 1.6 km away.
NO raw licensed N05 line geometries, pixel locations or maps are released.

Two BOA calibration interpretations are retained because Earth Search 2024+
reports 'earthsearch:boa_offset_applied=true' while raster band metadata
still reports '-0.1' offset. We do NOT designate one as the truth or
declare a harmonized causal railway effect. The fixed scene dates in
2017 versus 2024 differ by about a month; phenology is a confound.

This first pilot can fail due to remote S3 COG HTTP Range restrictions;
a failed/partial source is NOT a measured zero NDVI.
"""
from __future__ import annotations
import argparse
import collections
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.parse import urlparse
import zipfile

import rasterio
from pyproj import Transformer
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

from audit_n05_historical_railway import (
    bounded_zip_download, shape_reader, known_year, is_current_end
)
from audit_utsunomiya_lrt_photo_opportunity import (
    ROUTE, COMPANY, UTM_EPSG, OPENING_DATE
)
from audit_s2_utsunomiya_calibration import STAC,asset_meta
from audit_utsunomiya_sentinel_stac import BBOX
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.sentinel_ndvi import (
    CLEAR_SCL_PRIMARY,calculate_ndvi,ndvi_branch_values,metadata_ambiguous,
    paired_annual_ndvi,prepost_paired_change
)
import requests

# Fixed before reading any S2 ground pixels. Scenes are from source-only
# calibration run #42440904218; same-day two adjoining MGRS tiles per year.
SCENES={
    2017:("S2A_54SVF_20170707_0_L2A","S2A_54SUF_20170707_0_L2A"),
    2024:("S2B_54SVF_20240804_0_L2A","S2B_54SUF_20240804_0_L2A"),
    2025:("S2B_54SVF_20250730_0_L2A","S2B_54SUF_20250730_0_L2A"),
}
ANCHORED_TRANSECTS=9
NEAR_OFFSET_M=140.
CONTROL_OFFSET_M=1600.
WINDOW_FRAGMENT_SCL_VALID=(4,5)
MAX_IO_TIMEOUT_MINUTES=19
STAC_ITEM_BASE="https://earth-search.aws.element84.com/v1/collections/sentinel-2-l2a/items/"
HTTPS_ACCEPTED_HOST_SUFFIXES=("amazonaws.com",)


def get_lrt_lines(raw_n05):
    epsg=Transformer.from_crs(4326,UTM_EPSG,always_xy=True)
    lines=[]
    count=0
    with zipfile.ZipFile(io.BytesIO(raw_n05)) as z:
        reader,_=shape_reader(z,"RailroadSection2")
        fields=[x[0] for x in reader.fields[1:]]
        for shape_row in reader.iterShapeRecords():
            d=dict(zip(fields,shape_row.record))
            if (str(d.get("N05_002") or "").strip()!=ROUTE or
                str(d.get("N05_003") or "").strip()!=COMPANY or
                known_year(d.get("N05_004"))!=2023 or
                not is_current_end(d.get("N05_005e"))):
                continue
            count+=1
            shape=shape_row.shape
            breaks=list(shape.parts)+[len(shape.points)]
            for a,b in zip(breaks,breaks[1:]):
                if b-a<2: continue
                line=LineString([
                    epsg.transform(lon,lat) for lon,lat in shape.points[a:b]
                ])
                if line.length>1: lines.append(line)
        reader.close()
    if count<1 or not lines or sum(l.length for l in lines)<10000:
        raise ValueError("LRT 2023 source no longer identifies >10km surface railway")
    return lines, count


def build_fixed_pairs(lines, n_anchors=ANCHORED_TRANSECTS):
    """Generate matched 140m/1600m points on both sides of a route.

    Both sides and geography frozen independently of NDVI. The far
    point must be >=1km from ANY segment of the railway line. All pixels
    reflect a single correlated urban tram corridor.
    """
    longest=max(lines,key=lambda x:x.length)
    union=unary_union(lines)
    points={}
    pairs=[]
    anchor_meters=[]
    for i in range(n_anchors):
        f=(i+0.5)/n_anchors
        d=f*longest.length
        c=longest.interpolate(d)
        prev=longest.interpolate(max(0.,d-30.))
        nxt=longest.interpolate(min(longest.length,d+30.))
        dx,dy=(nxt.x-prev.x,nxt.y-prev.y)
        norm=math.hypot(dx,dy)
        if norm<=0:
            continue
        nx,ny=(-dy/norm,dx/norm)
        for sign in (-1,1):
            center=Point(c.x,c.y)
            near=Point(c.x+sign*nx*NEAR_OFFSET_M,
                       c.y+sign*ny*NEAR_OFFSET_M)
            far=Point(c.x+sign*nx*CONTROL_OFFSET_M,
                      c.y+sign*ny*CONTROL_OFFSET_M)
            if near.distance(union)>250 or far.distance(union)<1000:
                continue
            near_id=f"T{i:02d}_{'L' if sign<0 else 'R'}_near"
            far_id=f"T{i:02d}_{'L' if sign<0 else 'R'}_control"
            points[near_id]=(near.x,near.y)
            points[far_id]=(far.x,far.y)
            pairs.append((near_id,far_id))
        anchor_meters.append(d)
    if len(pairs)<4:
        raise ValueError("Too few geographically valid transect pairs")
    return points,pairs,{
        "planned_anchors":n_anchors,
        "near_offset_m":NEAR_OFFSET_M,
        "control_offset_m":CONTROL_OFFSET_M,
        "n_fixed_pairs":len(pairs),
        "n_fixed_measurement_pixels":len(points),
        "sampling_unit":"fixed pixel points within ONE correlated LRT corridor",
        "not_independent_railway_treatments":True,
    }


def fetch_stac_scene(client,item_id):
    response=client.get(STAC_ITEM_BASE+item_id,timeout=(20,60))
    response.raise_for_status()
    data=response.json()
    if data.get("id")!=item_id or data.get("type")!="Feature":
        raise ValueError("STAC item changed or not a valid feature")
    assets=data.get("assets") or {}
    for key in ("red","nir","scl"):
        if key not in assets:
            raise ValueError(f"STAC item {item_id} has no {key} COG")
        url=str(assets[key].get("href") or "")
        parts=urlparse(url)
        if parts.scheme!="https" or not parts.hostname or not any(
            parts.hostname==suffix or parts.hostname.endswith("."+suffix)
            for suffix in HTTPS_ACCEPTED_HOST_SUFFIXES
        ) or not parts.path.lower().endswith(".tif"):
            raise ValueError(f"Unexpected COG asset URL for {item_id}/{key}")
    return data


def pixel_samples(asset,coordinates):
    """Read local 1-pixel observations from remote Cloud Optimized GeoTIFF.

    rasterio performs HTTP range reads rather than downloading entire
    100km imagery tile. SCL 20m categorical; red and NIR 10m DN.
    """
    href=asset["href"]
    with rasterio.open(href) as src:
        if src.count!=1 or src.crs is None or src.width<1000:
            raise ValueError("Unexpected imagery band raster header")
        src_epsg=src.crs.to_epsg()
        if src_epsg not in (32654,):
            raise ValueError(f"Unexpected Sentinel tile CRS {src_epsg}")
        # Exact same projected N05 transect coordinates across all years.
        coords=[coordinates[key] for key in coordinates]
        vals={}
        for key,pixel in zip(coordinates,src.sample(coords,indexes=1,masked=True)):
            sample=pixel[0]
            if getattr(sample,"mask",False):
                vals[key]=None
                continue
            try:
                value=float(sample)
            except (TypeError,ValueError):
                value=float("nan")
            vals[key]=value if math.isfinite(value) else None
        return vals


def read_scene_pixels(item,points):
    all_data={}
    for band in ("red","nir","scl"):
        all_data[band]=pixel_samples(item["assets"][band],points)
    return all_data


def choose_year_pixel_items(year,points,items):
    """Choose first quality-valid tile for every fixed point per year.

    Each scene pair is acquired on the same date. Tie-break is fixed tile
    identifier order (SVF then SUF), not observed vegetation or NDVI.
    """
    results={}
    quality=collections.Counter()
    for item in items:
        observed=read_scene_pixels(item,points)
        for point_id in points:
            if point_id in results:
                continue
            red=observed["red"].get(point_id)
            nir=observed["nir"].get(point_id)
            scl=observed["scl"].get(point_id)
            if red is None or nir is None or scl is None:
                quality["off_tile_or_nodata"]+=1
                continue
            scl_int=int(scl)
            if scl_int not in CLEAR_SCL_PRIMARY:
                quality[f"scl_rejected_{scl_int}"]+=1
                continue
            vals=ndvi_branch_values(red,nir,item)
            if vals["provider_flag_branch"] is None or vals["raster_offset_branch"] is None:
                quality["radiometry_values_invalid_or_missing"]+=1
                continue
            results[point_id]={
                "provider_flag_branch":vals["provider_flag_branch"],
                "raster_offset_branch":vals["raster_offset_branch"]
            }
            quality["success_paired_pixel_source"]+=1
    return results,dict(quality)


def run():
    raw=bounded_zip_download()
    lines,n_lines=get_lrt_lines(raw)
    points,pairs,geo=build_fixed_pairs(lines)
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-lrt-public-COPERNICUS-NDVI/1.2"})
    output={
        "run_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_N05_SHA256":hashlib.sha256(raw).hexdigest(),
        "source_N05_license":"noncommercial-only, source geometry not redistributed",
        "rail_operator":COMPANY,
        "operator_opening_date":OPENING_DATE,
        "source_N05_qualified_line_rows":n_lines,
        "sampling_design":geo,
        "stac_source":STAC_ITEM_BASE,
        "SCL_primary_classes":sorted(CLEAR_SCL_PRIMARY),
        "selection":"scenes frozen by independent 2017/2024/2025 STAC calibration audit; not chosen by NDVI",
        "cross_era_processing_baseline_contradiction_expected":True,
        "calendar":{},
        "prepost_paired_differences":{},
        "errors":{},
        "scientific_disclaimers":[
            "Sentinel measured NDVI is a spectral greenness proxy, NOT weed species, abundance, flower/spore ROI, seed transport, plant height, regrowth mechanism or heritable wind adaptation.",
            "Same urban LRT corridor supplies all transects; pixels and offsets are not independent rail experiments. Far points are NOT road-matched.",
            "2017 before construction; 2024/25 after service, while 2018-22 urban construction is an independent environmental event.",
            "2017 scene July 7 while 2024 scene August 4 and 2025 July 30: seasonal mismatch, not a matched day-of-year phenology comparison.",
            "Earth Search post-2022 BOA application and raster:bands offsets conflict. Retain both calibrations; do NOT interpret raw cross-year contrast as a harmonized LRT treatment effect.",
            "Scene-level cloud pct cannot substitute for SCL pixel cloud/shadow filtering; only SCL 4 vegetation or 5 not-vegetated are accepted.",
            "All samples have mixed land use and a 10-20m spatial footprint; no verified railway-verge management, train wind exposure or true plant community survey."
        ]
    }
    yearly={}
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",
        GDAL_HTTP_RETRY_DELAY="2",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="30000000",
        GDAL_CACHEMAX=256,
    ):
        for year,scene_ids in SCENES.items():
            try:
                items=[fetch_stac_scene(client,x) for x in scene_ids]
                if len({(x.get("properties") or {}).get("datetime","")[:10]
                        for x in items})!=1:
                    raise ValueError("Two LRT tiles not from the same acquisition date")
                pixel_values,quality=choose_year_pixel_items(year,points,items)
                yearly[year]=pixel_values
                calibration_conflict=any(metadata_ambiguous(x) for x in items)
                output["calendar"][str(year)]={
                    "scene_ids":list(scene_ids),
                    "sensing_date":items[0]["properties"]["datetime"][:10],
                    "n_fixed_pixel_sources_valid":len(pixel_values),
                    "n_source_pixel_locations":len(points),
                    "radiometry_metadata_conflict":calibration_conflict,
                    "source_qa":quality,
                    "within_year_annual_paired_summaries":{
                        label:paired_annual_ndvi(pairs,{year:pixel_values},year,label)
                        for label in ("provider_flag_branch","raster_offset_branch")
                    },
                }
                print(json.dumps({
                    "year":year,"sensing_date":output["calendar"][str(year)]["sensing_date"],
                    "valid_n":len(pixel_values),"total":len(points),
                    "calibration_conflict":calibration_conflict,
                    "within_year_paired_ndvi_provider":
                        output["calendar"][str(year)]["within_year_annual_paired_summaries"]["provider_flag_branch"]["near_minus_control_ndvi"],
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                output["errors"][str(year)]=f"{type(exc).__name__}: {str(exc)[:320]}"
                print(json.dumps({"year":year,"failure":output["errors"][str(year)]}),flush=True)
    for year in (2024,2025):
        for label in ("provider_flag_branch","raster_offset_branch"):
            entry=prepost_paired_change(pairs,yearly,year,label,2017)
            output["prepost_paired_differences"][f"{year}_{label}"]=entry
    # No healthy three-year comparison unless SAME spatial point pairs
    # have SCL-accepted source red/NIR values for all years.
    for label in ("provider_flag_branch","raster_offset_branch"):
        all_valid=set.intersection(*[
            set(yearly.get(y,{})) for y in SCENES
        ]) if all(y in yearly for y in SCENES) else set()
        output[f"n_fixed_pixels_valid_all_three_years_{label}"]=len(all_valid)
    output["status"]=(
        "observed_ndvi_pixel_values_but_not_calibrated_causal_effect"
        if yearly else "no_real_ndvi_pixel_values"
    )
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/satellite/utsunomiya_lrt_s2_real_pixel_ndvi.json")
    args=p.parse_args()
    report=run()
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    if not report["calendar"]:
        raise SystemExit("Zero actual Sentinel2 imagery pixel source success; no NDVI result.")
    print(json.dumps({
        "status":report["status"],
        "successful_years":list(report["calendar"]),
        "failures":report["errors"],
        "cross_year_descriptive_branches":report["prepost_paired_differences"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
