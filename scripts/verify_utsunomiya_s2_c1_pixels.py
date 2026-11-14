#!/usr/bin/env python3
"""Real SAME-DATE Sentinel2 Collection-1 calibration validation on Utsunomiya.

C1 L2A is an independently reprocessed source with explicit -0.1
raster-band offset (DN -> scaled reflectance). Compare the same date,
tile and fixed pixel points against legacy L2A (whose STAC BOA flags
contradict stale per-band raster offset).

If physically agreeing, use C1-calibrated 2024 and 2025 values with
legacy 2017 baseline (2017 C1 unavailable in this source audit).
No N05 restricted geometry or original geo-coordinates are published.
This is NOT a causal rail greenness, seed, mowing or flowering effect.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import io
import json
import math
import statistics
import zipfile
from pathlib import Path
from urllib.parse import urlparse

import rasterio
import requests

from audit_n05_historical_railway import bounded_zip_download
from audit_utsunomiya_s2_c1 import COLLECTION
from measure_utsunomiya_s2_ndvi_pixels import (
    SCENES, STAC_ITEM_BASE, get_lrt_lines, build_fixed_pairs,
    fetch_stac_scene, read_scene_pixels, choose_year_pixel_items
)
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.sentinel_ndvi import (
    calculate_ndvi, summary, paired_annual_ndvi,prepost_paired_change,
    CLEAR_SCL_PRIMARY
)

C1_ITEM_ROOT="https://earth-search.aws.element84.com/v1/collections/sentinel-2-c1-l2a/items/"
C1_SCENES={
    2024:(
        "S2B_T54SVF_20240804T012656_L2A",
        "S2B_T54SUF_20240804T012656_L2A",
    ),
    2025:(
        "S2B_T54SVF_20250730T012655_L2A",
        "S2B_T54SUF_20250730T012655_L2A",
    )
}
MIN_CALIBRATION_COMPARABLE_PIXELS=4
MAX_SAME_SCENE_MEDIAN_ABS_NDVI_DIFF=.06


def fetch_collection1_item(client,id):
    r=client.get(C1_ITEM_ROOT+id,timeout=(20,75))
    r.raise_for_status()
    item=r.json()
    if item.get("id")!=id or item.get("collection")!=COLLECTION:
        raise ValueError("Unverified C1 STAC source item")
    for key in ("red","nir","scl"):
        a=(item.get("assets") or {}).get(key)
        if a is None:
            raise ValueError(f"C1 {id} missing {key}")
        part=urlparse(str(a.get("href") or ""))
        if part.scheme!="https" or not (part.hostname or "").endswith(".amazonaws.com") or not part.path.endswith(".tif"):
            raise ValueError("Unexpected C1 COG source URL")
    for key in ("red","nir"):
        spec=(item["assets"][key].get("raster:bands") or [{}])[0]
        if abs(float(spec.get("scale") or 0)-.0001)>1e-8 or abs(float(spec.get("offset") or 0)+.1)>1e-6:
            raise ValueError(f"Unexpected C1 raster calibration {id} {key}:{spec}")
    return item


def collection1_calibrated_pixels(points,items):
    entries={}
    qc={}
    for item in items:
        raw=read_scene_pixels(item,points)
        for point in points:
            if point in entries:continue
            red=raw["red"].get(point)
            nir=raw["nir"].get(point)
            scl=raw["scl"].get(point)
            if red is None or nir is None or scl is None:
                qc["out_of_tile_or_nodata"]=qc.get("out_of_tile_or_nodata",0)+1
                continue
            classification=int(scl)
            if classification not in CLEAR_SCL_PRIMARY:
                qc[f"scl_{classification}_excluded"]=qc.get(f"scl_{classification}_excluded",0)+1
                continue
            ndvi=calculate_ndvi(red,nir,subtract_dn=1000)
            if ndvi is None:
                qc["calibration_rejected"]=qc.get("calibration_rejected",0)+1
                continue
            entries[point]={
                "harmonized":ndvi,
                "scl_class":classification,
                "source_collection":COLLECTION
            }
            qc["source_pixels_valid"]=qc.get("source_pixels_valid",0)+1
    return entries,qc


def same_scene_calibration_differences(old,new):
    values=[]
    for name,pixel in new.items():
        reference=old.get(name,{}).get("provider_flag_branch")
        if reference is not None:
            values.append(pixel["harmonized"]-reference)
    details={
        "same_coordinate_same_date_comparable_valid_points":len(values),
        "collection1_minus_legacy_provider_ndvi":summary(values),
        "absolute_delta_ndvi":summary([abs(v) for v in values]),
    }
    median_abs=details["absolute_delta_ndvi"]["median"]
    details["calibration_cross_source_consistency_gate"]=(
        "pass" if len(values)>=MIN_CALIBRATION_COMPARABLE_PIXELS and
        median_abs is not None and median_abs<=MAX_SAME_SCENE_MEDIAN_ABS_NDVI_DIFF
        else "fail_insufficient_or_disagreeing_source"
    )
    return details


def run():
    raw_n05=bounded_zip_download()
    lines,_=get_lrt_lines(raw_n05)
    points,pairs,source_geometry=build_fixed_pairs(lines)
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-reprocessed-collection1-real-pixel-source/1.2"})
    result={
        "run_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "measurement":"real Sentinel-2 L2A Collection1-calibrated 2024/2025 vs legacy same acquisition and 2017 baseline",
        "source_collection1":COLLECTION,"source_legacy":"sentinel-2-l2a",
        "source_c1_no_2017_scenes_in_earlier_census":True,
        "source_2017_calibration":"legacy 2017 item baseline00.01, BOA offset=0 in source metadata",
        "source_N05_sha256":hashlib.sha256(raw_n05).hexdigest(),
        "licensed_geometry_not_published":True,
        "sample_frame":source_geometry,
        "source_same_day_checks":{},
        "matched_baseline_post_pixel_pairs":{},
        "source_failures":{},
        "important_limits":[
            "Near points are 140m from a single light-rail corridor, far points 1600m; NOT surveyed railway vegetation vs matched roads.",
            "One tram corridor, no independent railway replications or soil/land-cover/road/cutting exposure.",
            "2017 baseline July 7, 2024 August 4 and 2025 July 30: weather and season mismatch. Construction started 2018.",
            "2024/2025 C1 L2A uses raster:bands scale/offset; legacy data may have double-applied metadata offset.",
            "2017 Collection1 is absent in this first source screen; old baseline00.01 2017 is NOT an independently reprocessed baseline.",
            "All cloud/shadow/water/unclassified pixels discarded under SCL4/5 only.",
            "NDVI measured greenness, not seed dispersal, biodiversity, flower timing, spore-cone cover, mowing response or selection.",
            "Difference in NDVI is exploratory within one correlated corridor, NOT a causal effect of railway construction or operations.",
        ]
    }
    values={}
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="2",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        VSI_CACHE="TRUE",VSI_CACHE_SIZE="30000000",GDAL_CACHEMAX=256
    ):
        try:
            baseline_items=[fetch_stac_scene(client,x) for x in SCENES[2017]]
            pre_pixels,qc=choose_year_pixel_items(2017,points,baseline_items)
            values[2017]={
                k:{"harmonized":v["provider_flag_branch"]}
                for k,v in pre_pixels.items()
                if v["provider_flag_branch"] is not None
            }
            result["baseline2017_source"]={
                "legacy_scene_ids":list(SCENES[2017]),
                "n_baseline_clear_ndvi_pixels":len(values[2017]),
                "quality":qc
            }
        except Exception as ex:
            result["source_failures"]["2017"]=f"{type(ex).__name__}: {str(ex)[:300]}"
        for year,c1_scene_ids in C1_SCENES.items():
            try:
                c1_items=[fetch_collection1_item(client,x) for x in c1_scene_ids]
                old_items=[fetch_stac_scene(client,x) for x in SCENES[year]]
                c1dates={(item.get("properties") or {}).get("datetime","")[:10] for item in c1_items}
                olddates={(item.get("properties") or {}).get("datetime","")[:10] for item in old_items}
                if len(c1dates)!=1 or olddates!=c1dates:
                    raise ValueError(f"Legacy vs C1 acquisition date differs for {year}")
                c1_pixels,c1_qc=collection1_calibrated_pixels(points,c1_items)
                old_pixels,old_qc=choose_year_pixel_items(year,points,old_items)
                concordance=same_scene_calibration_differences(old_pixels,c1_pixels)
                values[year]=c1_pixels
                result["source_same_day_checks"][str(year)]={
                    "collection1_scene_ids":list(c1_scene_ids),
                    "legacy_scene_ids":list(SCENES[year]),
                    "acquisition_date":next(iter(c1dates)),
                    "collection1_source_pixel_count":len(c1_pixels),
                    "legacy_provider_flag_pixel_count":len(old_pixels),
                    "collection1_cloud_qa":c1_qc,
                    "legacy_cloud_qa":old_qc,
                    "source_calibration_concordance":concordance,
                    "collection1_within_year_paired_near_minus_far":
                        paired_annual_ndvi(pairs,{year:c1_pixels},year,"harmonized"),
                }
                print(json.dumps({
                    "year":year,"collection1_pixels":len(c1_pixels),
                    "legacy_pixels":len(old_pixels),
                    "crosssource":concordance,
                    "within_year_ndvi":
                        result["source_same_day_checks"][str(year)]["collection1_within_year_paired_near_minus_far"]["near_minus_control_ndvi"]
                },ensure_ascii=False),flush=True)
            except Exception as ex:
                result["source_failures"][str(year)]=f"{type(ex).__name__}: {str(ex)[:300]}"
    for year in C1_SCENES:
        result["matched_baseline_post_pixel_pairs"][str(year)]=(
            prepost_paired_change(pairs,values,year,"harmonized",2017)
        )
    result["source_quality_status"]=(
        "C1_reprocessed_calibrated_pixel_comparison_available"
        if result["source_same_day_checks"]
        else "no_actual_C1_pixel_values"
    )
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/satellite/utsunomiya_s2_c1_legacy_pixel_comparison.json")
    a=p.parse_args()
    result=run()
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "status":result["source_quality_status"],
        "crosssource":{y:x["source_calibration_concordance"]
                       for y,x in result["source_same_day_checks"].items()},
        "fixed_prepost":result["matched_baseline_post_pixel_pairs"],
        "errors":result["source_failures"]
    },ensure_ascii=False))
    if not result["source_same_day_checks"]:
        raise SystemExit("Collection 1 calibrated COG pixel access failed; no harmonized vegetation comparison")


if __name__=="__main__":
    main()
