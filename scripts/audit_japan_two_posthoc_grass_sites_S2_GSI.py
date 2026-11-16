#!/usr/bin/env python3
"""Posthoc 2-site 2021 ESA grassland -> 2025 S2 SCL + GSI image availability QA.

In the explicitly posthoc 370-normal WorldCover audit, only 2 samples
had a nearby 2024 road centreline >=50m away, an eligible 1.6km road-only
control, AND WorldCover 2021 class30 grassland at the nominal 35m rail
point. Selection used WorldCover 2021 landcover, NOT 2025 NDVI.
Independent 2025 L2A imagery and recent GSI orthophoto *coverage* are
then checked at exactly the same physical rail/road source site coordinates.

No attempt to identify actual grasses, at-grade tracks, tree canopy,
railway property, cutting/herbicide/wind or seed arrival is justified
from this posthoc source cross-check. No exact GIS positions/images
are copied to GitHub or emitted to public logs/results.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import rasterio
import requests
from PIL import Image

from audit_japan_ordinary_rail_road_separation import (
    N02,MAX_ZIP_RAIL,download_source,analyze_region
)
from audit_japan_ordinary_rail_worldcover_2021 import (
    open_public_ESA_COG,source_sample
)
from audit_japan_ordinary_rail_2025_s2_vegetation import (
    sample_region_2025
)
from pyproj import Transformer

REGION_CANDIDATES = {
    "niigata": ("信越線","東日本旅客鉄道"),
    "nagano": ("小海線","東日本旅客鉄道"),
}
ALLOWED_NEAR_CLASS = {30,90}
GIS_ROAD_CUTOFF_M = 50.
GIS_FAR_ROAD_MAX_M = 50.
GIS_FAR_RAIL_MIN_M = 1000.
GIS_NEAR_STATION_MIN_M = 300.
GSI_SEAMLESS = "https://cyberjapandata.gsi.go.jp/xyz/seamlessphoto/{z}/{x}/{y}.jpg"
GSI_ATTRIBUTION = "国土地理院 地理院地図 シームレス空中写真"
GSI_ZOOM = (17,18,19)


def geo_tile(lon,lat,z):
    if not (135<=lon<=142 and 34<=lat<=40 and z in GSI_ZOOM):
        raise ValueError("Unexpected Japanese regional photography sample")
    scale=2**z
    x=int((lon+180)/360*scale)
    r=math.radians(lat)
    y=int((1-math.asinh(math.tan(r))/math.pi)/2*scale)
    return z,x,y


def gsi_source_tile_check(client,lon,lat):
    """Transient JPEG byte content, never source image/coordinates exported."""
    results={}
    for z in GSI_ZOOM:
        zz,x,y=geo_tile(lon,lat,z)
        url=GSI_SEAMLESS.format(z=zz,x=x,y=y)
        try:
            r=client.get(url,timeout=(12,40))
            if r.status_code==404:
                results[str(z)]="no_tile_404"
                continue
            r.raise_for_status()
            if len(r.content)>2_000_000:
                raise ValueError("GSI image tile too large")
            with Image.open(io.BytesIO(r.content)) as im:
                if im.format!="JPEG" or im.size not in ((256,256),(512,512)):
                    raise ValueError("Unexpected GSI JPEG schema")
                rgb=im.convert("RGB")
                width,height=rgb.size
                sample=[
                    rgb.getpixel((
                        int((i+.5)*width/16),
                        int((j+.5)*height/16)
                    ))
                    for i in range(16) for j in range(16)
                ]
                if len(set(sample))<=3:
                    results[str(z)]="solid_or_placeholder"
                else:
                    results[str(z)]="nonblank_current_aerial_source"
        except (requests.RequestException,ValueError,OSError) as exc:
            results[str(z)]=f"source_error_{type(exc).__name__}"
    return results


def eligible_posthoc_worldcover_grass(site,raster):
    if not (
        site["near_road_m"]>=GIS_ROAD_CUTOFF_M and
        site["far_road_m"]<=GIS_FAR_ROAD_MAX_M and
        site["far_rail_m"]>=GIS_FAR_RAIL_MIN_M and
        site["near_station_m"]>=GIS_NEAR_STATION_MIN_M
    ):
        return False
    near=source_sample(raster,*site["near_xy"])
    return near["center_class"] in ALLOWED_NEAR_CLASS


def run():
    raw,n02_source=download_source(N02,MAX_ZIP_RAIL)
    if not raw or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("Failed official N02 passenger rail source")
    result={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "scope":"posthoc two-site ground-vegetation SOURCE validation, not independent discovery or railway ecological effect",
        "source_N02_SHA256":n02_source["source_sha256"],
        "WorldCover_class_year":2021,
        "S2_summer_class_year":2025,
        "GSI_source_year":"current seamless mosaic; actual photography acquisition year unverified",
        "2021_ESA_grassland_class_code":30,
        "2024_road_distance_threshold_m":GIS_ROAD_CUTOFF_M,
        "n_expected_posthoc_grassland_candidates":2,
        "site_reviews":{},
        "errors":{},
        "critical_limits":[
            "Two samples were chosen after WorldCover 2021 results and are posthoc, not a new confirmatory population or replicated rail experiments.",
            "2021 WorldCover grassland class30 is not proof of grasses on railway-owned ballast/embankment, and 2025 SCL4 vegetation could be tree cover or farm growth.",
            "Current GSI seamless aerial tiles can be from multiple acquisition dates; availability does NOT verify 2017/2021/2025 roads or a railway at grade.",
            "Rail and road source positions are GIS line-normal offsets, not actual safe public survey access points; sampling requires railway authorization.",
            "No raw GSI photo/image, N02 or N13 GIS source coordinates, station IDs, 1:1 pixel screenshots or copyrighted images exported.",
            "No seed/spore transport, mowing resilience, phenology, reproductive ROI, wind exposure, adaptive genotype, richness or causal railway cover effect."
        ]
    }
    client=requests.Session()
    client.headers.update({
        "User-Agent":"tetsu-ESA-2021-grassland-two-site-C1-GSI-source/1.6"
    })
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="2",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        VSI_CACHE="TRUE",VSI_CACHE_SIZE="40000000",GDAL_CACHEMAX=256,
    ):
        raster,source_uri=open_public_ESA_COG()
        result["ESA_source_uri"]=source_uri
        try:
            for region,route in REGION_CANDIDATES.items():
                try:
                    _,all_sites=analyze_region(
                        region,raw,private_site_positions=True,
                        private_candidate_scope="all_fixed_geo"
                    )
                    matches=[
                        site for site in all_sites
                        if (site["route_name"],site["operator_name"])==route
                        and eligible_posthoc_worldcover_grass(site,raster)
                    ]
                    if len(matches)!=1:
                        raise ValueError(
                            f"Original 2021 posthoc single-site candidate changed for {region}/{route[0]}: {len(matches)}"
                        )
                    near,far=matches[0]["near_xy"],matches[0]["road_only_xy"]
                    control=source_sample(raster,*far)
                    satellite=sample_region_2025(region,matches,client)
                    to_geo=Transformer.from_crs(
                        32654,4326,always_xy=True
                    )
                    lon,lat=to_geo.transform(*near)
                    gsi_near=gsi_source_tile_check(client,lon,lat)
                    lon,lat=to_geo.transform(*far)
                    gsi_far=gsi_source_tile_check(client,lon,lat)
                    result["site_reviews"][region]={
                        "N02_route":route[0],
                        "N02_operator":route[1],
                        "2021_ESA_near_landcover":"grassland",
                        "2021_ESA_road_control_class_code":control["center_class"],
                        "2021_ESA_road_control_plur_class_code":
                            control["window"]["plurality_class"],
                        "2025_S2_satellite":satellite,
                        "GSI_latest_seamless_aerial_near_zoom_coverage":gsi_near,
                        "GSI_latest_seamless_aerial_control_zoom_coverage":gsi_far,
                        "image_human_ground_validation_completed":False,
                        "actual_railway_verge_grass_confirmed":False,
                    }
                    print(json.dumps({
                        "region":region,"N02_route":route[0],
                        "near2021_WorldCover":"grassland",
                        "2025_near_SCL4_on_date":
                            satellite["sum_near_scl4_one_date"],
                        "2025_near_far_both_SCL4_same_date":
                            satellite["sum_both_near_and_control_vegetated_same_date"],
                        "2025_same_day_conditional_NDVI":
                            satellite["2025_conditional_same_day_near_vs_road_vegetated_NDVI"],
                        "GSI_aerial_source_availability":gsi_near
                    },ensure_ascii=False),flush=True)
                except Exception as exc:
                    result["errors"][region]=f"{type(exc).__name__}: {str(exc)[:320]}"
        finally:
            raster.close()
    result["status"]=(
        "two_posthoc_2021_grass_sites_2025_SCL_GSI_source_checked"
        if len(result["site_reviews"])==len(REGION_CANDIDATES)
        else "incomplete_posthoc_grass_site_source_check"
    )
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/japan_two_posthoc_grass_site_S2_GSI_v16.json")
    args=p.parse_args()
    result=run()
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "areas_checked":list(result["site_reviews"]),
        "errors":result["errors"]
    },ensure_ascii=False))
    if result["status"]!="two_posthoc_2021_grass_sites_2025_SCL_GSI_source_checked":
        raise SystemExit("Cannot claim complete 2-site cross-era landcover source screening")


if __name__=="__main__":
    main()
