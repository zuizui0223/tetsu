#!/usr/bin/env python3
"""Independent 2021 ESA WorldCover class audit for the SAME 10 N02/N13 rail sites.

No vegetation/NDVI outcome is used to choose the 10 GIS sites. Sample 2021
ESA WorldCover 10m 11-class landcover and its 3x3 neighborhood at fixed
35m rail and 1600m road-only candidate positions. This map distinguishes
tree canopy from mapped grassland/crops/built surfaces; a grassland map
pixel still does NOT prove railway-owned grass or at-grade exposure.

Use transient public 2021 ESA COG, N02-25 and N13-24; publish ONLY
aggregate class counts by independently mapped railway route. Never
republish source GIS geometries, sensitive sample locations or image tiles.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import zipfile

import rasterio
from rasterio.windows import Window
from pyproj import Transformer

from audit_japan_ordinary_rail_road_separation import (
    N02,MAX_ZIP_RAIL,REGIONS,analyze_region,download_source
)
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.worldcover_landcover import (
    CLASS_LABELS,DOI,REF_YEAR, SOURCE_BASE,SOURCE_ALT,
    tile_name,tile_url,normalize_class,summarize_neighborhood,
    summarize_route_source
)

SOURCE_TILE = "N36E138"  # all 3 fixed study boxes are 36..39N / 138..141E
EXPECTED_REGION_SOURCE_COUNTS={
    "sendai_miyagi":2,
    "niigata":4,
    "nagano":4
}
CLOUD_SCL_2025_RELATED_RUN="42448654924" # prior independent year-source result


def source_sample(dataset, x_utm54, y_utm54):
    """Read actual 2021 landcover class + 3x3 window at projected 2025 GIS point.

    Raster is in EPSG:4326, while source railway GIS points are in
    EPSG:32654. The source coordinate must be transformed exactly once.
    """
    if dataset.count!=1 or dataset.crs.to_epsg()!=4326:
        raise ValueError("ESA raster is not single-band WGS84")
    to_geo=Transformer.from_crs(32654,4326,always_xy=True)
    lon,lat=to_geo.transform(x_utm54,y_utm54)
    if tile_name(lon,lat)!=SOURCE_TILE:
        raise ValueError("Candidate outside official preselected WorldCover source tile")
    if not (
        dataset.bounds.left <= lon < dataset.bounds.right and
        dataset.bounds.bottom < lat <= dataset.bounds.top
    ):
        raise ValueError("ESA source data tile does not cover site")
    row,col=dataset.index(lon,lat)
    if not 1<=row<dataset.height-1 or not 1<=col<dataset.width-1:
        raise ValueError("Pixel too close to source tile edge for 3x3 audit")
    image=dataset.read(
        1,window=Window(col-1,row-1,3,3),
        boundless=True,fill_value=0
    )
    if image.shape!=(3,3):
        raise ValueError("ESA 3x3 pixel source read incomplete")
    main=normalize_class(int(image[1,1]))
    patch=summarize_neighborhood(image.flatten().tolist())
    return {"center_class":main,"window":patch}


def open_public_ESA_COG():
    errors=[]
    for root in (SOURCE_BASE,SOURCE_ALT):
        url=tile_url(SOURCE_TILE,root)
        try:
            src=rasterio.open(url)
            if (src.count!=1 or src.crs.to_epsg()!=4326 or
                src.width<10000 or src.height<10000 or
                src.nodata not in (0,None)):
                src.close()
                raise ValueError("Unexpected ESA WorldCover map schema")
            # Band class-domain and tile bounds will be confirmed at each
            # read; no optical SCL values are mistaken for landcover codes.
            return src,url
        except Exception as ex:
            errors.append(f"{type(ex).__name__}: {str(ex)[:160]}")
    raise RuntimeError(f"WorldCover 2021 public COG unavailable: {errors}")


def run():
    raw,source=download_source(N02,MAX_ZIP_RAIL)
    if not raw or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("Official N02 2025 source is not ZIP")
    result={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "kind":"verified 2021 land cover at road-separated ordinary rail sample opportunities, NOT botanical railway treatment effects",
        "N02_2025_source_sha256":source["source_sha256"],
        "N13_reference_year":2024,
        "prior_2025_Sentinel2_source_run":CLOUD_SCL_2025_RELATED_RUN,
        "WorldCover_product":"ESA WorldCover 2021 v200",
        "WorldCover_source_doi":DOI,
        "WorldCover_reference_year":REF_YEAR,
        "WorldCover_nominal_pixel_m":10,
        "WorldCover_source_tile":SOURCE_TILE,
        "WorldCover_EPSG":4326,
        "2021_map_class_codes":{str(k):v for k,v in CLASS_LABELS.items()},
        "reference_2021_to_2025_years_apart":4,
        "regions":{},"source_errors":{},
        "interpretation_limits":[
            "2021 WorldCover class is an algorithmic land-cover MAP, not a 2025/26 survey, actual train wind, grass cutting, herbaceous trait or plant taxon.",
            "Tree cover can obscure herbaceous plants below, and grassland class does not prove railway-owned vegetation, exposed trackbed, or at-grade railway operation.",
            "ESA WorldCover 2021 Sentinel1/2 and 2025 Sentinel-2 scene-classification share satellite inputs in a broad sense and are not biologically independent field measurements.",
            "N02 35m from line may be under tree canopy above tunnels, along bridges, or farmland and forests outside railway right-of-way.",
            "A 3x3 WorldCover window covers nominally 30x30m and may straddle multiple land uses; per-site class majority is descriptive not ground truth.",
            "No 2020-to-2021 landcover-change inference because WorldCover v100 and v200 algorithms differ.",
            "GIS road-separation eligibility was fixed BEFORE 2021 landcover or 2025 greenness was checked. Do not select on current vegetation for a future retrospective rail treatment-effect estimate.",
            "No ecological railway seed-spore transport, taxon richness, flower visibility, mowing regrowth, heritable adaptation or railway-caused cover change is measured.",
            "10 source samples are nested within five candidate route groups and three preselected regions; never use pair counts as independent experimental rail replicates."
        ]
    }
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",
        GDAL_HTTP_RETRY_DELAY="2",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="40000000",
        GDAL_CACHEMAX=256,
    ):
        source_map,actual_url=open_public_ESA_COG()
        result["WorldCover_verified_source_url"]=actual_url
        result["WorldCover_raster_metadata"]={
            "width":source_map.width,"height":source_map.height,
            "band_count":source_map.count,"nodata":source_map.nodata,
            "data_type":source_map.dtypes[0],
            "source_bounds_round_2dec":[round(x,2) for x in source_map.bounds]
        }
        try:
            for region in REGIONS:
                try:
                    report,sites=analyze_region(
                        region,raw,private_site_positions=True
                    )
                    if len(sites)!=EXPECTED_REGION_SOURCE_COUNTS[region]:
                        raise ValueError(
                            f"Source GIS sample cohort changed in {region}: "
                            f"{len(sites)} != {EXPECTED_REGION_SOURCE_COUNTS[region]}"
                        )
                    by_route={}
                    sampled=[]
                    for source_row in sites:
                        near=source_sample(source_map,*source_row["near_xy"])
                        control=source_sample(source_map,*source_row["road_only_xy"])
                        group=(source_row["route_name"],
                               source_row["operator_name"])
                        entry={
                            "near_class":near["center_class"],
                            "control_class":control["center_class"],
                            "near_window":near["window"],
                            "control_window":control["window"],
                        }
                        by_route.setdefault(group,[]).append(entry)
                        sampled.append(entry)
                    region_output={
                        "n_source_geographic_pairs":len(sites),
                        "n_source_routes":len(by_route),
                        "region_classes":summarize_route_source(sampled),
                        "route_classes":[
                            {"route_name":key[0],"operator_name":key[1],
                             **summarize_route_source(values)}
                            for key,values in sorted(by_route.items())
                        ],
                        "source_N13_2024":report["source_roads"]["source_url"],
                        "source_N13_sha256":report["source_roads"]["source_sha256"],
                    }
                    result["regions"][region]=region_output
                    print(json.dumps({
                        "region":region,"n_geo":len(sites),
                        "near_class_counts":
                            region_output["region_classes"]["near_2021_worldcover_class_counts"],
                        "road_class_counts":
                            region_output["region_classes"]["road_control_2021_worldcover_class_counts"],
                        "n_near_grassland_or_wetland":
                            region_output["region_classes"]["n_near_2021_worldcover_grassland_or_herbaceous_wetland"],
                        "n_near_tree_cover":
                            region_output["region_classes"]["n_near_2021_worldcover_tree_cover"],
                    },ensure_ascii=False),flush=True)
                except Exception as exc:
                    result["source_errors"][region]=f"{type(exc).__name__}: {str(exc)[:360]}"
        finally:
            source_map.close()
    totals=Counter()
    for entry in result["regions"].values():
        x=entry["region_classes"]
        totals["n_rail_road_geo_pairs"]+=x["n_rail_road_geometry_selected_pairs"]
        totals["near_tree_cover"]+=x["n_near_2021_worldcover_tree_cover"]
        totals["near_grassland_or_wetland"]+=x["n_near_2021_worldcover_grassland_or_herbaceous_wetland"]
        totals["near_and_control_ground_herbaceous"]+=x["n_near_and_control_potential_ground_herbaceous"]
        totals["near_2021_class_unclassified"]+=x["near_2021_worldcover_class_counts"].get("unknown",0)
    result["pooled_source_counts_NOT_independent_effect_units"]=dict(totals)
    result["source_status"]=(
        "2021_landcover_real_COG_verified_at_all_10_fixed_sites"
        if len(result["regions"])==len(REGIONS) and totals["n_rail_road_geo_pairs"]==10
        else "incomplete_regional_source_no_10_site_ecological_conclusion"
    )
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="results/audit/japan_worldcover_2021_ordinary_rail_landclass_v16.json")
    args=parser.parse_args()
    result=run()
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "status":result["source_status"],
        "pooled":result["pooled_source_counts_NOT_independent_effect_units"],
        "errors":result["source_errors"]
    },ensure_ascii=False),flush=True)
    if result["source_status"]!="2021_landcover_real_COG_verified_at_all_10_fixed_sites":
        raise SystemExit("ESA WorldCover source incomplete; do not infer 2021 grass or forest at all candidate locations")


if __name__=="__main__":
    main()
