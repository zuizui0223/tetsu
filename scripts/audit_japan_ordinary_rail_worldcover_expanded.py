#!/usr/bin/env python3
"""POSTHOC source-opportunity expansion to all 370 fixed rail/road transects.

The initial 10 strict GIS candidates proved to be 9 ESA 2021 tree-cover
pixels plus 1 permanent-water pixel, NOT mapped grassland. Rather than
claim these are railside grass patches, audit the ORIGINAL fixed 370
N02/N13 normal transects at 35m with identical WorldCover 2021 v200.
This is exploratory site discovery following results; never an
independent confirmatory test or biological railway treatment estimate.

Nothing in the source report contains the licensed GIS coordinates,
precise train-adjacent locations, or raw ESA 2021 image data.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import zipfile

import rasterio

from audit_japan_ordinary_rail_road_separation import (
    N02,MAX_ZIP_RAIL,REGIONS,analyze_region,download_source
)
from audit_japan_ordinary_rail_worldcover_2021 import (
    open_public_ESA_COG,source_sample
)
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.worldcover_landcover import (
    CLASS_LABELS,DOI,POTENTIAL_GROUND_HERBS,
    summarize_route_source
)

EXPECTED_ALL_FIXED_GEO={
    "sendai_miyagi":81,"niigata":154,"nagano":135
}
EXPECTED_ORIGINAL_STRICT={
    "sendai_miyagi":2,"niigata":4,"nagano":4
}
ROAD_NEAR_THRESHOLDS=(25,50,75,100)
ROAD_CONTROL_MAX_M=50.
FAR_RAIL_MIN_M=1000.
NEAR_STATION_MIN_M=300.


def source_comparator_eligibility(site,near_min_road_m):
    if near_min_road_m not in ROAD_NEAR_THRESHOLDS:
        raise ValueError("Unapproved exploratory road distance cutoff")
    return (
        site["near_road_m"] >= near_min_road_m
        and site["far_road_m"] <= ROAD_CONTROL_MAX_M
        and site["far_rail_m"] >= FAR_RAIL_MIN_M
        and site["near_station_m"] >= NEAR_STATION_MIN_M
    )


def audit_region(region,sites,raster):
    if len(sites)!=EXPECTED_ALL_FIXED_GEO[region]:
        raise ValueError("Original fixed normal-transect sampling changed")
    rows=[]
    for site in sites:
        near=source_sample(raster,*site["near_xy"])
        control=source_sample(raster,*site["road_only_xy"])
        class_near=near["center_class"]
        class_far=control["center_class"]
        rows.append({
            "route_name":site["route_name"],
            "operator_name":site["operator_name"],
            "near_class":class_near,
            "control_class":class_far,
            "near_window":near["window"],
            "control_window":control["window"],
            "road_flags":{
                str(thresh):source_comparator_eligibility(site,thresh)
                for thresh in ROAD_NEAR_THRESHOLDS
            },
            "strict_original_GIS_source":site["source_gates"][
                "strict_near_road75_far_road50_far_rail1000_station300"
            ],
        })
    if sum(x["strict_original_GIS_source"] for x in rows)!=EXPECTED_ORIGINAL_STRICT[region]:
        raise ValueError("Original strict 10-site subset failed reproduction")
    if any(r["road_flags"]["75"]!=r["strict_original_GIS_source"] for r in rows):
        raise ValueError("Original GIS strict road75 gate has changed")

    by_cutoff={}
    all_routes=sorted({(x["route_name"],x["operator_name"]) for x in rows})
    for threshold in ROAD_NEAR_THRESHOLDS:
        eligible=[x for x in rows if x["road_flags"][str(threshold)]]
        grouped=defaultdict(list)
        for sample in eligible:
            grouped[(sample["route_name"],sample["operator_name"])].append(sample)
        candidates=[]
        for route in all_routes:
            entries=grouped[route]
            if not entries:continue
            source=summarize_route_source(entries)
            candidates.append({
                "route_name":route[0],
                "operator_name":route[1],
                "n_source_road_eligible":len(entries),
                "near_WorldCover_2021_lc_class_counts":
                    source["near_2021_worldcover_class_counts"],
                "n_near_grassland_or_herbaceous_wetland":
                    source["n_near_2021_worldcover_grassland_or_herbaceous_wetland"],
                "n_near_tree_cover":
                    source["n_near_2021_worldcover_tree_cover"],
                "n_both_road_and_near_WC_grassland":
                    source["n_near_and_control_potential_ground_herbaceous"],
            })
        candidates.sort(key=lambda row:(
            -row["n_near_grassland_or_herbaceous_wetland"],
            -row["n_source_road_eligible"],
            row["route_name"]
        ))
        classified=Counter(
            CLASS_LABELS.get(x["near_class"],"unknown") for x in eligible
        )
        by_cutoff[str(threshold)]={
            "n_2024_road_control_eligible_pairs":len(eligible),
            "near_2021_WorldCover_class_counts":dict(classified),
            "n_near_grassland_or_herbaceous_wetland":
                sum(x["near_class"] in POTENTIAL_GROUND_HERBS for x in eligible),
            "n_non_tree_near_not_crops_not_built_not_water":
                sum(x["near_class"] in (20,30,60,90) for x in eligible),
            "n_current_2024_road_source_route_groups":len(grouped),
            "route_source_coverage_ranked":
                candidates[:30]
        }

    all_class=summarize_route_source(rows)
    return {
        "region":region,
        "n_original_GIS_fixed_source_pairs":len(rows),
        "all_fixed_source_2021_WorldCover_class_counts":
            all_class["near_2021_worldcover_class_counts"],
        "all_fixed_source_tree_cover_n":
            all_class["n_near_2021_worldcover_tree_cover"],
        "all_fixed_source_grassland_or_herb_wetland_n":
            all_class["n_near_2021_worldcover_grassland_or_herbaceous_wetland"],
        "n_strict_original_road75_source_pairs":
            sum(x["strict_original_GIS_source"] for x in rows),
        "2024_road_distance_threshold_sensitivities":
            by_cutoff,
        "posthoc_discovery_NOT_confirmatory":True,
        "source_location_classifications_not_published":True
    }


def run():
    raw,rail_source=download_source(N02,MAX_ZIP_RAIL)
    if not raw or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("MLIT N02 2025 source missing")
    output={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "scientific_scope":"Exploratory potential rail-grass site discovery across preexisting fixed 370 railway/road geometry samples, NOT ecological effect",
        "selection_warning":"Expanded AFTER finding 9 tree and 1 water in original strict 10; report as posthoc reconnaissance, never independent positive treatment validation.",
        "N02_2025_SHA256":rail_source["source_sha256"],
        "N13_2024_road_centerline_not_historical":True,
        "ESA_WorldCover_2021_reference_doi":DOI,
        "ESA_2021_classes":{str(k):v for k,v in CLASS_LABELS.items()},
        "WorldCover_tile":"N36E138",
        "source_comparison_controls":{
            "near_rail_offset_m":35,
            "road_control_far_offset_m":1600,
            "near_min_distance_to_2024_roads_m":list(ROAD_NEAR_THRESHOLDS),
            "far_max_distance_to_2024_road_m":ROAD_CONTROL_MAX_M,
            "far_min_distance_to_any_passenger_rail_m":FAR_RAIL_MIN_M,
            "near_min_distance_to_passenger_station_m":NEAR_STATION_MIN_M,
        },
        "regions":{},"source_errors":{},
        "important_limits":[
            "The 2021 class at near35m does not establish an at-grade track or a vegetation patch on rail property; a pixel can be forest over a tunnel, floodplain water under bridge, farm or garden.",
            "WorldCover class30 grassland and class90 herbaceous wetlands are LANDCOVER MAP categories, not independently collected plant flora. Class10 tree cover can hide herbaceous understory.",
            "Lowering near-road threshold from the original 75m to 50/25m is explicitly POSTHOC GIS-ecological reconnaissance, not preregistered or an effect estimator.",
            "A control within 50m of the road is not necessarily a comparable road-grass verge; actual paved width and agricultural/forest habitat vary.",
            "No original polygons, individual site coordinates, railway properties or ESA raster imagery redistributed; only aggregated class tallies and public route names.",
            "Each fixed 2021 WorldCover map pixel uses Sentinel 1/2, so its class is not wholly sensor-independent from Sentinel2 2025 SCL.",
            "A route with several nearby eligible transects is not several independent railway treatments.",
            "No dispersal/seed/spore/plant species, wind gust, flowering, cutting, herbicide or genetically adaptive response can be inferred."
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
        raster,uri=open_public_ESA_COG()
        output["ESA_verified_source_url"]=uri
        try:
            for region in REGIONS:
                try:
                    meta,sites=analyze_region(
                        region,raw,private_site_positions=True,
                        private_candidate_scope="all_fixed_geo"
                    )
                    stats=audit_region(region,sites,raster)
                    output["regions"][region]=stats
                    print(json.dumps({
                        "region":region,
                        "all":stats["n_original_GIS_fixed_source_pairs"],
                        "all_near_class":stats["all_fixed_source_2021_WorldCover_class_counts"],
                        "road25":{
                            k:v for k,v in stats["2024_road_distance_threshold_sensitivities"]["25"].items()
                            if k!="route_source_coverage_ranked"
                        },
                        "road50":{
                            k:v for k,v in stats["2024_road_distance_threshold_sensitivities"]["50"].items()
                            if k!="route_source_coverage_ranked"
                        },
                        "road75":{
                            k:v for k,v in stats["2024_road_distance_threshold_sensitivities"]["75"].items()
                            if k!="route_source_coverage_ranked"
                        }
                    },ensure_ascii=False),flush=True)
                except Exception as exc:
                    output["source_errors"][region]=f"{type(exc).__name__}: {str(exc)[:320]}"
        finally:
            raster.close()
    totals=Counter()
    for region,d in output["regions"].items():
        totals["fixed_normal_pairs"]+=d["n_original_GIS_fixed_source_pairs"]
        totals["strict_original_road75_pairs"]+=d["n_strict_original_road75_source_pairs"]
        totals["all_fixed_grassland_or_wetland"]+=d["all_fixed_source_grassland_or_herb_wetland_n"]
        for thr in ROAD_NEAR_THRESHOLDS:
            subset=d["2024_road_distance_threshold_sensitivities"][str(thr)]
            totals[f"2024_road{thr}_eligible"]+=subset["n_2024_road_control_eligible_pairs"]
            totals[f"2024_road{thr}_grassland_or_herbaceous_wetland"]+=subset["n_near_grassland_or_herbaceous_wetland"]
            totals[f"2024_road{thr}_low_ground_candidate"]+=subset["n_non_tree_near_not_crops_not_built_not_water"]
    output["pooled_source_coverage_NOT_railway_effect"]=dict(totals)
    output["source_status"]=(
        "verified_all_370_posthoc_source_sites"
        if len(output["regions"])==3 and
        totals["fixed_normal_pairs"]==370 and
        totals["strict_original_road75_pairs"]==10
        else "incomplete_source_no_biological_comparison"
    )
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/japan_370_rail_normal_worldcover_2021_posthoc_v16.json")
    args=p.parse_args()
    result=run()
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "source_status":result["source_status"],
        "total":result["pooled_source_coverage_NOT_railway_effect"],
        "errors":result["source_errors"]
    },ensure_ascii=False),flush=True)
    if result["source_status"]!="verified_all_370_posthoc_source_sites":
        raise SystemExit("Incomplete prior-source 370 site universe; do not make site-ranking inference")


if __name__=="__main__":
    main()
