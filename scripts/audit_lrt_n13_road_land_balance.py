#!/usr/bin/env python3
"""REAL MLIT N13 road proximity + 2017 preconstruction SCL/NDVI comparability.

This is a diagnostic *postconstruction road vs preconstruction spectral*
source audit for one 2023 Utsunomiya LRT corridor. Road centreline N13 is
dated September 2024; using it as 2017 baseline infrastructure would be
post-treatment bias. This audit DOES NOT fit a railway causal effect.
Original N05 noncommercial geometry and N13 road GIS are fetched
transiently, neither feature coordinates nor source files are published.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import statistics
import tempfile

import rasterio
import requests

from audit_n05_historical_railway import bounded_zip_download
from measure_utsunomiya_s2_ndvi_pixels import (
    get_lrt_lines,build_fixed_pairs,fetch_stac_scene,
    choose_year_pixel_items
)
from measure_lrt_2017_pixel_clear_composite import (
    BASELINE_DATES,baseline_pair_first_clear
)
from verify_utsunomiya_s2_c1_pixels import (
    C1_SCENES,fetch_collection1_item,collection1_calibrated_pixels
)

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_road_landbalance import (
    build_n13_road_index, pair_balance_flags, count_balance,blocked_summary,
    ROAD_DISTANCE_MAX_FOR_BALANCE_M,PRE_NDVI_MAX_DIFF
)


def pct_bin(v):
    if v<=25:return "0-25m"
    if v<=50:return "25-50m"
    if v<=100:return "50-100m"
    if v<=200:return "100-200m"
    return "gt200m"


def only_aggregates(points, pairs,road_at_point,pre_data,per_date):
    """Summarize N13 road exposure and pre-S2 variation without leaking geometry."""
    near=[road_at_point[a] for a,b in pairs]
    far=[road_at_point[b] for a,b in pairs]
    c=Counter()
    flags={}
    gaps=[]
    for near_id,far_id in pairs:
        # Baseline values have already been normalized to a mapping of
        # physical point ID -> scalar harmonized NDVI (or absent/None).
        # Do not access them as if they were nested dicts.
        one=pre_data.get(near_id)
        two=pre_data.get(far_id)
        scl_first=None;scl_second=None
        # Only one original acquisition date per pair. The source
        # selection is SCL-only, NOT picked to improve any contrast sign.
        for ds,_ in BASELINE_DATES:
            row=per_date.get(ds,{})
            if (row.get(near_id,{}).get("provider_flag_branch") is not None and
                row.get(far_id,{}).get("provider_flag_branch") is not None):
                scl_first=row[near_id].get("scl_class")
                scl_second=row[far_id].get("scl_class")
                break
        flags[(near_id,far_id)] = pair_balance_flags(
            road_at_point[near_id], road_at_point[far_id],
            scl_first, scl_second, one,two
        )
        if one is not None and two is not None:
            gaps.append(abs(one-two))
    for group,values in (("rail_near_offset_140m",near),
                         ("1p6km_far_control",far)):
        for x in values:
            c[(group,pct_bin(x))]+=1
    return {
        "road_distance_median_near":round(statistics.median(near),2),
        "road_distance_median_far":round(statistics.median(far),2),
        "road_distance_bins":{
            name:{bin:c[(name,bin)] for bin in
                  ("0-25m","25-50m","50-100m","100-200m","gt200m")}
            for name in ("rail_near_offset_140m","1p6km_far_control")
        },
        "2017_abs_ndvi_between_exposed_and_far_pairs_median":(
            round(statistics.median(gaps),5) if gaps else None
        ),
        "2017_abs_ndvi_between_exposed_and_far_pairs_n":len(gaps),
        "predefined_balance_gates":count_balance(flags)
    },flags


def run(n_anchors=48):
    raw_n05=bounded_zip_download()
    lines,n_lines=get_lrt_lines(raw_n05)
    pts,pairs,spatial=build_fixed_pairs(lines,n_anchors=n_anchors)
    with tempfile.TemporaryDirectory() as tmp:
        nearest_road,road_meta=build_n13_road_index(download_dir=Path(tmp))
        road_at_point={key:nearest_road(coordinate)
                       for key,coordinate in pts.items()}

    client=requests.Session()
    client.headers.update({
        "User-Agent":"tetsu-LRT-N13-road-baseline-landcover-control-v1.3"
    })
    per_date={}
    yearly={}
    errors={}
    source_checks={}
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="2",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        VSI_CACHE="TRUE",VSI_CACHE_SIZE="30000000",GDAL_CACHEMAX=256,
    ):
        for date,ids in BASELINE_DATES:
            try:
                scenes=[fetch_stac_scene(client,x) for x in ids]
                p,q=choose_year_pixel_items(2017,pts,scenes)
                per_date[date]=p
                source_checks[date]={
                    "n_scl_clear_pixels":len(p),
                    "n_near_far_pairs_both_scl_clear":sum(
                        a in p and b in p for a,b in pairs
                    ),
                    "QA":q
                }
            except Exception as ex:
                errors[date]=f"{type(ex).__name__}: {str(ex)[:240]}"

        base,chosen_dates,_=baseline_pair_first_clear(pairs,per_date)
        yearly[2017]={
            x:y["harmonized"] for x,y in base.items()
        }
        source_checks["2017_chosen_baseline_date_counts"]=chosen_dates
        for year in (2024,2025):
            try:
                scenes=[fetch_collection1_item(client,x)
                        for x in C1_SCENES[year]]
                p,q=collection1_calibrated_pixels(pts,scenes)
                yearly[year]={x:v["harmonized"] for x,v in p.items()}
                source_checks[str(year)]={
                    "collection1_calibrated_real_pixel_count":len(p),
                    "QA":q
                }
            except Exception as ex:
                errors[str(year)]=f"{type(ex).__name__}: {str(ex)[:240]}"
    balance,flags=only_aggregates(
        pts,pairs,road_at_point,yearly[2017],per_date
    )
    results={}
    for year in (2024,2025):
        if year not in yearly:
            continue
        results[str(year)]={
            gate:blocked_summary(pairs,yearly,flags,year,gate=gate)
            for gate in (
                "all",
                "near_road_le100_and_far_road_le100",
                "road_abs_diff_le50_both_le100",
                "2017_same_scl_class",
                "2017_both_vegetation_scl4",
                "2017_abs_NDVI_gap_le015",
                "road50_and_2017_same_class",
                "road50_sameclass_preNDVI_similar",
            )
        }
        print(json.dumps({
            "year":year,
            "source_availability":source_checks[str(year)],
            "all_valid":results[str(year)]["all"]["paired_complete"],
            "road50_and_scl_matched":
                results[str(year)]["road50_and_2017_same_class"]["paired_complete"],
            "strict_road_landbaseline_preNDVI":
                results[str(year)]["road50_sameclass_preNDVI_similar"]["paired_complete"]
        },ensure_ascii=False),flush=True)
    report={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "N05_source_sha256":hashlib.sha256(raw_n05).hexdigest(),
        "N05_license":"noncommercial-only; source geometry not copied to result",
        "N13_road_provenance":road_meta,
        "N13_official_year":2024,
        "N13_is_post_2018_rail_construction_and_NOT_valid_2017_control":True,
        "sampling_frame":spatial,
        "fixed_road_distance_thresholds_m":{
            "both_road_within":ROAD_DISTANCE_MAX_FOR_BALANCE_M,
            "road_absolute_difference":50,
            "sensitivity_road_absolute_difference":100
        },
        "fixed_2017_spectral_QC_max_abs_ndvi_difference":PRE_NDVI_MAX_DIFF,
        "2017_preconstruction_baseline_source_dates":[x[0] for x in BASELINE_DATES],
        "2024_2025_imagery":"reprocessed C1 2017 baseline legacy zero-offset, actual same fixed ground coordinates",
        "source_QA":source_checks,
        "balance":balance,
        "observed_results_NOT_causal":results,
        "source_errors":errors,
        "interpretation_limits":[
            "2024 N13 road distance is a postconstruction variable; must not claim to have balanced 2017 road exposure.",
            "No actual verified vegetation/right-of-way polygons; 140m near rail may be road, parking or garden, not trackside plants.",
            "Historical SCL vegetation/bare from two different 2017 dates, cloud/sensor effects and early July phenology cause selection bias.",
            "2017 reflectance is legacy baseline00.01; post 2024/25 reprocessed C1 with proven correct offset. Slight seasonal mismatch remains.",
            "Matching on baseline NDVI is a posthoc QC restriction which can induce regression-to-mean and does not identify causality.",
            "Same LRT corridor source, spatially correlated anchor pairs, NOT independent transport/operation experimental units.",
            "Road N13 centreline has no operational traffic and no guarantee of road existence in 2017; post-road differences may itself be an effect of railway construction.",
            "No weed species, seeded/sprouted propagules, mowing, herbicide, train wind, or heritable adaptation measured by 10m S2 NDVI.",
            "If complete road/baseline matched pairs fail replication gates, report no qualified ecological contrast instead of effect direction.",
        ],
    }
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--anchors",type=int,default=48)
    p.add_argument("--output",default="results/satellite/lrt_n13_2024_road_baseline_2017_balance_v13.json")
    args=p.parse_args()
    if args.anchors!=48:
        p.error("v1.3 fixed sensitivity uses exactly 48 anchors from v1.2")
    result=run(args.anchors)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "source_road_parts":result["N13_road_provenance"]["n_line_parts_in_spatial_index"],
        "road_pre_SCL_pair_gates":result["balance"]["predefined_balance_gates"],
        "yearly_results":{
            y:{
                "all":v["all"]["paired_complete"],
                "road_SCL":v["road50_and_2017_same_class"]["paired_complete"],
                "road_SCL_preNDVI":
                    v["road50_sameclass_preNDVI_similar"]["paired_complete"]
            }
            for y,v in result["observed_results_NOT_causal"].items()
        },
        "errors":result["source_errors"],
    },ensure_ascii=False,indent=2))
    if not result["observed_results_NOT_causal"]:
        raise SystemExit("No actual N13+S2 source overlap, cannot infer land-cover effects")


if __name__=="__main__":
    main()
