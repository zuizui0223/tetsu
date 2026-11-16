#!/usr/bin/env python3
"""Improve actual preconstruction S2 baseline without cherry-picking NDVI.

2017-07-22 (closest to Aug 1, scene clouds high) and 2017-07-07
(secondary) are the only verified 2017 July10-Aug25 L2A dates in
Earth Search's exhaustively paginated 2017 source query. Require BOTH
near/far pixels clear under SCL 4/5 in ONE acquisition. Choose the
nearest source date to Aug1 by SCL coverage only, before viewing NDVI.

Use source-validated reprocessed C1 L2A red/NIR with correct raster
offset on 2024-08-04 and 2025-07-30. The baseline uses old 2017
baseline 00.01 with zero BOA offset, because C1 was unavailable in 2017.

Even a change at fixed sample points is mixed urban land cover of
one non-independent tram route, not evidence that train wind, railway
construction, grass cutting or dispersal altered plant fitness.
"""
from __future__ import annotations
import argparse
from collections import Counter
import datetime as dt
import hashlib
import json
from pathlib import Path

import rasterio
import requests

from audit_n05_historical_railway import bounded_zip_download
from measure_utsunomiya_s2_ndvi_pixels import (
    SCENES,get_lrt_lines,build_fixed_pairs,choose_year_pixel_items,
    fetch_stac_scene
)
from verify_utsunomiya_s2_c1_pixels import (
    C1_SCENES,fetch_collection1_item,collection1_calibrated_pixels
)
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.sentinel_ndvi import paired_annual_ndvi,prepost_paired_change

BASELINE_DATES=(
    ("2017-07-22",(
        "S2B_54SVF_20170722_0_L2A",
        "S2B_54SUF_20170722_0_L2A",
    )),
    ("2017-07-07",SCENES[2017]),
)
EXPECTED_BASELINE_BOA_ZERO=True
CLOUD_SCL_ALLOWED=(4,5)


def baseline_pair_first_clear(pairs,pixel_by_date):
    """Choose closest-to-Aug1 source date with BOTH pair pixels valid.

    Unlike pixel-first selection, prevents near 7/22 vs far 7/07
    being mistaken for a same-date near/far vegetation difference.
    """
    source={}
    selections=Counter()
    pair_dates={}
    for near,far in pairs:
        for ds in (d for d,_ in BASELINE_DATES):
            records=pixel_by_date.get(ds,{})
            va=records.get(near,{}).get("provider_flag_branch")
            vb=records.get(far,{}).get("provider_flag_branch")
            if va is None or vb is None:
                continue
            source[near]={"harmonized":va}
            source[far]={"harmonized":vb}
            selections[ds]+=1
            pair_dates[(near,far)]=ds
            break
    return source,dict(selections),pair_dates


def run(n_anchors=9):
    source=bounded_zip_download()
    lines,_=get_lrt_lines(source)
    points,pairs,selection=build_fixed_pairs(lines,n_anchors=n_anchors)
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-C1-S2-true-same-date-preconstruction-LRT-pixel/1.3"})
    output={
        "run_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "MLIT_N05_sha256":hashlib.sha256(source).hexdigest(),
        "MLIT_N05_not_redistributed":True,
        "sampling_frame":selection,
        "posthoc_coverage_density_sensitivity_not_preregistered":n_anchors!=9,
        "source_2017_scenes":[{"date":d,"scene_ids":list(ids)}
                               for d,ids in BASELINE_DATES],
        "source_2017_scene_sort":"closest Julian day to Aug 1 first, then fixed July 7 fallback; never NDVI-selected",
        "source_2024_2025":"Reprocessed Sentinel2 C1 L2A, pixel scale=0.0001 offset=-0.1, crosschecked with legacy on 49 exact same-date pixels in run #42441974063",
        "2017_reprocessed_C1_unavailable_in_this_source":True,
        "quality":{},
        "results":{},
        "errors":{},
        "disclaimers":[
            "Baseline July 22 has scene-cloud cover >60% and only local SCL 4/5 valid pixels are retained.",
            "July 7 fallback is 25 days earlier than target Aug 1, July 22 is 10 days earlier; 2024 and 2025 post observations are Aug 4 and July 30.",
            "A near/far site pair is always observed on the SAME baseline date, but different pairs may have different baseline dates.",
            "All pairs belong to ONE urban tram line; no independent railway treatment replication or road-only comparator.",
            "10-20m S2 pixels mix sidewalk, roofs, roads, ballast, gardens and grasses; NDVI is not plant species, seed arrival, fitness or railway mowing.",
            "2017 and 2024/25 sensors/processing history differ; C1 reprocessed calibration is used post-2022 but 2017 legacy was not available in C1.",
            "Construction began 2018, 2023 passenger operations; 2018-22 cannot serve as an untreated baseline.",
            "No p-values, train wind adaptive growth or seed dispersal effect can be supported by the source alone."
        ]
    }
    per_date={}
    data={}
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="2",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        VSI_CACHE="TRUE",VSI_CACHE_SIZE="30000000",GDAL_CACHEMAX=256
    ):
        for ds,ids in BASELINE_DATES:
            try:
                items=[fetch_stac_scene(client,id) for id in ids]
                for item in items:
                    props=item.get("properties") or {}
                    if props.get("earthsearch:boa_offset_applied") is not False:
                        raise ValueError("2017 source has nonzero/unknown BOA flags")
                    for key in ("red","nir"):
                        rb=(item["assets"][key].get("raster:bands") or [{}])[0]
                        if abs(float(rb.get("offset") or 0))>1e-6:
                            raise ValueError("2017 baseline not zero-offset calibrated")
                pixels,quality=choose_year_pixel_items(2017,points,items)
                per_date[ds]=pixels
                output["quality"][ds]={
                    "valid_pixel_locations":len(pixels),
                    "same_date_valid_pairs":sum(
                        a in pixels and b in pixels for a,b in pairs
                    ),
                    "source_qa":quality
                }
                print(json.dumps({
                    "2017_date":ds,"valid_pixels":len(pixels),
                    "valid_pairs":output["quality"][ds]["same_date_valid_pairs"]
                }),flush=True)
            except Exception as exc:
                output["errors"][ds]=f"{type(exc).__name__}: {str(exc)[:260]}"
        base,counts,_=baseline_pair_first_clear(pairs,per_date)
        data[2017]=base
        output["baseline_2017_selected_same_date_pairs"]=counts
        output["baseline_2017_total_selected_pairs"]=sum(counts.values())
        output["2017_spectral_paired_summary"]=paired_annual_ndvi(
            pairs,{2017:base},2017,"harmonized"
        )
        for yr in (2024,2025):
            try:
                items=[fetch_collection1_item(client,id) for id in C1_SCENES[yr]]
                entry,q=collection1_calibrated_pixels(points,items)
                data[yr]=entry
                output["quality"][str(yr)]={
                    "valid_pixel_locations":len(entry),
                    "source_quality":q,
                }
                outcome=prepost_paired_change(pairs,data,yr,"harmonized",2017)
                output["results"][str(yr)]={
                    "source_annual_C1_paired_summary":paired_annual_ndvi(
                        pairs,{yr:entry},yr,"harmonized"
                    ),
                    "same_pixel_near_minus_control_NDVI_change_since_2017":
                        outcome,
                    "biological_effect_not_identified":True,
                }
                print(json.dumps({
                    "year":yr,
                    "C1_valid_pixels":len(entry),
                    "n_2017_post_valid_pairs":
                        outcome["paired_near_minus_control_delta_NDVI"]["n"],
                    "descriptive_delta":
                        outcome["paired_near_minus_control_delta_NDVI"]["median"]
                }),flush=True)
            except Exception as exc:
                output["errors"][str(yr)]=f"{type(exc).__name__}: {str(exc)[:260]}"
    output["status"]=(
        "measured_pixel_change_with_strict_multi_date_baseline_NOT_causal"
        if output["results"] and output["baseline_2017_total_selected_pairs"]
        else "no_comparable_real_pixel_source"
    )
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/satellite/utsunomiya_2017_clear_composite_C1_NDVI.json")
    p.add_argument("--anchors",type=int,default=9)
    a=p.parse_args()
    if not (8<=a.anchors<=80):
        p.error("--anchors must be between 8 and 80")
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    result=run(n_anchors=a.anchors)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "status":result["status"],"baseline":result["baseline_2017_selected_same_date_pairs"],
        "real_2024_2025":result["results"],"source_errors":result["errors"]
    },ensure_ascii=False),flush=True)
    if result["status"]=="no_comparable_real_pixel_source":
        raise SystemExit("No actual same-day/COG-calibrated near-far baseline and postsource")


if __name__=="__main__":
    main()
