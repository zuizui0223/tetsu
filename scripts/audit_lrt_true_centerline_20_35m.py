#!/usr/bin/env python3
"""Real 20/35m official railway centreline satellite/road suitability audit.

N05 route geometry is noncommercial and fetched transiently, N13 is
2024 road geometry (AFTER 2018 tram construction), S2 L2A 2017 baseline
and C1 reprocessed 2024/25. 10m NDVI/20m SCL pixel at 20m from a
tram centreline is not independently verified railway habitat.
No ecological causal effect, seed migration, mowing or wind inference.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
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
from tetsu.rail_road_landbalance import build_n13_road_index
from tetsu.rail_verge_pixel_eligibility import (
    near_points_from_140m,analyze_pairs,OFFSETS
)

N_ANCHORS=48


def audit():
    n05=bounded_zip_download()
    lines,n_line_rows=get_lrt_lines(n05)
    original,original_pairs,spatial=build_fixed_pairs(lines,n_anchors=N_ANCHORS)
    points={}
    pairs_by_offset={}
    for offset in OFFSETS:
        mapped,pairs=near_points_from_140m(original,original_pairs,offset)
        points.update(mapped)
        pairs_by_offset[int(offset)]=pairs
    if len(points)!=len(original_pairs)*3:
        raise ValueError("unexpected or duplicated 20m/35m point coordinates")
    result={
        "source_N05_SHA256":hashlib.sha256(n05).hexdigest(),
        "source_N05_noncommercial_geometry_reused_transiently":True,
        "source_N05_LRT_qualified_line_rows":n_line_rows,
        "source_N13_2024_is_post_construction":True,
        "sampling_frame":spatial,
        "sampling_design":"48 even anchors, both available sides, exact 20m and 35m perpendicular offsets, 1600m old far sites; one correlated urban tram corridor",
        "n_geographic_fixed_side_pairs":len(original_pairs),
        "n_physical_sample_points":len(points),
        "2017_preconstruction_acquisitions":[d for d,_ in BASELINE_DATES],
        "post_C1_acquisitions":{"2024":"2024-08-04","2025":"2025-07-30"},
        "source_qa":{},"errors":{},"offsets":{},
        "scientific_caveats":[
            "The measured 20m/35m points are near an N05 track centreline, but DO NOT establish an actual vegetated railway right-of-way: 10m red/NIR and 20m SCL pixels mix paved track, asphalt, buildings, gardens and trees.",
            "N13 road distance is measured in 2024, six years after rail construction commenced in 2018; 2017 road existence is NOT identified by this source.",
            "Road-free >=25m is only centreline distance and is not proof of absence of road surface; wide roads/rights-of-way remain possible.",
            "Baseline source days Jul7/Jul22 differ from Aug4 2024 and Jul30 2025 and are nonrandomly missing due to cloud.",
            "All spatial points are correlated observations within ONE light-rail line, not independent treatments or replicated railway corridors.",
            "No train-driven seed dispersal, spore production, species diversity, flowering, regrowth after cutting, train gust or genetic adaptation measured.",
        ]
    }
    with tempfile.TemporaryDirectory() as root:
        road_dist,road_sources=build_n13_road_index(download_dir=Path(root))
        road_at={key:road_dist(point) for key,point in points.items()}
    result["N13_road_sources"]=road_sources
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-rail20m35m-N13-S2-source-v14"})
    baseline_days={}
    post={}
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="2",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        VSI_CACHE="TRUE",VSI_CACHE_SIZE="30000000",
        GDAL_CACHEMAX=256
    ):
        for day,scene_ids in BASELINE_DATES:
            try:
                items=[fetch_stac_scene(client,scene_id) for scene_id in scene_ids]
                for item in items:
                    if (item.get("properties") or {}).get("earthsearch:boa_offset_applied") is not False:
                        raise ValueError("2017 BOA flag not verified")
                    for band in ("red","nir"):
                        spec=(item["assets"][band].get("raster:bands") or [{}])[0]
                        if abs(float(spec.get("offset") or 0))>1e-8:
                            raise ValueError("2017 L2A nonzero BOA offset")
                vals,qa=choose_year_pixel_items(2017,points,items)
                baseline_days[day]=vals
                result["source_qa"][day]={"n_valid_pixels":len(vals),"QA":qa}
                print(json.dumps({
                    "2017_day":day,"scl_clear":len(vals),
                    "20m_same_day_pairs":sum(a in vals and b in vals for a,b in pairs_by_offset[20]),
                    "35m_same_day_pairs":sum(a in vals and b in vals for a,b in pairs_by_offset[35])
                }),flush=True)
            except Exception as e:
                result["errors"][day]=f"{type(e).__name__}: {str(e)[:250]}"
        for year,scene_ids in C1_SCENES.items():
            try:
                items=[fetch_collection1_item(client,scene_id) for scene_id in scene_ids]
                vals,qa=collection1_calibrated_pixels(points,items)
                post[year]={key:row["harmonized"] for key,row in vals.items()}
                result["source_qa"][str(year)]={"n_valid_pixels":len(vals),"QA":qa}
            except Exception as e:
                result["errors"][str(year)]=f"{type(e).__name__}: {str(e)[:250]}"
    for offset,pairs in pairs_by_offset.items():
        pre,dates,pair_dates=baseline_pair_first_clear(pairs,baseline_days)
        pre_values={key:value["harmonized"] for key,value in pre.items()}
        baseline_classes={}
        for (near,far),day in pair_dates.items():
            baseline_classes[near]=baseline_days[day][near]["scl_class"]
            baseline_classes[far]=baseline_days[day][far]["scl_class"]
        analysis=analyze_pairs(pairs,post,pre_values,baseline_classes,road_at)
        analysis["n_baseline_same_day_valid_pairs"]=sum(dates.values())
        analysis["2017_selected_date_counts"]=dates
        result["offsets"][f"{offset}m"]=analysis
        print(json.dumps({
            "offset":offset,
            "median_current_road_distance_m":analysis["2024_road_distance_near_median_m"],
            "eligible":analysis["eligibility_counts"],
            "NDVI":{
                year:{
                    gate:item[gate] for gate in (
                        "all","2017_both_vegetation",
                        "strict_veg_road_isolated_and_road_matched")
                } for year,item in analysis.get("pixel_NDVI_contrasts",{}).items()
            }
        },ensure_ascii=False),flush=True)
    result["created_utc"]=datetime.now(timezone.utc).isoformat()
    result["status"]=(
        "actually_observed_near_centerline_road_and_2017_scl"
        if any(x["n_baseline_same_day_valid_pairs"] for x in result["offsets"].values())
        else "no_valid_2017_rail_centerline_source_pairs"
    )
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/satellite/lrt_rail_centerline_20m_35m_road_scl_v14.json")
    a=p.parse_args()
    result=audit()
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"status":result["status"],
                      "road_parts":result["N13_road_sources"]["n_line_parts_in_spatial_index"],
                      "offset_counts":{k:v["eligibility_counts"] for k,v in result["offsets"].items()},
                      "source_errors":result["errors"]},ensure_ascii=False))
    if result["status"]=="no_valid_2017_rail_centerline_source_pairs":
        raise SystemExit("No real baseline source overlap for 20m/35m, no railway greenness inference")


if __name__=="__main__":
    main()
