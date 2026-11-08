#!/usr/bin/env python3
"""2025 real Sentinel-2 L2A SCL/NDVI for independently road-separated rail sites.

Select sites FIRST, using only MLIT N02 conventional railway and N13
road/rail/station geometry in 3 distinct JP regions, never plant/NDVI.
Then assess current Sentinel-2 Collection1 band pixel availability and
SCL=4 vegetation on the near-track 35m pixel and the 1600m matched
road-only comparator. Output aggregate counts per public railway line,
not private site coordinates/rasters/road shapes. Not a biological
ground-verified verge, mowing, seed dispersal, adaptative or causal effect.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import date,datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import math
import requests
import zipfile
from urllib.parse import urlparse

import rasterio

from audit_japan_ordinary_rail_road_separation import (
    REGIONS,N02,MAX_ZIP_RAIL,download_source,analyze_region
)
from audit_lrt_2017_july_aug_s2_dates import collect_pages
from verify_utsunomiya_s2_c1_pixels import (
    collection1_calibrated_pixels
)

STAC="https://earth-search.aws.element84.com/v1/search"
COLLECTION="sentinel-2-c1-l2a"
YEAR=2025
DATE_START="2025-07-15T00:00:00Z"
DATE_END="2025-08-31T23:59:59Z"
TARGET_DOY=213 # August 1, independently of plant greenness
MAX_CLOUD_SCENE=85.
N_ACQUISITION_DATES=4
MAX_STAC_PAGES=5


def validate_stac_tile(item, year):
    props=item.get("properties") or {}
    day=str(props.get("datetime") or "")[:10]
    if not day.startswith(str(year)+"-") or not (
        date(year,7,15)<=date.fromisoformat(day)<=date(year,8,31)
    ):
        raise ValueError("STAC 2025 seasonal window filter ignored")
    if item.get("type")!="Feature":
        raise ValueError("STAC item schema drift")
    assets=item.get("assets") or {}
    if not {"red","nir","scl"}.issubset(assets):
        return False
    for band in ("red","nir","scl"):
        uri=urlparse(str(assets[band].get("href") or ""))
        if (uri.scheme!="https" or not uri.hostname or
            not uri.hostname.endswith(".amazonaws.com") or
            not uri.path.lower().endswith(".tif")):
            raise ValueError("Unexpected remote Sentinel COG asset host/path")
    for band in ("red","nir"):
        entries=assets[band].get("raster:bands") or []
        if len(entries)<1:
            return False
        b=entries[0]
        if (abs(float(b.get("scale") or 0)-.0001)>1e-8 or
            abs(float(b.get("offset") or 0)+.1)>1e-5):
            raise ValueError(
                f"Unverified C1 BOA scale/offset {item.get('id')}"
            )
    cloud=props.get("eo:cloud_cover")
    if not isinstance(cloud,(int,float)) or not 0<=cloud<MAX_CLOUD_SCENE:
        raise ValueError("STAC cloud cover source filter ignored")
    return True


def choose_source_dates(geojson_items, year=YEAR):
    by_day=defaultdict(list)
    skipped=Counter()
    for item in geojson_items:
        try:
            usable=validate_stac_tile(item,year)
        except Exception:
            raise
        if not usable:
            skipped["not_all_red_nir_scl_assets"]+=1
            continue
        ds=item["properties"]["datetime"][:10]
        by_day[ds].append(item)
    day_order=sorted(by_day,key=lambda ds:(
        abs((date.fromisoformat(ds)-date(year,1,1)).days+1-TARGET_DOY),
        min(item["properties"]["eo:cloud_cover"] for item in by_day[ds]),
        ds
    ))
    chosen=day_order[:N_ACQUISITION_DATES]
    return {
        ds:sorted(by_day[ds],key=lambda item:(
            item["properties"]["eo:cloud_cover"],item["id"]
        ))[:8]
        for ds in chosen
    },dict(skipped)


def sample_region_2025(region,private_sites,client):
    bbox=REGIONS[region]["bbox"]
    payload={
        "collections":[COLLECTION],"bbox":list(bbox),
        "datetime":f"{DATE_START}/{DATE_END}",
        "limit":100,"query":{"eo:cloud_cover":{"lt":MAX_CLOUD_SCENE}}
    }
    res=client.post(STAC,json=payload,timeout=(20,100))
    res.raise_for_status()
    complete,pages=collect_pages(client,res.json(),max_pages=MAX_STAC_PAGES)
    dates,skips=choose_source_dates(complete["features"])
    if not dates:
        raise ValueError("No verified Collection1 L2A imagery in region")
    points={}
    pairs=[]
    names={}
    for i,site in enumerate(private_sites):
        n=f"C{i:03d}n";f=f"C{i:03d}c"
        # Each site is a fixed source candidate, not observation-selected.
        points[n]=site["near_xy"];points[f]=site["road_only_xy"]
        pairs.append((n,f))
        names[n]=(site["route_name"],site["operator_name"])
    per_day={}
    failures={}
    for ds,items in dates.items():
        try:
            observations,qa=collection1_calibrated_pixels(points,items)
            per_day[ds]=observations
            print(json.dumps({
                "region":region,"date":ds,"available_spectral_pixels":len(observations),
                "n_candidate_pairs":len(pairs),
                "n_near_scl4":sum(
                    v.get("scl_class")==4 for k,v in observations.items()
                    if k.endswith("n")
                ),
                "n_road_control_scl4":sum(
                    v.get("scl_class")==4 for k,v in observations.items()
                    if k.endswith("c")
                ),
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            failures[ds]=f"{type(exc).__name__}: {str(exc)[:250]}"
    # Fail closed: if ALL acquisition dates fail at I/O/CRS, an
    # empty list of observations is UNKNOWN, NOT proof of no vegetation.
    if not per_day:
        raise ValueError(
            f"No successful 2025 Sentinel C1 dates in {region}; "
            f"source failures={failures}"
        )
    # For source eligibility, a pixel is observed *vegetation* if SCL=4
    # on ANY of the four independently selected source dates; '0' is
    # conditional on actual 2025 source coverage, not biological absence.
    outputs=defaultdict(lambda:{
        "n_geo_candidates":0,
        "near_scl4_on_at_least_1_preselected_date":0,
        "near_scl4_on_at_least_2_preselected_dates":0,
        "near_and_control_scl4_on_same_day":0,
        "near_and_control_cloud_clear_same_day":0,
        "no_quality_valid_near_on_any_date":0
    })
    for near,control in pairs:
        group=names[near]
        result=outputs[group]
        result["n_geo_candidates"]+=1
        nearveg=sum(
            data.get(near,{}).get("scl_class")==4 for data in per_day.values()
        )
        bothveg=sum(
            data.get(near,{}).get("scl_class")==4 and
            data.get(control,{}).get("scl_class")==4
            for data in per_day.values()
        )
        sameclear=sum(near in data and control in data
                       for data in per_day.values())
        nnear=sum(near in data for data in per_day.values())
        result["near_scl4_on_at_least_1_preselected_date"]+=int(nearveg>=1)
        result["near_scl4_on_at_least_2_preselected_dates"]+=int(nearveg>=2)
        result["near_and_control_scl4_on_same_day"]+=int(bothveg>=1)
        result["near_and_control_cloud_clear_same_day"]+=int(sameclear>=1)
        result["no_quality_valid_near_on_any_date"]+=int(nnear==0)
    summaries=[
        {"route_name":key[0],"operator_name":key[1],**v}
        for key,v in sorted(outputs.items(),key=lambda kv:(
            -kv[1]["near_scl4_on_at_least_2_preselected_dates"],
            -kv[1]["near_and_control_scl4_on_same_day"],
            kv[0]
        ))
    ]
    return {
        "region":region,
        "n_road_geo_candidate_pairs":len(private_sites),
        "n_stac_source_items_all_dates":len(complete["features"]),
        "n_stac_pages_fetched":pages,
        "selected_acquisition_dates":list(dates),
        "selected_source_scene_ids_by_date":{
            ds:[x["id"] for x in items] for ds,items in dates.items()
        },
        "n_selected_dates_real_pixel_success":len(per_day),
        "date_source_failures":failures,
        "invalid_asset_reasons":skips,
        "route_SCL_vegetation_coverage_no_raw_pixels":summaries,
        "sum_near_scl4_one_date":sum(
            x["near_scl4_on_at_least_1_preselected_date"] for x in summaries
        ),
        "sum_near_scl4_two_dates":sum(
            x["near_scl4_on_at_least_2_preselected_dates"] for x in summaries
        ),
        "sum_both_near_and_control_vegetated_same_date":sum(
            x["near_and_control_scl4_on_same_day"] for x in summaries
        ),
        "data_quality":"One-point Sentinel SCL20m and NDVI10m are mixed pixels; vegetation class != verified railway grass plots."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/japan_ordinary_rail_current_s2_vegetation_v15.json")
    a=p.parse_args()
    raw,n02=download_source(N02,MAX_ZIP_RAIL)
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise SystemExit("Official 2025 N02 invalid ZIP")
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-ordinary-rail-2025-C1-vegetation-source/1.5"})
    result={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "N02_source_sha256":n02["source_sha256"],
        "N02_reference":2025,
        "N13_reference":2024,
        "STAC":STAC,"Collection":COLLECTION,"spectral_year":YEAR,
        "fixed_selection":"Only N02 type 11/12 ordinary rail, N13 road >75m near; N02 station >300m, 1600m control road≤50m and all passenger rail >1000m; no NDVI selected sites",
        "fixed_scenes":"2025 July15-Aug31, sorted by date near August1 then scene cloud, first four acquisition dates independent of SCL and NDVI",
        "results":{},"errors":{},
        "limitations":[
            "Mere near35m rail/road-isolated C1 SCL4 is coarse 20m pixel vegetation, not survey-confirmed railway right-of-way plant cover.",
            "A road-only source site can have low SCL4 because its 10-20m mixed pixel is paved; no absence of plants is implied.",
            "2024 road map and 2025 Sentinel imagery do not verify vegetation or road locations historically.",
            "Multiple normal transects along one railway are not independent railway replication; current positive SCL status is not an ecological effect estimate.",
            "Conventional passenger N02 has no freight-only tracks, operational volume, train type or actual seed-spore transportation.",
            "No wind, mowing, herbicide, blooming, propagule arrival, mortality, genetic selection or NDVI changes estimated."
        ]
    }
    for region in REGIONS:
        try:
            summary,sites=analyze_region(
                region,raw,private_site_positions=True
            )
            if not sites:
                raise ValueError("No N02/N13 strict road-source candidate pairs")
            source=sample_region_2025(region,sites,client)
            result["results"][region]={
                "geo_n_N02_N13_candidates":len(sites),
                "n_ordinary_rail_route_groups":summary["n_ordinary_route_groups_over_2p5km"],
                "satellite":source
            }
            print(json.dumps({
                "region":region,
                "source_geo_road_separated":len(sites),
                "near_SCL4_2025_on_1_date":source["sum_near_scl4_one_date"],
                "near_SCL4_2025_on_2_dates":source["sum_near_scl4_two_dates"],
                "near_and_road_control_both_SCL4_same_date":
                    source["sum_both_near_and_control_vegetated_same_date"],
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            result["errors"][region]=f"{type(exc).__name__}: {str(exc)[:350]}"
            print(json.dumps({
                "region":region,"source_error":result["errors"][region]
            }),flush=True)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    if not result["results"]:
        raise SystemExit("No road-separated site with actual 2025 Sentinel source audit")


if __name__=="__main__":
    main()
