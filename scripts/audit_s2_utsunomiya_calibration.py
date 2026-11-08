#!/usr/bin/env python3
"""Inspect real S2 COG radiometric calibration BEFORE comparing time periods.

Earth Search Sentinel-2 L2A STAC metadata can expose both per-band
raster:bands offset and earthsearch:boa_offset_applied. These may be
incompatible. Do NOT silently assume raw DN/10000 is harmonized or
that apparent seasonal NDVI is a true ecological before/after change.

Runs on public source, no images downloaded, no proprietary N05 geometry.
"""
from __future__ import annotations
import argparse
import datetime as dt
import json
from pathlib import Path
from urllib.parse import urlparse
import requests

from audit_utsunomiya_sentinel_stac import STAC,stac_query_params

YEARS=(2017,2024,2025)
TARGET_DOY=213


def calendar_distance(when: str) -> int:
    d=dt.date.fromisoformat(when[:10])
    return abs((d-dt.date(d.year,1,1)).days+1 - TARGET_DOY)


def pick_candidates(features, *, top=4):
    """Rank by seasonal distance first, then cloud, then stable ID.

    Do not select on NDVI or observed greenness.
    """
    candidates=[]
    for item in features:
        assets=item.get("assets") or {}
        if not {"red","nir","scl"}.issubset(assets):
            continue
        props=item.get("properties") or {}
        cloud=props.get("eo:cloud_cover")
        date=props.get("datetime")
        if not isinstance(cloud,(float,int)) or not date:
            continue
        candidates.append((calendar_distance(date),float(cloud),str(item["id"]),item))
    candidates.sort(key=lambda x:x[:3])
    return [x[3] for x in candidates[:top]]


def asset_meta(asset):
    url=str(asset.get("href") or "")
    parts=urlparse(url)
    if parts.scheme!="https" or not parts.netloc.endswith(("amazonaws.com","element84.com")):
        raise ValueError("Unexpected public Sentinel COG asset host")
    return {
        "url_host":parts.netloc,
        "asset_file_basename":parts.path.rsplit("/",1)[-1],
        "raster_bands":asset.get("raster:bands"),
        "eo_bands":asset.get("eo:bands"),
        "proj_shape":asset.get("proj:shape"),
        "proj_transform":asset.get("proj:transform"),
        "type":asset.get("type"),
    }


def summarize_item(item):
    p=item.get("properties") or {}
    metadata={
        "item_id":item["id"],
        "date":p.get("datetime"),
        "platform":p.get("platform"),
        "cloud_cover_pct":p.get("eo:cloud_cover"),
        "s2_processing_baseline":p.get("s2:processing_baseline"),
        "earthsearch_boa_offset_applied":p.get("earthsearch:boa_offset_applied"),
        "proj_epsg":p.get("proj:epsg"),
        "mgrs_tile":p.get("s2:mgrs_tile"),
        "footprint_bbox":item.get("bbox"),
        "band_metadata":{
            key:asset_meta(item["assets"][key]) for key in ("red","nir","scl")
        }
    }
    return metadata


def run():
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-rail-ndvi-calibration-source/1.2"})
    reports={}
    for year in YEARS:
        response=client.post(STAC,json=stac_query_params(year),timeout=(20,90))
        response.raise_for_status()
        data=response.json()
        items=pick_candidates(data.get("features") or [])
        reports[str(year)]={
            "returned_item_count":len(data.get("features") or []),
            "candidate_item_calibration":[summarize_item(x) for x in items],
            "has_more_stac_pages":any(x.get("rel")=="next" for x in data.get("links") or []),
        }
        print(json.dumps({
            "year":year,
            "candidate_items":[{
                "id":x["item_id"],
                "date":x["date"],
                "tile":x["mgrs_tile"],
                "baseline":x["s2_processing_baseline"],
                "boa_offset_applied":x["earthsearch_boa_offset_applied"],
                "red_raster_bands":x["band_metadata"]["red"]["raster_bands"],
                "nir_raster_bands":x["band_metadata"]["nir"]["raster_bands"],
            } for x in reports[str(year)]["candidate_item_calibration"]]
        },ensure_ascii=False),flush=True)
    return {
        "source":STAC,"collection":"sentinel-2-l2a",
        "target_doy":TARGET_DOY,
        "selection":"closest to day 213, then lower scene cloud, then stable item ID",
        "years":reports,
        "source_only":True,
        "interpretation":"Asset scale/offset and BOA application may disagree; verify before any cross-era NDVI comparison.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/s2_utsunomiya_radiometry_metadata.json")
    a=p.parse_args()
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    report=run()
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")


if __name__=="__main__":
    main()
