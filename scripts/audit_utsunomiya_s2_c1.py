#!/usr/bin/env python3
"""Audit Element84 reprocessed Sentinel-2 Collection 1 before choosing NDVI source.

The legacy sentinel-2-l2a COG collection has contradictory BOA correction
metadata on 2024/25 scenes. Element84 recommends the replacement
sentinel-2-c1-l2a with consistent radiometry. This audit searches all source
items on the frozen Utsunomiya LRT footprint, years 2017, 2024, 2025,
without downloading pixels or selecting for greenness.

2017 Collection 1 reprocessing is incomplete; absence is source missing,
NOT absence of vegetation.
"""
from __future__ import annotations
import argparse
import datetime as dt
import json
from pathlib import Path
import requests

from audit_utsunomiya_sentinel_stac import STAC,BBOX

COLLECTION="sentinel-2-c1-l2a"
YEARS=(2017,2024,2025)
CLOUD_LIMIT=20
MAX_ITEMS=100
TARGET_DOY={2017:188,2024:217,2025:211} # original legacy fixed scene dates


def query(year):
    return {
        "collections":[COLLECTION],"bbox":BBOX,
        "datetime":f"{year}-06-01T00:00:00Z/{year}-09-30T23:59:59Z",
        "limit":MAX_ITEMS,"query":{"eo:cloud_cover":{"lt":CLOUD_LIMIT}}
    }


def summary(year,data):
    if data.get("type")!="FeatureCollection":
        raise ValueError("Collection 1 STAC source invalid")
    features=data.get("features") or []
    items=[]
    for f in features:
        if not all(x in (f.get("assets") or {}) for x in ("red","nir","scl")):
            continue
        props=f.get("properties") or {}
        date=str(props.get("datetime") or "")
        if not date.startswith(str(year)+"-"):
            raise ValueError("C1 server returned wrong year")
        try:
            doy=(dt.date.fromisoformat(date[:10])-dt.date(year,1,1)).days+1
        except ValueError as e:
            raise ValueError("Invalid C1 sensing date") from e
        assets=f["assets"]
        band_specs={}
        for key in ("red","nir"):
            bands=assets[key].get("raster:bands") or []
            band_specs[key]=bands[0] if bands else None
        items.append({
            "item_id":f["id"],"date":date[:10],
            "cloud_cover_pct":props.get("eo:cloud_cover"),
            "processing_baseline":props.get("s2:processing_baseline"),
            "boa_offset_applied":props.get("earthsearch:boa_offset_applied"),
            "bands":band_specs,
            "doy_distance_from_original_legacy_sensing_date":abs(doy-TARGET_DOY[year])
        })
    items.sort(key=lambda x:(
        x["doy_distance_from_original_legacy_sensing_date"],
        x["cloud_cover_pct"] if isinstance(x["cloud_cover_pct"],(int,float)) else 1000,
        x["item_id"]
    ))
    return {
        "year":year,"n_stac_items_returned":len(features),
        "n_items_with_red_nir_scl":len(items),
        "candidate_scenes_nearest_original_date":items[:12],
        "has_additional_paging":any(x.get("rel")=="next" for x in data.get("links") or []),
        "no_2017_reprocessed_data_is_not_a_biological_absence":True,
    }


def run():
    s=requests.Session()
    s.headers.update({"User-Agent":"tetsu-reprocessed-s2-collection1-calibration/1.2"})
    result={
        "created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "STAC":STAC,"collection":COLLECTION,
        "source":"reprocessed and harmonized Sentinel-2 Collection 1 L2A",
        "prior_collection":"sentinel-2-l2a",
        "years":{},"failures":{},
        "source_caveat":"2017 C1 reprocessing may have very sparse items, do not fill missing observations with newer photos.",
        "scientific_caveat":"No source scenes were selected by NDVI or railway flowering outcome; no observed NDVI output.",
    }
    for year in YEARS:
        try:
            r=s.post(STAC,json=query(year),timeout=(20,90))
            r.raise_for_status()
            out=summary(year,r.json())
            result["years"][str(year)]=out
            print(json.dumps({
                "year":year,
                "returned":out["n_stac_items_returned"],
                "ready":out["n_items_with_red_nir_scl"],
                "nearest":out["candidate_scenes_nearest_original_date"][:5]
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            result["failures"][str(year)]=f"{type(exc).__name__}: {str(exc)[:320]}"
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/utsunomiya_s2_collection1_source.json")
    args=p.parse_args()
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    report=run()
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    if not report["years"]:
        raise SystemExit("No C1 L2A response; no harmonized satellite result")


if __name__=="__main__":
    main()
