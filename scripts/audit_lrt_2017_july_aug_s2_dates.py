#!/usr/bin/env python3
"""Discover candidate true preconstruction July–August 2017 S2 L2A scenes.

The previous scene-cloud<20% filter excluded all except a few 2017 dates.
Scene-level cloud percentage can be high while the precise LRT pixel is
clear. Expand the fixed 2017 acquisition window 10 July–25 August and
scene-cloud<85%, but STILL require physical SCL cloud filtering per pixel.

This audit does not read NDVI or optimize any railway effect; fixed
date/quality screen prevents cherry-picking ecological response.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from datetime import date
import json
from pathlib import Path
import requests
from urllib.parse import urlparse

from audit_utsunomiya_sentinel_stac import BBOX,STAC

YEAR=2017
TARGET_DOY=213
FIRST="2017-07-10T00:00:00Z"
LAST="2017-08-25T23:59:59Z"
SCENE_CLOUD_LT=85
MAX_STAC_ITEMS=100


def summary_source(data):
    if data.get("type")!="FeatureCollection":
        raise ValueError("STAC query returned non-FeatureCollection")
    features=data.get("features") or []
    by_date=defaultdict(list)
    skip=[]
    for item in features:
        date_string=str((item.get("properties") or {}).get("datetime") or "")[:10]
        try:
            day=date.fromisoformat(date_string)
        except ValueError:
            skip.append("invalid_day")
            continue
        if not (date(2017,7,10)<=day<=date(2017,8,25)):
            raise ValueError("STAC date filtering ignored")
        assets=item.get("assets") or {}
        if not {"red","nir","scl"}.issubset(assets):
            skip.append("missing_red_nir_scl")
            continue
        meta=item.get("properties") or {}
        baseline=meta.get("s2:processing_baseline")
        band=(assets["red"].get("raster:bands") or [{}])[0]
        off=float(band.get("offset") or 0)
        applied=meta.get("earthsearch:boa_offset_applied")
        if abs(off)>0.000001 or applied is not False:
            skip.append("2017_offset_not_frozen_to_zero")
            continue
        by_date[date_string].append({
            "scene_id":item.get("id"),
            "cloud_cover_pct":meta.get("eo:cloud_cover"),
            "processing_baseline":baseline,
            "boa_offset_applied":applied,
        })
    groups=[]
    for ds,scenes in by_date.items():
        day=date.fromisoformat(ds)
        groups.append({
            "date":ds,
            "doy_distance_from_aug1":abs((day-date(2017,1,1)).days+1-TARGET_DOY),
            "n_tiles":len(scenes),
            "scene_ids":[x["scene_id"] for x in sorted(scenes,key=lambda x:x["scene_id"])],
            "cloud_pct":[x["cloud_cover_pct"] for x in scenes],
        })
    groups.sort(key=lambda x:(x["doy_distance_from_aug1"],x["date"]))
    return {
        "n_returned_source_items":len(features),
        "n_source_acquisition_dates":len(groups),
        "scene_dates_in_outcome_independent_order":groups,
        "excluded_source_reasons":dict((k,skip.count(k)) for k in set(skip)),
        "pagination_uncertain":any(x.get("rel")=="next" for x in data.get("links") or []),
        "conclusion":"Source date availability only, NO NDVI or biological railway effect measured",
    }


def collect_pages(client,initial_data,max_pages=20):
    """Follow exact STAC pagination links, respecting GET or POST method."""
    items=[]
    page=initial_data
    page_count=0
    visited=set()
    while True:
        if page.get("type")!="FeatureCollection":
            raise ValueError("invalid STAC pagination result")
        page_count+=1
        if page_count>max_pages:
            raise ValueError("More than 20 STAC source pages; incomplete inventory")
        for item in page.get("features") or []:
            sid=item.get("id")
            if sid in visited:
                raise ValueError("repeated STAC ID across page boundary")
            visited.add(sid)
            items.append(item)
        nxt=next((x for x in page.get("links") or []
                  if x.get("rel")=="next"),None)
        if nxt is None:
            break
        href=nxt.get("href") or STAC
        url=urlparse(href)
        if url.scheme!="https" or url.netloc!="earth-search.aws.element84.com":
            raise ValueError("untrusted STAC pagination link host")
        method=str(nxt.get("method") or "GET").upper()
        if method=="POST":
            body=nxt.get("body")
            if not isinstance(body,dict):
                raise ValueError("POST STAC pagination without body")
            response=client.post(href,json=body,timeout=(20,90))
        elif method=="GET":
            response=client.get(href,timeout=(20,90))
        else:
            raise ValueError("unsupported STAC pagination method")
        response.raise_for_status()
        page=response.json()
    return {"type":"FeatureCollection","features":items,"links":[]},page_count


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/lrt_2017_july_aug_sentinel_scl_candidate_dates.json")
    a=p.parse_args()
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-LRT-2017-preconstruction-SCL-source/1.2"})
    payload={
        "collections":["sentinel-2-l2a"],"bbox":BBOX,
        "datetime":f"{FIRST}/{LAST}",
        "limit":MAX_STAC_ITEMS,
        "query":{"eo:cloud_cover":{"lt":SCENE_CLOUD_LT}}
    }
    r=client.post(STAC,json=payload,timeout=(20,90))
    r.raise_for_status()
    complete,n_pages=collect_pages(client,r.json())
    results=summary_source(complete)
    results["n_source_pages_fetched"]=n_pages
    results["source"]=STAC
    results["window"]=[FIRST,LAST]
    results["cloud_cover_lt_pct"]=SCENE_CLOUD_LT
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "returned":results["n_returned_source_items"],
        "n_dates":results["n_source_acquisition_dates"],
        "top_dates":results["scene_dates_in_outcome_independent_order"][:12],
        "source_pages":n_pages,
        "paging":results["pagination_uncertain"],
    },ensure_ascii=False),flush=True)
    if results["n_source_acquisition_dates"]==0:
        raise SystemExit("No preconstruction July/August scenes meeting fixed source gates")


if __name__=="__main__":
    main()
