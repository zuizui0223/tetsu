#!/usr/bin/env python3
"""Availability of S2 surface-reflectance scenes around 2023 Utsunomiya LRT.

Uses public Element84 Earth Search STAC catalogue; NO imagery downloaded.
Preconstruction 2017, construction 2018/2022, operating 2024/25.
A scene-wide cloud percentage DOES NOT guarantee clear pixels over LRT.
Even repeated NDVI would show vegetative cover, NOT species, seeds, wind
or mowing mechanism. GIS N05 is noncommercial; no licensed line pixels
or source geometry leave the pipeline.

This is a source-coverage audit only.
"""
from __future__ import annotations
import argparse
import collections
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import requests

STAC="https://earth-search.aws.element84.com/v1/search"
STAC_COLLECTION="sentinel-2-l2a"
# Coarse derived study envelope of actual 2023 Utsunomiya LRT from
# N05 railway GIS, audited in run #42438821083 (not source line geometry).
BBOX=[139.84,36.49,140.07,36.63]
YEAR_PHASES={
    2017:"preconstruction_source_screen",
    2018:"construction_started_june",
    2019:"construction",
    2020:"construction",
    2021:"construction",
    2022:"construction_plus_trials_november",
    2023:"opening_august_transition_excluded",
    2024:"operating",
    2025:"operating"
}
MONTHS=(6,7,8,9)
MAX_ITEMS_PER_QUERY=100
SCENE_CLOUD_PCT_LT=20


def stac_query_params(year):
    if year not in YEAR_PHASES:
        raise ValueError("year not in predeclared S2 season window")
    return {
        "collections":[STAC_COLLECTION],
        "bbox":BBOX,
        "datetime":f"{year}-06-01T00:00:00Z/{year}-09-30T23:59:59Z",
        "limit":MAX_ITEMS_PER_QUERY,
        "query":{"eo:cloud_cover":{"lt":SCENE_CLOUD_PCT_LT}},
    }


def summarize_items(year, response):
    if not isinstance(response,dict) or response.get("type")!="FeatureCollection":
        raise ValueError("Earth Search response not a STAC FeatureCollection")
    rows=response.get("features")
    if not isinstance(rows,list):
        raise ValueError("Earth Search missing list of STAC items")
    scenes=[]
    source_keys=collections.Counter()
    dates=set()
    for row in rows:
        identity=row.get("id")
        props=row.get("properties") or {}
        when=str(props.get("datetime") or "")
        if not identity or not when.startswith(str(year)+"-"):
            raise ValueError("STAC returned unexpected item year or missing ID")
        if len(when)<10 or not 6<=int(when[5:7])<=9:
            raise ValueError("STAC returned item outside fixed summer season")
        cloud=props.get("eo:cloud_cover")
        if not isinstance(cloud,(int,float)) or not 0<=cloud<SCENE_CLOUD_PCT_LT:
            raise ValueError(f"Unverified STAC cloud fraction {cloud}")
        assets=row.get("assets") or {}
        source_keys.update(assets.keys())
        dates.add(when[:10])
        # We do not download COGs; verify required spectral and cloud assets.
        keys=set(assets)
        compatible=(
            ({"red","nir"}<=keys or {"red","nir08"}<=keys or
             {"B04","B08"}<=keys) and
            (bool(keys&{"scl","SCL","cloud-mask"}))
        )
        scenes.append({"scene_id":identity,"date":when[:10],
                       "cloud_fraction":cloud,"available_for_ndvi":compatible})
    matched=response.get("context",{}).get("matched")
    # The page may only contain the first 100 scenes: never declare it
    # an exhaustive seasonal inventory when next-link or matched>limit.
    next_links=response.get("links") or []
    has_next=any(x.get("rel")=="next" for x in next_links)
    return {
        "year":year,"phase":YEAR_PHASES[year],
        "returned_stac_items":len(scenes),
        "number_of_distinct_acquisition_dates":len(dates),
        "n_items_with_red_nir_and_cloud_assets":sum(x["available_for_ndvi"] for x in scenes),
        "is_truncated_or_has_next":bool(has_next or (isinstance(matched,int) and matched>len(scenes))),
        "minimum_cloud_cover_scene_pct":min(
            [x["cloud_fraction"] for x in scenes],default=None
        ),
        "stac_asset_keys_observed":dict(source_keys.most_common(18)),
        "scene_quality_caveat":"Scene-level cloud fraction not local pixel cloud clearance; actual ROI analysis requires SCL pixel mask, equal season, acquisition geometry and background roads."
    }


def run():
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-railway-satellite-source-audit/1.1 (https://github.com/zuizui0223/tetsu)"})
    by_year=[]
    errors={}
    for year in YEAR_PHASES:
        try:
            q=stac_query_params(year)
            res=client.post(STAC,json=q,timeout=(20,75))
            res.raise_for_status()
            summary=summarize_items(year,res.json())
            by_year.append(summary)
            print(json.dumps({
                "year":year,"phase":summary["phase"],
                "scenes":summary["returned_stac_items"],
                "days":summary["number_of_distinct_acquisition_dates"],
                "red_nir_scl_ready":summary["n_items_with_red_nir_and_cloud_assets"],
                "capped":summary["is_truncated_or_has_next"]
            }),flush=True)
        except Exception as exc:
            errors[str(year)]=f"{type(exc).__name__}: {str(exc)[:240]}"
    return {
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "STAC_source":STAC,
        "collection":STAC_COLLECTION,
        "source_spatial_envelope_round_2decimal":BBOX,
        "source_bbox_is_derived_from_N05_with_noncommercial_rights":True,
        "cloud_cover_lt_pct":SCENE_CLOUD_PCT_LT,
        "months":list(MONTHS),
        "season_years":list(YEAR_PHASES),
        "years":by_year,
        "source_errors":errors,
        "scope":"satellite IMAGE AVAILABILITY only, no NDVI measured, no claim of vegetation cover or railway effect",
        "ecological_warnings":[
            "LRT construction began 2018; 2018-2022 is not an untreated causal baseline.",
            "Tram is mostly road-adjacent so Sentinel 10-m mixed pixels may capture buildings, pavement and road trees, not railway weeds.",
            "Cloud percentage applies to broad satellite scene, not exactly the corridor; physical SCL pixel masking essential.",
            "NDVI (if analyzed) measures greenness/plant cover, not species, seeds, spore cones, flowering, grass-cutting dates, train gusts or local genetic adaptation.",
            "Even reproducible NDVI change requires climate, land-cover, development, illumination, acquisition-season and matched road-only control."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/utsunomiya_lrt_sentinel_stac_2017_25.json")
    a=p.parse_args()
    result=run()
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    if not result["years"]:
        raise SystemExit("No Sentinel2 source results; do not claim vegetation change")
    print(json.dumps({"successful_years":len(result["years"]),
                      "errors":result["source_errors"]},ensure_ascii=False))


if __name__=="__main__":
    main()
