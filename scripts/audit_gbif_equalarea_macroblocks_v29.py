#!/usr/bin/env python3
"""v2.9: independent macro-geographic 50 km EQUAL-AREA plant-data coverage.

Pre-frozen 24 geographic anchor blocks (3 each in 8 countries) chosen by
broad latitude/longitude spread, NEVER from plant or railway outcomes.
Each block consists of a 4x4 EPSG:6933 equal-area 50km square grid
(16 cells of 2500 km2); 384 cells in complete frame. Uses official
GBIF source COUNT + top dataset category FACET to assess source
opportunity for 2010-17 vs 2018-25.

This is a source-coverage audit, not a railway-ecological effect, species
incidence, or an official ecoregion classification. Country polygon and
land/ocean fraction are not yet verified, so low-count cells cannot be
treated as botanical absences or true non-rail controls.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import time

from pyproj import Transformer
from shapely.geometry import box
from shapely.ops import transform

from audit_gbif_macro_country_coverage_v18 import get_json

GBIF="https://api.gbif.org/v1/occurrence/search"
PERIODS=((2010,2017),(2018,2025))
GRID_CRS="EPSG:6933"
EACH_CELL_METERS=50000
BLOCK_CELLS_PER_SIDE=4
MIN_RECORDS=100
MIN_DATASET_CATEGORIES=3
XMIN_DIVISIONS=4

# Geographic source stratification, not claim of actual ecoregion or biome.
# Coordinates here are WGS84 (longitude, latitude), centered on broad
# geographical zones of eight countries with independent GRIIS checklists.
BLOCKS={
    "jp_tohoku":("JP",141.15,39.70),
    "jp_chubu":("JP",138.20,36.60),
    "jp_kyushu":("JP",130.70,32.80),
    "us_midwest":("US",-89.40,43.10),
    "us_southeast":("US",-78.60,35.80),
    "us_pacific":("US",-123.09,44.05),
    "ca_ontario":("CA",-75.70,45.42),
    "ca_prairie":("CA",-97.15,49.90),
    "ca_pacific":("CA",-123.12,49.26),
    "gb_yorkshire":("GB",-1.10,53.96),
    "gb_scotland":("GB",-4.22,57.48),
    "gb_south":("GB",-1.80,51.07),
    "fr_burgundy":("FR",5.04,47.32),
    "fr_mediterranean":("FR",3.88,43.61),
    "fr_atlantic":("FR",-1.55,47.22),
    "au_southeast":("AU",149.13,-35.28),
    "au_east":("AU",153.03,-27.47),
    "au_southwest":("AU",115.86,-31.95),
    "de_saxony":("DE",12.37,51.34),
    "de_blackforest":("DE",7.84,48.00),
    "de_hamburg":("DE",10.00,53.55),
    "nz_north":("NZ",175.30,-37.78),
    "nz_wellington":("NZ",174.78,-41.29),
    "nz_south":("NZ",170.50,-45.87),
}
# Fixed pilot (10/24 blocks): includes 8 countries, north/south climates and
# urban+rural surrounding cells, NOT optimized by the results.
PILOT=(
    "jp_tohoku","jp_kyushu","us_midwest","us_southeast",
    "ca_ontario","gb_yorkshire","fr_burgundy","au_southeast",
    "de_saxony","nz_north"
)


def validate_frame():
    if len(BLOCKS)!=24 or len({v[0] for v in BLOCKS.values()})!=8:
        raise ValueError("Frozen macro source block count changed")
    if set(PILOT)-set(BLOCKS) or len(PILOT)!=10:
        raise ValueError("Macro pilot block names changed")
    geo_to_proj=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    areas={}
    for name,(country,lon,lat) in BLOCKS.items():
        if len(country)!=2 or not (-179<lon<179 and -80<lat<80):
            raise ValueError(f"Invalid country / geography {name}")
        x,y=geo_to_proj.transform(lon,lat)
        half=EACH_CELL_METERS*BLOCK_CELLS_PER_SIDE/2.
        areas[name]=box(x-half,y-half,x+half,y+half)
    # Every freeze-frame block has independent, NON-OVERLAPPING grid
    # footprint within each country. Do not count neighboring blocks as
    # independent if their actual 200km grid windows overlap.
    for a,(country,_,_) in BLOCKS.items():
        for b,(other,_,_) in BLOCKS.items():
            if a>=b or country!=other:continue
            if areas[a].intersection(areas[b]).area>0:
                raise ValueError(f"Overlapping country macro blocks {a} {b}")
    return areas


def dense_polygon_in_epsg6933(cell):
    # Densify rectangle edges at 10km resolution *in equal-area CRS*
    # before WGS84 reprojected WKT used for GBIF polygon search.
    l,b,r,t=cell.bounds
    n=5
    xy=[]
    for i in range(n):
        xy.append((l+(r-l)*i/n,b))
    for i in range(n):
        xy.append((r,b+(t-b)*i/n))
    for i in range(n):
        xy.append((r-(r-l)*i/n,t))
    for i in range(n):
        xy.append((l,t-(t-b)*i/n))
    xy.append(xy[0])
    projection=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    lonlat=[projection.transform(x,y) for x,y in xy]
    if not all(math.isfinite(x) and math.isfinite(y) and -180<x<180 and -90<y<90 for x,y in lonlat):
        raise ValueError("Invalid projected geographic source cells")
    pts=",".join(f"{x:.6f} {y:.6f}" for x,y in lonlat)
    return "POLYGON(("+pts+"))"


def make_grid(block):
    if block not in BLOCKS:
        raise ValueError("Unregistered source block")
    _,lon,lat=BLOCKS[block]
    x,y=Transformer.from_crs(4326,GRID_CRS,always_xy=True).transform(lon,lat)
    step=EACH_CELL_METERS
    half=step*BLOCK_CELLS_PER_SIDE/2
    grid=[]
    for yy in range(BLOCK_CELLS_PER_SIDE):
        for xx in range(BLOCK_CELLS_PER_SIDE):
            cell=box(x-half+xx*step,y-half+yy*step,
                     x-half+(xx+1)*step,y-half+(yy+1)*step)
            grid.append({
                "cell_id":f"{block}_x{xx}_y{yy}",
                "country":BLOCKS[block][0],
                "epsg6933_equal_area_km2":round(cell.area/1e6,3),
                "gbif_geometry_wgs84":dense_polygon_in_epsg6933(cell),
            })
    if len(grid)!=16 or any(r["epsg6933_equal_area_km2"]!=2500. for r in grid):
        raise ValueError("Equal-area macro cell geometry changed")
    return grid


def query_params(cell,period):
    if period not in PERIODS:
        raise ValueError("Nonpredefined GBIF source time period")
    return {
        "country":cell["country"],"geometry":cell["gbif_geometry_wgs84"],
        "kingdomKey":6,"taxonRank":"SPECIES","occurrenceStatus":"PRESENT",
        "year":f"{period[0]},{period[1]}",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "facet":"datasetKey","facetLimit":4,
        "datasetKey.facetLimit":4,"limit":0,
    }


def parse_coverage(payload):
    if not isinstance(payload,dict):
        raise ValueError("Missing actual GBIF response JSON")
    count=payload.get("count")
    if isinstance(count,bool) or not isinstance(count,int) or count<0:
        raise ValueError("Missing or bad indexed GBIF occurrence count")
    entries=payload.get("facets")
    if not isinstance(entries,list):
        raise ValueError("Missing dataset facet source response")
    facets=[v for v in entries if isinstance(v,dict) and
            "".join(c for c in str(v.get("field","")).lower() if c.isalnum())=="datasetkey"]
    if count==0 and not facets:
        data=[]
    elif len(facets)==1 and isinstance(facets[0].get("counts"),list):
        data=facets[0]["counts"]
    else:
        raise ValueError("No source-verified GBIF DATASET_KEY facet")
    if len(data)>4:
        raise ValueError("Returned more source categories than requested")
    names=[]
    for r in data:
        if not isinstance(r,dict) or not str(r.get("name","")).strip():
            raise ValueError("Invalid GBIF dataset source facet")
        if isinstance(r.get("count"),bool) or not isinstance(r.get("count"),int) or r["count"]<0:
            raise ValueError("Invalid GBIF source dataset category count")
        names.append(str(r["name"]))
    if len(set(names))!=len(names):
        raise ValueError("Duplicate dataset category keys")
    return {
        "indexed_present_plant_record_count":count,
        "n_dataset_categories_at_least":len(data),
        "sufficient_records":count>=MIN_RECORDS,
        "sufficient_datasets":len(data)>=MIN_DATASET_CATEGORIES,
        "source_eligible":count>=MIN_RECORDS and len(data)>=MIN_DATASET_CATEGORIES,
    }


def run(blocks,source_api=get_json,sleep=time.sleep):
    validate_frame()
    if not blocks or any(b not in BLOCKS for b in blocks) or len(set(blocks))!=len(blocks):
        raise ValueError("Invalid macro block selection")
    import requests
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-macro-equal-area-GBIF-source-coverage-v2.9"})
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "study":"fixed 50km equal-area EPSG6933 regional source opportunity, NOT railway effect",
        "frozen_24_blocks":list(BLOCKS),
        "frozen_10_block_pilot":list(PILOT),
        "request_blocks":list(blocks),
        "country_origin":"predefined study windows, not random global or verified WWF ecoregions",
        "cell_area_km2":2500,
        "periods":[list(x) for x in PERIODS],
        "blocks":{},"source_errors":{},
    }
    for block in blocks:
        rows={}
        for cell in make_grid(block):
            item={"country":cell["country"],"source_cell_id":cell["cell_id"],
                  "cell_area_km2":cell["epsg6933_equal_area_km2"],
                  "period_coverage":{}}
            for period in PERIODS:
                key=f"{period[0]}_{period[1]}"
                try:
                    result=parse_coverage(source_api(client,GBIF,query_params(cell,period)))
                    item["period_coverage"][key]=result
                except Exception as exc:
                    out["source_errors"][f"{cell['cell_id']}_{key}"]=f"{type(exc).__name__}: {str(exc)[:170]}"
                    item["period_coverage"][key]={"source_error":True}
                sleep(.04)
            rows[cell["cell_id"]]=item
        complete=[c for c in rows.values() if len(c["period_coverage"])==2
                  and all(v.get("source_eligible") for v in c["period_coverage"].values())]
        def counts(period):
            return [c["period_coverage"][period].get("indexed_present_plant_record_count")
                    for c in rows.values()]
        out["blocks"][block]={
            "country":BLOCKS[block][0],
            "n_equal_area_cells":len(rows),
            "n_cells_both_period_recording_gate":len(complete),
            "n_early_records_sum_in_cells":sum(v or 0 for v in counts("2010_2017")),
            "n_late_records_sum_in_cells":sum(v or 0 for v in counts("2018_2025")),
            "n_period_complete_cells":sum(all("indexed_present_plant_record_count" in c["period_coverage"][p]
                                             for p in ("2010_2017","2018_2025"))
                                           for c in rows.values()),
            "cells":rows,
        }
        print(json.dumps({
            "block":block,"country":BLOCKS[block][0],
            "n_cells":len(rows),
            "eligible_both_periods":len(complete),
            "early_indexed_records":out["blocks"][block]["n_early_records_sum_in_cells"],
            "late_indexed_records":out["blocks"][block]["n_late_records_sum_in_cells"],
            "n_source_errors":len(out["source_errors"])
        },ensure_ascii=False),flush=True)
    out["complete"]=not out["source_errors"] and len(out["blocks"])==len(blocks)
    out["summary"]={
        "n_requested_geographical_blocks":len(blocks),
        "n_country_sources":len({BLOCKS[b][0] for b in blocks}),
        "n_equal_area_source_cells":len(blocks)*16,
        "n_cells_with_both_period_coverage_gate":sum(b["n_cells_both_period_recording_gate"] for b in out["blocks"].values()),
        "n_complete_source_period_bins":sum(b["n_period_complete_cells"]*2 for b in out["blocks"].values()),
    }
    out["interpretation_gates"]=[
        "Country and block names chosen by geographic anchors, not by train or plant effects; not yet independent ecoregion labels.",
        "Equal-area cell area is 2500 km2 in EPSG6933, but ocean cells and country-border fragments are not yet removed.",
        "All counts include kingdom Plantae not strictly vascular taxa and not introduced-only.",
        "Dataset count is a source category count, not number of independent survey protocols or observers.",
        "No speciesKey composition or rail/road topology measured for these larger cells in this source-only pilot.",
        "All geographic blocks must be checked for land fraction, true ecoregion, railroad/road active status before ecological effect tests."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--blocks",nargs="+",choices=list(BLOCKS),default=list(PILOT))
    p.add_argument("--output",default="results/audit/gbif_50km_10macroblocks_source_v29.json")
    a=p.parse_args()
    out=run(a.blocks)
    target=Path(a.output)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],
                      "errors":out["source_errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("50-km equal-area source audit incomplete; no ecological inference")


if __name__=="__main__":
    main()
