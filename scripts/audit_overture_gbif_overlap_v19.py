#!/usr/bin/env python3
"""v1.9 Overture rail-road × GBIF plant-taxon spatial source JOINT audit.

This is a source gate for an eventual MACROECOLOGICAL CELL ANALYSIS.
- Overture October-2026 global transportation is ODbL, not train counts.
- GBIF returns an actual area-specific indexed Plantae occurrence count
  and top-species facet (not a species-level complete assemblage).
- Neither source can be interpreted as train-carried seed dispersal.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import requests
import time

RELEASE="2026-10-21.0"
OVERTURE="s3://overturemaps-us-west-2/release/"+RELEASE+"/theme=transportation/type=segment/*"
GBIF="https://api.gbif.org/v1/occurrence/search"
# Outcome-blind, named first pilots spanning East Asia / Europe; city windows
# purposely distinct from their full countries and NOT representative.
REGIONS={
    "jp_sendai": ("JP", 140.77,38.16,141.01,38.37),
    "de_leipzig": ("DE",12.24,51.25,12.52,51.45),
    "kr_daejeon":("KR",127.27,36.24,127.57,36.48),
    "tw_hsinchu":("TW",120.84,24.68,121.14,24.95),
    "us_madison":("US",-89.57,42.98,-89.25,43.22),
    "ca_ottawa":("CA",-75.90,45.28,-75.52,45.53),
    "gb_york":("GB",-1.26,53.87,-0.92,54.08),
    "fr_dijon":("FR",4.86,47.20,5.18,47.45),
    "au_canberra":("AU",148.93,-35.48,149.29,-35.16),
    "nz_hamilton":("NZ",175.12,-37.89,175.42,-37.68),
}
RAIL_CLASSES=("standard_gauge","narrow_gauge","light_rail")
ROAD_CLASSES=("motorway","trunk","primary","secondary","tertiary","residential","unclassified")
SOURCE_GEOGRAPHY_WARNING="City windows are a source-availability pilot, NOT independent ecoregions or a rail-network connectivity test"


def validate_regions(regions=REGIONS):
    if len(regions)!=len(set(regions)) or len(regions)!=10:
        raise ValueError("Expected fixed 10 regions")
    country_set=set()
    for name,frame in regions.items():
        country,west,south,east,north=frame
        if not (len(country)==2 and country.isupper() and west<east and south<north
                and abs(west)<180 and abs(east)<180 and abs(south)<90 and abs(north)<90):
            raise ValueError(f"Bad region {name}")
        country_set.add(country)
    if len(country_set)!=10:
        raise ValueError("Country replication frame collapsed")


def sql_bbox(where_class,frame):
    """Use only frozen numeric bbox values; no raw user-supplied SQL."""
    c,w,s,e,n=frame
    for v in (w,s,e,n):
        if not math.isfinite(v):
            raise ValueError("Non-finite region extent")
    sql=f"""SELECT subtype, "class" AS route_class, COUNT(*) AS n_source_segments
       FROM read_parquet('{OVERTURE}', hive_partitioning=1)
       WHERE bbox.xmin <= {e} AND bbox.xmax >= {w}
       AND bbox.ymin <= {n} AND bbox.ymax >= {s}
       AND subtype IN ('rail','road')
       GROUP BY subtype, "class" ORDER BY subtype, route_class"""
    return sql


def gbif_query_params(frame,period=(2018,2025)):
    c,w,s,e,n=frame
    poly=f"POLYGON(({w} {s},{e} {s},{e} {n},{w} {n},{w} {s}))"
    return {
        "country":c,"geometry":poly,"kingdomKey":6,
        "taxonRank":"SPECIES","year":f"{period[0]},{period[1]}",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "facet":"speciesKey","speciesKey_facetLimit":50,
        "speciesKey_facetOffset":0,"limit":0
    }


def parse_gbif_facet(payload):
    if not isinstance(payload,dict) or not isinstance(payload.get("count"),int):
        raise ValueError("Missing GBIF actual indexed occurrence count")
    if payload["count"]<0:
        raise ValueError("Negative source count")
    facets=payload.get("facets")
    if not isinstance(facets,list):
        raise ValueError("Missing GBIF speciesKey facet list; no source incidence possible")
    fs=[f for f in facets if f.get("field")=="speciesKey"]
    if len(fs)!=1 or not isinstance(fs[0].get("counts"),list):
        raise ValueError("Malformed GBIF speciesKey facet")
    rows=fs[0]["counts"]
    parsed=[]
    for row in rows:
        if not isinstance(row,dict):
            raise ValueError("Malformed facet category")
        k=row.get("name")
        v=row.get("count")
        if k is None or isinstance(v,bool) or not isinstance(v,int) or v<0:
            raise ValueError("Invalid speciesKey facet")
        parsed.append({"gbif_speciesKey":str(k),"n_indexed_occurrences":v})
    if len(parsed)>50:
        raise ValueError("Unexpected GBIF facet limit overflow")
    # Facet top 50 is outcome/effort-biased; no species list completeness.
    return {
        "source_record_count":payload["count"],
        "top50_species_key_categories":parsed,
        "n_top_categories_returned":len(parsed),
        "exhaustive_species_inventory":False,
    }


def open_network_connection():
    import duckdb
    db=duckdb.connect(database=":memory:")
    db.execute("INSTALL httpfs; LOAD httpfs;")
    db.execute("SET s3_region='us-west-2';")
    db.execute("SET threads=2;")
    db.execute("SET memory_limit='4GB';")
    db.execute("SET enable_progress_bar=false;")
    return db


def query_network(db,frame):
    rows=db.execute(sql_bbox(None,frame)).fetchall()
    classes=[]
    for subtype,route_class,n in rows:
        if subtype not in ("road","rail") or not isinstance(n,int) or n<0:
            raise ValueError("Unexpected Overture source class")
        classes.append({"subtype":subtype,"class":route_class,"n_bbox_intersecting_segments":n})
    return {
        "source_classes":classes,
        "n_rail_segments_eligible_class":sum(x["n_bbox_intersecting_segments"]
            for x in classes if x["subtype"]=="rail" and x["class"] in RAIL_CLASSES),
        "n_road_segments_major_class":sum(x["n_bbox_intersecting_segments"]
            for x in classes if x["subtype"]=="road" and x["class"] in ROAD_CLASSES),
        "not_a_network_graph":True,
        "segments_only_not_clipped_length":True,
        "no_active_service_or_freight_inferred":True,
    }


def get_gbif(session,params):
    for attempt in range(4):
        try:
            r=session.get(GBIF,params=params,timeout=(20,80))
            r.raise_for_status()
            return parse_gbif_facet(r.json())
        except (requests.RequestException,ValueError):
            if attempt==3:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("GBIF query exhausted retry budget")


def run(names):
    validate_regions()
    if not names or any(n not in REGIONS for n in names) or len(names)!=len(set(names)):
        raise ValueError("Unexpected source region selection")
    report={
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "release":RELEASE,"overture_data_path":OVERTURE,
        "source_license":"ODbL, OpenStreetMap contributors / Overture Maps",
        "gbif_source":GBIF,
        "source_type":"geographical rail+road class and spatial GBIF taxon facet; not biological association",
        "all_region_names_frozen":list(REGIONS),
        "regions_requested":list(names),"results":{},"source_errors":{},
        "warning":SOURCE_GEOGRAPHY_WARNING,
    }
    db=None
    try:
        db=open_network_connection()
        sess=requests.Session()
        sess.headers.update({"User-Agent":"tetsu-macroecology-cross-source-rail-gbif-v1.9"})
        for name in names:
            outcome={"country":REGIONS[name][0]}
            try:
                outcome["gbif"]=get_gbif(sess,gbif_query_params(REGIONS[name]))
            except Exception as exc:
                report["source_errors"][f"{name}_gbif"]=f"{type(exc).__name__}: {str(exc)[:250]}"
            try:
                outcome["overture"]=query_network(db,REGIONS[name])
            except Exception as exc:
                report["source_errors"][f"{name}_overture"]=f"{type(exc).__name__}: {str(exc)[:250]}"
            report["results"][name]=outcome
            print(json.dumps({
                "region":name,"country":outcome["country"],
                "gbif_record_count":outcome.get("gbif",{}).get("source_record_count"),
                "top_species_facets":outcome.get("gbif",{}).get("n_top_categories_returned"),
                "rail_segments":outcome.get("overture",{}).get("n_rail_segments_eligible_class"),
                "major_road_segments":outcome.get("overture",{}).get("n_road_segments_major_class"),
                "n_source_errors":len(report["source_errors"])
            },ensure_ascii=False),flush=True)
    finally:
        if db:
            db.close()
    report["source_complete"]=(
        len(report["results"])==len(names) and not report["source_errors"]
        and all("gbif" in report["results"][n] and "overture" in report["results"][n] for n in names)
    )
    report["do_not_infer"]=[
        "GBIF occurrence record counts are not species richness, alien status, standardized sampling, or verified railway habitat.",
        "Top 50 GBIF species facet is incomplete and recording-effort biased; not the plant assemblage.",
        "Overture bbox-intersecting segments include portions outside the bounds, bridges and potentially nonoperating rails.",
        "No rail graph or network connectivity has been estimated in this pilot.",
        "Association and dispersal require independently checked biological sampling, geographic/effort controls, graph topology and replication."
    ]
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--regions",nargs="+",choices=list(REGIONS),default=["jp_sendai","de_leipzig"])
    p.add_argument("--output",default="results/audit/overture_gbif_two_region_source_v19.json")
    a=p.parse_args()
    r=run(a.regions)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"source_complete":r["source_complete"],
                      "n_regions":len(r["results"]),
                      "errors":r["source_errors"]},ensure_ascii=False),flush=True)
    if not r["source_complete"]:
        raise SystemExit("Cross-source railway × GBIF gate incomplete")


if __name__=="__main__":
    main()
