#!/usr/bin/env python3
"""v7.5: check FRA CROSSING Form 71 train traffic opportunity in 3 FROZEN
US 200km macro blocks, 48 identical EPSG:6933 50km cells.

Source: public ArcGIS hosted mirror of FRA Crossing Inventory, with live
source schema and page-count conservation. Report crossing-level
number of trains/day and actual source year only as OBSERVATION-SITE
summaries. This is not line-level 2026 freight train count, freight
tonne-km, actual seed movement, slipstream wind, or mowing.

NO raw crossing IDs, geolocations, personal details, or precise railway
positions are retained in the output. An inventory crossing can include
both passenger and freight; distinguish categorical service from total
through train counts, never subtract a passenger count to invent freight.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from statistics import median
import json,math,time,requests
from shapely.geometry import Point, box
from shapely import wkt
from shapely.ops import unary_union
from audit_gbif_equalarea_macroblocks_v29 import make_grid

SOURCE = "https://services1.arcgis.com/4yjifSiIG17X0gW4/arcgis/rest/services/FRA_Crossing_Inventory_Form_71_Current/FeatureServer/0"
BLOCKS = ("us_midwest","us_southeast","us_pacific")
EXPECTED_CELLS = 16
MAX_SOURCE_ROWS_PER_BLOCK = 20000
PAGE_SIZE = 1000
MAX_PAGES = 21
FIELDS = (
    "ObjectId","Crossing_ID","Crossing_Closed","ReportStatus",
    "Latitude","Longitude","Total_Daylight_Thru_Trains",
    "Total_Nighttime_Thru_Trains","Total_Switching_Trains",
    "Total_Thru_Trains","Trains_Per_Week_Captured_Year",
    "Type_Of_Train_Service_IDs",
    "Maximum_Timetable_Speed","Typical_Maximum_Speed_Over_Crossing",
    "Reporting_Railroad_Freight","Reporting_Railroad_Passenger",
    "State_Name"
)

def api_get(session,suffix,params):
    for retry in range(3):
        try:
            r=session.get(SOURCE+suffix,params={"f":"json",**params},timeout=(15,90))
            r.raise_for_status()
            out=r.json()
            if not isinstance(out,dict) or out.get("error"):
                raise ValueError("ArcGIS source returned incomplete/invalid payload")
            return out
        except requests.RequestException:
            if retry==2:raise
            time.sleep(1+retry)
    raise ValueError("Unreachable source retry state")

def frozen_grid(block):
    rows=make_grid(block)
    if block not in BLOCKS or len(rows)!=EXPECTED_CELLS:
        raise ValueError("Not a predeclared US 50km macro source block")
    cells={v["cell_id"]:wkt.loads(v["gbif_geometry_wgs84"]) for v in rows}
    footprint=unary_union(list(cells.values()))
    return cells,footprint.bounds

def query_box(bounds):
    return {"xmin":bounds[0],"ymin":bounds[1],"xmax":bounds[2],
            "ymax":bounds[3],"spatialReference":{"wkid":4326}}

def query_params(bounds):
    return {"where":"1=1",
            "geometry":json.dumps(query_box(bounds),separators=(",",":")),
            "geometryType":"esriGeometryEnvelope",
            "spatialRel":"esriSpatialRelIntersects","inSR":4326}

def fetch_features(session,block):
    cells,bounds=frozen_grid(block)
    q=query_params(bounds)
    meta=api_get(session,"",{})
    known={v.get("name") for v in meta.get("fields",[])}
    if meta.get("objectIdField")!="ObjectId" or any(f not in known for f in FIELDS):
        raise ValueError("Public 2026 FRA mirror source field contract changed")
    if meta.get("geometryType")!="esriGeometryPoint":
        raise ValueError("FRA crossing record source not spatial point features")
    count=api_get(session,"/query",{**q,"returnCountOnly":"true"})
    total=count.get("count")
    if type(total)!=int or total<0 or total>MAX_SOURCE_ROWS_PER_BLOCK:
        raise ValueError(f"Source bounding crossing count above 20k resource contract: {block}, {total}")
    records=[]
    ids=set()
    for page in range(math.ceil(total/PAGE_SIZE)):
        response=api_get(session,"/query",{**q,"outFields":",".join(FIELDS),
                            "returnGeometry":"true","outSR":4326,
                            "orderByFields":"ObjectId ASC",
                            "resultOffset":page*PAGE_SIZE,"resultRecordCount":PAGE_SIZE})
        batch=response.get("features")
        if not isinstance(batch,list) or len(batch)>PAGE_SIZE:
            raise ValueError("FRA spatial source page incomplete/oversized")
        for item in batch:
            attrs=item.get("attributes")
            if not isinstance(attrs,dict) or any(k not in attrs for k in FIELDS):
                raise ValueError("Incomplete FRA source feature attributes")
            oid=attrs["ObjectId"]
            if type(oid)!=int or oid in ids:
                raise ValueError("Duplicate or invalid ArcGIS crossing OBJECTID")
            ids.add(oid)
        records.extend(batch)
    if len(records)!=total or len(ids)!=total:
        raise ValueError("Actual FRA query page rows differ from independent source count")
    return cells, records, {
        "n_source_mirror_rows_in_200km_bbox":total,
        "n_source_arcgis_pages":math.ceil(total/PAGE_SIZE),
        "source_metadata_lastEditDate_ms":meta.get("editingInfo",{}).get("lastEditDate"),
        "source_schema_n_fields":len(meta.get("fields",[])),
        "source_page_counts_reconciled":True,
        "bounded_source_region_complete":True,
    }

def numeric_count(value,maximum):
    if type(value)==int and 0<=value<=maximum:return value
    return None

def year_of_count(value):
    if type(value) in (int,str) and str(value).isdigit():
        n=int(value)
        return n if 1970<=n<=2026 else None
    return None

def source_categories(value):
    if value is None:return ()
    if not isinstance(value,str):return ()
    import re
    return tuple(sorted(set(re.findall(r"(?<![0-9])[0-9]{2}(?![0-9])",value))))

def assign_cell(cells,attrs,geometry):
    lat,lon=attrs.get("Latitude"),attrs.get("Longitude")
    if any(type(z) not in (int,float) or not math.isfinite(z) for z in (lat,lon)):
        return None,"missing_coordinates"
    if not (-180<lon<180 and -90<lat<90):
        return None,"invalid_coordinates"
    if geometry and all(type(geometry.get(k)) in (int,float) for k in ("x","y")):
        # ArcGIS requested outSR=4326. Declared original degrees must agree.
        if abs(geometry["x"]-lon)>.02 or abs(geometry["y"]-lat)>.02:
            return None,"coordinate_source_disagreement"
    p=Point(lon,lat)
    hits=[k for k,v in cells.items() if v.covers(p)]
    if len(hits)!=1:return None,"outside_or_boundary_ambiguous"
    return hits[0],None

def summarise_block(cells,records):
    by_cell={cid:[] for cid in cells}
    excluded=Counter();seen=set();raw_types=Counter();closed=Counter()
    for item in records:
        attr=item["attributes"]
        cid=attr.get("Crossing_ID")
        if not isinstance(cid,str) or not cid.strip():
            excluded["missing_crossing_id"]+=1
            continue
        # FRA Crossing_ID is source-unique in the hosted current layer.
        if cid in seen:
            excluded["duplicate_crossing_id"]+=1
            continue
        seen.add(cid)
        which,reason=assign_cell(cells,attr,item.get("geometry") or {})
        if reason:
            excluded[reason]+=1
            continue
        by_cell[which].append(attr)
        raw_types[str(attr.get("Type_Of_Train_Service_IDs"))]+=1
        closed[str(attr.get("Crossing_Closed"))]+=1
    output={}
    for cell,rows in by_cell.items():
        yr=Counter();nreport=0;through=[];freight_tag=0;pass_tag=0
        speed=[];supported=[];unverified=0
        for r in rows:
            year=year_of_count(r.get("Trains_Per_Week_Captured_Year"))
            if year is None:yr["unknown"]+=1
            else:yr[str(year)]+=1
            # Neither "freight operator" nor service code indicates
            # an observed freight count on this route.
            tags=source_categories(r.get("Type_Of_Train_Service_IDs"))
            if "11" in tags: freight_tag+=1
            if any(x in tags for x in ("12","13","14","15")):pass_tag+=1
            if year is None or not 2018<=year<=2026:
                unverified+=1
                continue
            day=numeric_count(r.get("Total_Daylight_Thru_Trains"),500)
            night=numeric_count(r.get("Total_Nighttime_Thru_Trains"),500)
            if day is not None and night is not None:
                nreport+=1;through.append(day+night)
                supported.append(int(day+night>0))
            sp=numeric_count(r.get("Typical_Maximum_Speed_Over_Crossing"),200)
            if sp is not None:speed.append(sp)
        output[cell]={
            "n_distinct_crossing_inventory_records":len(rows),
            "n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields":nreport,
            "n_crossings_with_positive_reported_daily_through_train_counts":sum(supported),
            "median_reported_total_daily_through_trains_at_eligible_crossings":median(through) if through else None,
            "median_reported_typical_maximum_train_speed_mph_at_year_eligible_crossings":median(speed) if speed else None,
            "n_crossings_reporting_freight_service_type_NOT_freight_train_counts":freight_tag,
            "n_crossings_reporting_passenger_service_type_NOT_passenger_only_control":pass_tag,
            "reported_train_count_year_histogram":dict(yr),
            "n_missing_or_outside_2018to2026_year":unverified,
            "NOT_active_freight_volume_or_measured_gust_intensity":True
        }
    return {
        "by_original_50km_cells":output,
        "n_in_cell_source_crossing_records":sum(len(v) for v in by_cell.values()),
        "n_cells_with_any_record":sum(bool(rows) for rows in by_cell.values()),
        "n_cells_with_year_eligible_daily_train_counts":sum(x["n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields"]>0 for x in output.values()),
        "source_exclusion_diagnostics":dict(excluded),
        "source_type_service_raw_category_histogram":dict(raw_types),
        "source_crossing_closed_raw_category_histogram":dict(closed),
        "DO_NOT_TREAT_CROSSINGS_AS_INDEPENDENT_RAIL_SEGMENTS":True,
    }

def run(session):
    result={
        "complete":True,"run_utc":datetime.now(timezone.utc).isoformat(),
        "source_service":SOURCE,"source":"2026 publicly ArcGIS-hosted FRA Form 71 Crossing Inventory mirror",
        "frozen_US_geographic_blocks":list(BLOCKS),"n_frozen_equal_area_50km_cells":len(BLOCKS)*16,
        "study":"Actual reported crossing-level train movement/source-year coverage opportunity; not macro ecological railway effect",
        "regions":{},
        "CANNOT_CLAIM_WIND_STRESS_SEED_FLOW_FREIGHT_TRANSPORT_OR_LOCAL_ADAPTATION":True}
    for block in BLOCKS:
        cells,records,audit=fetch_features(session,block)
        s=summarise_block(cells,records)
        result["regions"][block]={**audit,**s}
        print(json.dumps({"block":block,"source_bbox_records":audit["n_source_mirror_rows_in_200km_bbox"],
                          "records_in_frozen_cells":s["n_in_cell_source_crossing_records"],
                          "cells_with_eligible_train_counts":s["n_cells_with_year_eligible_daily_train_counts"],
                          "source_exclusions":s["source_exclusion_diagnostics"]}),flush=True)
    result["limitations"]=[
        "Hosted public ArcGIS source is a copy of FRA Crossing Inventory, and version/date/schema source completeness and true official-source synchronization must be verified independently.",
        "FRA reporting year can differ by crossing; current snapshot cannot be treated as 2010-17 railway activity, annual panel or historical freight train frequency.",
        "Crossing locations are road-rail intersections and their density is strongly road-system dependent; not a representative route-wide sample.",
        "Total train counts include passenger and freight, with no freight-only count; 'freight' type is categorical not a tonnage or movement rate.",
        "Typical speed is not train-induced gust speed; mowing/herbicide regime and plant regrowth are unmeasured.",
        "This is source availability across 3 US macro geographic blocks; GBIF plant relationships must be tested separately with source/geographic controls and independent holdouts."
    ]
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True)
    args=p.parse_args()
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-US-FRA-crossing-source-audit-v7.5"})
    out=run(session)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":True,"n_regions":len(out["regions"]),
                      "n_total_source_rows":sum(x["n_source_mirror_rows_in_200km_bbox"] for x in out["regions"].values()),
                      "no_plant_or_freight_or_adaptation_claim":True}),flush=True)

if __name__=="__main__":main()
