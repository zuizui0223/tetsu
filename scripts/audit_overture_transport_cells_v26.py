#!/usr/bin/env python3
"""v2.6 Overture line length per fixed geographic plant/rail source cell.

Exactly clip real 2026-10-21.0 Overture transport line geometry to each
of the 90 fixed city cells. Measure major ROAD and selected RAIL length
(km), per-cell km/km2 to control infrastructure opportunity before
attributing observed plant similarity to rail connectivity.

This is a spatial DATA-SOURCE covariate, not road traffic or rail service.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import json
import math
from pathlib import Path

from pyproj import Geod
from shapely import wkb
from shapely.geometry import box

from audit_overture_gbif_overlap_v19 import (
    REGIONS, OVERTURE, RAIL_CLASSES, ROAD_CLASSES, open_network_connection
)
from harvest_gbif_cell_species_facets_v20 import grid_cells

GEOD=Geod(ellps="WGS84")
VALID_CLASSES={"rail":set(RAIL_CLASSES),"road":set(ROAD_CLASSES)}
TRANSPORT_TYPES=("rail","road")


def query_transport(frame):
    _,w,s,e,n=frame
    road=",".join("'"+x+"'" for x in ROAD_CLASSES)
    rail=",".join("'"+x+"'" for x in RAIL_CLASSES)
    return f"""SELECT id, subtype, "class" AS transport_class,
               ST_AsWKB(geometry) AS geometry_wkb
       FROM read_parquet('{OVERTURE}',hive_partitioning=1)
       WHERE bbox.xmin <= {e} AND bbox.xmax >= {w}
         AND bbox.ymin <= {n} AND bbox.ymax >= {s}
         AND ((subtype='rail' AND "class" IN ({rail}))
           OR (subtype='road' AND "class" IN ({road})))
       ORDER BY id"""


def geo_line_length_km(geom):
    """Length of clipped lines only, not area or point intersections."""
    if geom.is_empty:return 0.
    if geom.geom_type in ("LineString","MultiLineString"):
        d=abs(float(GEOD.geometry_length(geom)))
        if not math.isfinite(d) or d<0:
            raise ValueError("Nonfinite geographic line length")
        return d/1000
    if geom.geom_type=="GeometryCollection":
        return sum(geo_line_length_km(g) for g in geom.geoms)
    if geom.geom_type in ("Point","MultiPoint"):
        return 0.
    raise ValueError(f"Unexpected source transport geometry {geom.geom_type}")


def cell_area_km2(geometry):
    area,_=GEOD.geometry_area_perimeter(geometry)
    if not math.isfinite(area) or area==0:
        raise ValueError("Degenerate target source cell area")
    return abs(area)/1e6


def build_cell_geography(name):
    if name not in REGIONS:
        raise ValueError("Unregistered geographic country window")
    cells=grid_cells(name)
    geoms={r["cell_id"]:box(*r["bbox"]) for r in cells}
    areas={cid:cell_area_km2(geometry) for cid,geometry in geoms.items()}
    if len(cells)!=9 or len(areas)!=9 or any(a<=0 for a in areas.values()):
        raise ValueError("Invalid frozen 3x3 city-cell grid")
    return geoms,areas


def process_rows(name,source_rows):
    geometries,areas=build_cell_geography(name)
    metrics={cid:{
        "area_km2":round(area,6),
        "road_class_length_km":defaultdict(float),
        "rail_class_length_km":defaultdict(float),
        "n_intersecting_road_segments":0,
        "n_intersecting_rail_segments":0,
    } for cid,area in areas.items()}
    n_rows=0
    by_subtype=Counter()
    seen=set()
    for sid,subtype,klass,raw in source_rows:
        if not isinstance(sid,str) or not sid:
            raise ValueError("Missing original Overture segment ID")
        if sid in seen:
            raise ValueError("Duplicate segment ID found in geographic source")
        seen.add(sid)
        if subtype not in TRANSPORT_TYPES or klass not in VALID_CLASSES[subtype]:
            raise ValueError("Unexpected source transport class")
        geom=wkb.loads(bytes(raw))
        if geom.is_empty or geom.geom_type not in ("LineString","MultiLineString"):
            raise ValueError("Transport source not a valid mapped polyline")
        if not geom.is_valid:
            raise ValueError("Invalid transport polyline source")
        if any(not math.isfinite(v) for v in geom.bounds):
            raise ValueError("Invalid transport geometry bounds")
        n_rows+=1
        by_subtype[subtype]+=1
        xmin,ymin,xmax,ymax=geom.bounds
        for cid,cell in geometries.items():
            a,b,c,d=cell.bounds
            if xmax<=a or xmin>=c or ymax<=b or ymin>=d:
                continue
            piece=geom.intersection(cell)
            km=geo_line_length_km(piece)
            if km<=0:continue
            m=metrics[cid]
            key=f"{subtype}_class_length_km"
            m[key][klass]+=km
            m[f"n_intersecting_{subtype}_segments"]+=1
    for cid,m in metrics.items():
        rail=round(sum(m["rail_class_length_km"].values()),6)
        road=round(sum(m["road_class_length_km"].values()),6)
        m["road_class_length_km"]={k:round(v,6) for k,v in sorted(m["road_class_length_km"].items())}
        m["rail_class_length_km"]={k:round(v,6) for k,v in sorted(m["rail_class_length_km"].items())}
        m["road_km"]=road
        m["rail_km"]=rail
        m["road_length_density_km_per_km2"]=round(road/m["area_km2"],8)
        m["rail_length_density_km_per_km2"]=round(rail/m["area_km2"],8)
        m["area_degrees_not_equal_area"]:False
    return {
        "country":REGIONS[name][0],"n_geographic_cells":9,
        "n_source_segment_ids":n_rows,
        "source_segment_classes":dict(by_subtype),
        "cells":metrics,
        "all_road_length_km_within_target_cells":round(sum(x["road_km"] for x in metrics.values()),6),
        "all_rail_length_km_within_target_cells":round(sum(x["rail_km"] for x in metrics.values()),6),
        "source_length_clipped_to_each_cell":True,
        "not_historical_active_traffic_or_maintenance":True,
    }


def run(names, db=None):
    if not names or len(names)!=len(set(names)) or any(n not in REGIONS for n in names):
        raise ValueError("Invalid outcome-independent source region selection")
    result={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "overture_source":OVERTURE,
        "country_windows":list(names),
        "source_data_scope":"real 2026 Overture segment-length geodesic clipped per source cell",
        "regions":{},"errors":{},
    }
    owned=db is None
    if owned:
        db=open_network_connection()
        db.execute("INSTALL spatial; LOAD spatial")
    try:
        for name in names:
            try:
                rows=db.execute(query_transport(REGIONS[name])).fetchall()
                v=process_rows(name,rows)
                result["regions"][name]=v
                print(json.dumps({
                    "region":name,"country":v["country"],
                    "source_segments":v["n_source_segment_ids"],
                    "rail_km":v["all_rail_length_km_within_target_cells"],
                    "road_km":v["all_road_length_km_within_target_cells"],
                    "n_cells_with_road":sum(x["road_km"]>0 for x in v["cells"].values()),
                    "n_cells_with_rail":sum(x["rail_km"]>0 for x in v["cells"].values()),
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                result["errors"][name]=f"{type(exc).__name__}: {str(exc)[:250]}"
                print(json.dumps({"source_error":name,"error":result["errors"][name]},ensure_ascii=False),flush=True)
    finally:
        if owned:db.close()
    result["complete"]=len(result["regions"])==len(names) and not result["errors"]
    result["scientific_limits"]=[
        "Physical map kilometers, not roads actually driven or rail trains running or freight dispatches.",
        "Source segment line pieces are clipped to city study rectangles but incomplete beyond study windows.",
        "Map capture and network update years can differ from 2010-2017/2018-2025 GBIF periods.",
        "Cell line length is a land use proxy, not actual planted trackbed habitats, mowing, wind or seed deposition.",
        "No botanical effect estimate can be inferred from infrastructure length alone."
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--regions",nargs="+",choices=list(REGIONS),
                   default=["jp_sendai","de_leipzig"])
    p.add_argument("--output",default="results/audit/overture_transport_cell_lengths_v26.json")
    a=p.parse_args()
    result=run(a.regions)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":result["complete"],"regions":len(result["regions"]),
                     "errors":result["errors"]},ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("Incomplete Overture road+rail per-cell length source")


if __name__=="__main__":
    main()
