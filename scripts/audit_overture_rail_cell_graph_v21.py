#!/usr/bin/env python3
"""v2.1 mapped rail physical-connector components by 3x3 source cells.

Read actual Overture 2026-10 transport segment connectors and line
geometry. A shared connector ID defines physical map connectivity, NOT
train service, freight flows or seed transport. Spatial boundary
truncation may split otherwise connected national railway networks.
"""
from __future__ import annotations
import argparse
from collections import defaultdict,Counter
from datetime import datetime,timezone
import json
from pathlib import Path
import sys

from pyproj import Geod
from shapely import wkb

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_graph import build_rail_connector_graph,graph_components
from audit_overture_gbif_overlap_v19 import (
    REGIONS,OVERTURE,RAIL_CLASSES,open_network_connection
)
from harvest_gbif_cell_species_facets_v20 import grid_cells

GEOD=Geod(ellps="WGS84")


def rail_query(frame):
    _,w,s,e,n=frame
    return f"""SELECT id, "class" AS rail_class, connectors,
               ST_AsWKB(geometry) AS geometry_wkb
       FROM read_parquet('{OVERTURE}',hive_partitioning=1)
       WHERE subtype='rail'
         AND "class" IN ('standard_gauge','narrow_gauge','light_rail')
         AND bbox.xmin <= {e} AND bbox.xmax >= {w}
         AND bbox.ymin <= {n} AND bbox.ymax >= {s}
       ORDER BY id"""


def cell_for_point(cells,lon,lat):
    if not (-180<=lon<=180 and -90<=lat<=90):
        raise ValueError("Invalid source coordinate")
    for c in cells:
        w,s,e,n=c["bbox"]
        if w<=lon<e and s<=lat<n:
            return c["cell_id"]
    return None


def graph_component_ids(graph):
    groups=[]
    seen=set()
    for src in sorted(graph):
        if src in seen:continue
        todo=[src]
        seen.add(src)
        group=set()
        while todo:
            x=todo.pop()
            group.add(x)
            for v in graph[x]:
                if v not in seen:
                    seen.add(v);todo.append(v)
        groups.append(group)
    return sorted(groups,key=lambda g:-len(g))


def source_region(db,name):
    if name not in REGIONS:
        raise ValueError("Unknown frozen city source window")
    rows=db.execute(rail_query(REGIONS[name])).fetchall()
    if not rows:
        raise ValueError("No mapped eligible-class source rail segments")
    cells=grid_cells(name)
    segments=[]
    connector_cells=defaultdict(set)
    n_geometry_errors=0
    n_schema_errors=0
    source_error_samples=[]
    for sid,klass,connectors,geom_wkb in rows:
        try:
            if not isinstance(connectors,list) or len(connectors)<2:
                raise ValueError("Missing connectors array")
            geo=wkb.loads(bytes(geom_wkb))
            if geo.is_empty or geo.geom_type!="LineString":
                raise ValueError("Unexpected Overture segment geometry")
            length=abs(float(GEOD.geometry_length(geo)))
            if length<=0:raise ValueError("Nonpositive source geodesic length")
            cont=[]
            for c in connectors:
                if not isinstance(c,dict):raise ValueError("Unrecognized connector element")
                cid=c.get("connector_id")
                at=c.get("at")
                cont.append({"connector_id":str(cid),"at":at})
                if isinstance(at,(float,int)) and 0<=at<=1:
                    p=geo.interpolate(float(at),normalized=True)
                    found=cell_for_point(cells,p.x,p.y)
                    if found:connector_cells[str(cid)].add(found)
            segments.append({"id":str(sid),"class":klass,"subtype":"rail",
                             "length_m":length,"connectors":cont})
        except Exception as exc:
            n_geometry_errors+=1
            if len(source_error_samples)<5:
                source_error_samples.append(f"{type(exc).__name__}: {str(exc)[:100]}")
    if n_geometry_errors:
        raise ValueError(f"{n_geometry_errors}/{len(rows)} invalid rail geospatial segments; e.g. {source_error_samples}")
    graph,stats=build_rail_connector_graph(segments,allowed_classes=RAIL_CLASSES)
    parts=graph_component_ids(graph)
    mapping=defaultdict(set)
    for i,component in enumerate(parts):
        for connector in component:
            for cell in connector_cells.get(connector,()):
                mapping[cell].add(i)
    ids=[c["cell_id"] for c in cells]
    pair_counts=[]
    for i,a in enumerate(ids):
        for b in ids[i+1:]:
            conn=bool(mapping[a]&mapping[b])
            pair_counts.append({"a":a,"b":b,"within_clipped_graph_connected":conn})
    return {
        "city":name,"country":REGIONS[name][0],
        "n_eligible_source_rail_segments":len(rows),
        "n_mapped_source_rail_segments_in_connector_graph":stats["n_eligible_segments"],
        "n_connectors_in_source_graph":len(graph),
        "n_distinct_components_in_source_window":len(parts),
        "largest_component_n_connectors":max(map(len,parts),default=0),
        "n_cells_total":len(cells),
        "n_cells_with_eligible_connector":sum(bool(mapping[i]) for i in ids),
        "cell_component_cover":{
            i:{"n_source_components":len(mapping[i]),
               "connected_to_at_least_one_other_cell":any(
                   i==p["a"] or i==p["b"] for p in pair_counts if p["within_clipped_graph_connected"]
               )}
            for i in ids
        },
        "n_cell_pairs":len(pair_counts),
        "n_cell_pairs_connected_in_clipped_source_graph":sum(p["within_clipped_graph_connected"] for p in pair_counts),
        "cell_pair_connection_descriptors":pair_counts,
        "n_source_errors":n_geometry_errors,
        "inference_limits":[
            "Global Overture release includes mapped physical line geometry, not current passenger/freight operations.",
            "The extracted rail graph is geographically truncated at the study rectangle; absence of a connection is not proof of national disconnection.",
            "Connector linear references are used to interpolate approximate locations within 3x3 WGS84 city cells, not verified public boarding access or seed release.",
            "No vascular plant species, alien status, botanical observations, causal plant adaptation or rail propagule flux is tested in this graph-only source audit.",
        ]
    }


def run(names):
    if not names or any(x not in REGIONS for x in names) or len(names)!=len(set(names)):
        raise ValueError("Unexpected region frame")
    out={"run_utc":datetime.now(timezone.utc).isoformat(),
         "overture_2026_10_source":OVERTURE,
         "source_license":"ODbL / Overture Maps and OpenStreetMap contributors",
         "regions":{},"source_errors":{},
         "type":"Physical mapped source rail connectors × geographic cell components; not ecology"}
    conn=open_network_connection()
    try:
        conn.execute("INSTALL spatial; LOAD spatial;")
        for name in names:
            try:
                out["regions"][name]=source_region(conn,name)
                v=out["regions"][name]
                print(json.dumps({"city":name,"rail_segments":v["n_eligible_source_rail_segments"],
                    "graph_connectors":v["n_connectors_in_source_graph"],
                    "components":v["n_distinct_components_in_source_window"],
                    "cells_covered":v["n_cells_with_eligible_connector"],
                    "cells_pairs_connected":v["n_cell_pairs_connected_in_clipped_source_graph"]},ensure_ascii=False),flush=True)
            except Exception as exc:
                out["source_errors"][name]=f"{type(exc).__name__}: {str(exc)[:300]}"
                print(json.dumps({"source_error":name,"description":out["source_errors"][name]},ensure_ascii=False),flush=True)
    finally:
        conn.close()
    out["complete"]=not out["source_errors"] and len(out["regions"])==len(names)
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--regions",nargs="+",choices=list(REGIONS),default=["jp_sendai","de_leipzig"])
    p.add_argument("--output",default="results/audit/overture_rail_connectors_jp_de_v21.json")
    a=p.parse_args()
    res=run(a.regions)
    path=Path(a.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":res["complete"],"regions":len(res["regions"]),
                      "source_errors":res["source_errors"]},ensure_ascii=False),flush=True)
    if not res["complete"]:
        raise SystemExit("Rail connector × source-cell graph invalid; no connectivity claim")


if __name__=="__main__":
    main()
