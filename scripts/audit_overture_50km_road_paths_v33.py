#!/usr/bin/env python3
"""v3.3 physical ROAD route length on frozen 50km equal-area macro blocks.

Use the SAME outcome-frozen 200-km blocks, 16-cell EPSG:6933 map and
75-km surrounding halo used for RAIL network v3.2. Overture ROAD
physical connectors, not geometric line crossings, determine routes.

Two a priori anchors in each cell, identical definitions to v3.2 rail:
 nearest source connector to area-cell centroid IN THAT CELL, and
 nearest such connector of the block's largest physical component.

Real road vehicle access, one-way restrictions, freight permissibility,
traffic and historical 2010/2018 routes are NOT observed. This is an
optimistic geometric distance proxy and an ECOL. NEGATIVE CONTROL.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
from heapq import heappop,heappush
from itertools import combinations
import json,math
from pathlib import Path
from statistics import median

from pyproj import Transformer
from shapely import wkb

from audit_overture_50km_macro_rail_v30 import (
    BLOCKS,PILOT,HALO_M,GRID_CRS,macro_geometry,bbox_with_halo,
    open_network_connection,GEOD
)
from audit_overture_gbif_overlap_v19 import OVERTURE,ROAD_CLASSES
from audit_overture_road_network_v28 import _valid_road_segment

SNAP_MAX_M=35000.0
CELL_M=50000.0
MAX_SKIPPED_INCOMPLETE_ENDPOINT_ROAD_FRACTION=0.02


def road_query(block):
    w,s,e,n=bbox_with_halo(block)
    classes=",".join("'"+x+"'" for x in ROAD_CLASSES)
    return f"""SELECT id, "class" AS road_class, connectors,
             ST_AsWKB(geometry) AS geometry_wkb
      FROM read_parquet('{OVERTURE}',hive_partitioning=1)
      WHERE subtype='road' AND "class" IN ({classes})
      AND bbox.xmin <= {e} AND bbox.xmax >= {w}
      AND bbox.ymin <= {n} AND bbox.ymax >= {s}
      ORDER BY id"""


def _components(g):
    membership={};sizes=[]
    for node in g:
        if node in membership:continue
        i=len(sizes);sizes.append(0);todo=[node];membership[node]=i
        while todo:
            curr=todo.pop();sizes[i]+=1
            for other in g[curr]:
                if other not in membership:
                    membership[other]=i;todo.append(other)
    return membership,sizes


def _dijkstra(g,start,targets):
    if start not in g or not set(targets).issubset(g):
        raise ValueError("Requested road source vertex absent from graph")
    needed=set(targets);dist={start:0.};q=[(0.,start)];done={}
    while q and needed:
        d,u=heappop(q)
        if d!=dist.get(u):continue
        if u in needed:
            done[u]=d;needed.remove(u)
        for v,weight in g[u].items():
            if not math.isfinite(weight) or weight<=0:
                raise ValueError("Bad road network edge length")
            nd=d+weight
            if nd<dist.get(v,float("inf")):
                dist[v]=nd;heappush(q,(nd,v))
    return done


def network_source(block,rows):
    footprint,polys=macro_geometry(block)
    to_proj=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    to_geo=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    x0,y0,x1,y1=footprint.bounds
    pos={}
    inside=defaultdict(set)
    graph=defaultdict(dict)
    found=set();class_counts=Counter();same_at=0;disagree=0
    n_skipped_incomplete_endpoint_rows=0
    if not rows:raise ValueError("Missing Overture ROAD source rows")
    for row in rows:
        try:
            sid,cls,geo,length,links=_valid_road_segment(row)
        except ValueError as exc:
            if str(exc)!="Road endpoint connectors missing":
                raise
            # Only exclude incomplete physical source row; do not invent
            # an unknown junction or interpret missing road as plant absence.
            n_skipped_incomplete_endpoint_rows+=1
            continue
        if sid in found:raise ValueError("Duplicate road segment source ID")
        found.add(sid);class_counts[cls]+=1
        for at,cid in links:
            point=geo.interpolate(at,normalized=True)
            lon,lat=point.x,point.y
            if cid not in pos:
                pos[cid]=(lon,lat)
            else:
                _,_,err=GEOD.inv(*pos[cid],lon,lat)
                if err>100:disagree+=1
            # O(1) equal-area cell lookup, rather than 16 shapely checks
            x,y=to_proj.transform(lon,lat)
            col=math.floor((x-x0)/CELL_M)
            rowid=math.floor((y-y0)/CELL_M)
            if 0<=col<4 and 0<=rowid<4:
                name=f"{block}_x{col}_y{rowid}"
                if name in polys:inside[name].add(cid)
        for (a,u),(b,v) in zip(links,links[1:]):
            if a==b:same_at+=1
            weight=max(.001,length*(b-a))
            prev=graph[u].get(v,float("inf"))
            if weight<prev:
                graph[u][v]=weight;graph[v][u]=weight
    if n_skipped_incomplete_endpoint_rows/len(rows)>MAX_SKIPPED_INCOMPLETE_ENDPOINT_ROAD_FRACTION:
        raise ValueError(f"Road source missing-endpoint fraction exceeds fixed 2% QA gate: {n_skipped_incomplete_endpoint_rows}/{len(rows)}")
    if not graph:raise ValueError("Source Overture road graph had no edges")
    comp,sizes=_components(graph)
    biggest=max(range(len(sizes)),key=sizes.__getitem__)
    anchors={"nearest_in_cell":{},"nearest_largest_component_in_cell":{}}
    for cid,poly in polys.items():
        center=to_geo.transform(poly.centroid.x,poly.centroid.y)
        candidates=[]
        for v in inside[cid]:
            if v not in graph:continue
            _,_,meters=GEOD.inv(*center,*pos[v])
            candidates.append((meters,v))
        candidates.sort()
        choices={
            "nearest_in_cell":candidates,
            "nearest_largest_component_in_cell":[x for x in candidates if comp[x[1]]==biggest],
        }
        for strategy,options in choices.items():
            if options and options[0][0]<=SNAP_MAX_M:
                m,v=options[0]
                anchors[strategy][cid]={"node":v,"snap_m":m,"component":comp[v]}
    return graph,anchors,{
        "country":BLOCKS[block][0],
        "n_source_road_segments_within_bbox_plus_halo":len(found),
        "n_original_source_road_segments_returned":len(rows),
        "n_excluded_source_road_segments_missing_endpoints":n_skipped_incomplete_endpoint_rows,
        "excluded_incomplete_road_fraction":round(n_skipped_incomplete_endpoint_rows/len(rows),7),
        "incomplete_road_source_NOT_interpreted_as_biological_absence":True,
        "n_source_road_connectors":len(graph),
        "n_source_graph_components":len(sizes),
        "n_largest_component_connectors":sizes[biggest],
        "road_source_classes":dict(class_counts),
        "n_road_equal_at_connector_links_epsilon":same_at,
        "n_same_connector_positions_disagree_gt100m":disagree,
        "halo_km":HALO_M/1000.,
    }


def route_summary(block,graph,anchors):
    _,polys=macro_geometry(block)
    to_geo=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    centers={cid:to_geo.transform(g.centroid.x,g.centroid.y) for cid,g in polys.items()}
    result={}
    for strat,selected in anchors.items():
        paths={}
        for cid,a in selected.items():
            target=[v["node"] for name,v in selected.items() if name!=cid and v["component"]==a["component"]]
            if target:paths[cid]=_dijkstra(graph,a["node"],target)
        pairs=[]
        ratios=[]
        for x,y in combinations(sorted(centers),2):
            p=selected.get(x);q=selected.get(y)
            _,_,straight=GEOD.inv(*centers[x],*centers[y])
            reachable=bool(p and q and p["component"]==q["component"])
            line=(paths.get(x,{}).get(q["node"]) if reachable else None)
            total=(line+p["snap_m"]+q["snap_m"] if line is not None else None)
            ratio=(total/straight if total is not None and straight>0 else None)
            if ratio is not None:ratios.append(ratio)
            pairs.append({
                "a":x,"b":y,
                "both_have_source_road_anchor":bool(p and q),
                "road_physically_connected_in_halo":reachable and line is not None,
                "road_connector_path_km":round(line/1000,6) if line is not None else None,
                "centroid_snap_augmented_road_path_km":round(total/1000,6) if total is not None else None,
                "road_to_geodesic_ratio":round(ratio,6) if ratio is not None else None,
            })
        result[strat]={
            "n_anchored_cells":len(selected),
            "n_road_route_pairs":len(ratios),
            "n_road_route_circuity_below_point95":sum(t<.95 for t in ratios),
            "median_road_to_geodesic_ratio":round(median(ratios),6) if ratios else None,
            "min_road_to_geodesic_ratio":round(min(ratios),6) if ratios else None,
            "max_road_to_geodesic_ratio":round(max(ratios),6) if ratios else None,
            "pair_paths":pairs
        }
    return result


def run(blocks):
    if not blocks or len(set(blocks))!=len(blocks) or any(b not in BLOCKS for b in blocks):
        raise ValueError("Non-frozen regional sampling frame")
    result={"utc":datetime.now(timezone.utc).isoformat(),
            "source":"2026-10 Overture mapped physical roads",
            "study":"Road network ROUTE length on 50km equal-area cells (not traffic)",
            "regions":{},"errors":{}}
    db=open_network_connection()
    db.execute("INSTALL spatial; LOAD spatial")
    try:
        for name in blocks:
            try:
                rows=db.execute(road_query(name)).fetchall()
                g,anchors,meta=network_source(name,rows)
                paths=route_summary(name,g,anchors)
                result["regions"][name]={"source":meta,"routes":paths}
                v=paths["nearest_in_cell"]
                m=paths["nearest_largest_component_in_cell"]
                print(json.dumps({
                    "block":name,"road_segments":meta["n_source_road_segments_within_bbox_plus_halo"],
                    "roads_anchors":v["n_anchored_cells"],
                    "road_pairs":v["n_road_route_pairs"],
                    "median_road_circuity":v["median_road_to_geodesic_ratio"],
                    "main_component_pairs":m["n_road_route_pairs"],
                    "main_component_median_circuity":m["median_road_to_geodesic_ratio"],
                    "suspect_source_routes":v["n_road_route_circuity_below_point95"],
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                result["errors"][name]=f"{type(exc).__name__}: {str(exc)[:250]}"
                print(json.dumps({"source_error":name,"detail":result["errors"][name]},ensure_ascii=False),flush=True)
    finally:
        db.close()
    result["complete"]=len(result["regions"])==len(blocks) and not result["errors"]
    result["limitations"]=[
        "Physical roads with no turn restrictions, one-way access, freight vehicle limits, actual traffic, speeds or ferry transfers.",
        "75km around block is an incomplete regional graph; paths may leave it and return, creating artificially long routes.",
        "Road cell-centroid anchor is not an actual location of a GBIF botanical occurrence.",
        "Source only, no use of alien-plant response data for regional selection or route anchoring.",
        "Road route length times normalized source at fractions approximates rather than exactly clips geodesic polylines.",
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--blocks",nargs="+",choices=list(BLOCKS),
                   default=["us_midwest","fr_burgundy"])
    p.add_argument("--output",default="results/audit/overture_50km_macro_road_paths_v33.json")
    a=p.parse_args()
    d=run(a.blocks)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(d,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":d["complete"],"regions":len(d["regions"]),
                      "errors":d["errors"]},ensure_ascii=False),flush=True)
    if not d["complete"]:
        raise SystemExit("Macro road network source not complete")
if __name__=="__main__":
    main()
