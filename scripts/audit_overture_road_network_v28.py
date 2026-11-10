#!/usr/bin/env python3
"""v2.8: ROAD graph route-distance source from Overture explicit connectors.

For the same fixed 3x3 city windows used by GBIF speciesKey records,
measure 2026 Overture ROAD network *physical* connectedness and approximate
shortest graph length. Unlike prior road-km density, this distinguishes
topological path length from Euclidean distance.

NO roads are inferred at crossings without shared connector_id. No
turn/access restrictions, motor-vehicle eligibility, traffic speed or
freight journey is modeled; paths are optimistic physical connectivity,
not legal or real observed vehicle journeys. Source geography is clipped
at the city window: a disconnected or long path can leave that window.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from heapq import heappop, heappush
from itertools import combinations
import json
import math
from pathlib import Path

from pyproj import Geod
from shapely import wkb

from audit_overture_gbif_overlap_v19 import (
    REGIONS, OVERTURE, ROAD_CLASSES, open_network_connection
)
from harvest_gbif_cell_species_facets_v20 import grid_cells

GEOD=Geod(ellps="WGS84")
SOURCE_SNAP_MAX_M=5000.0
LENGTH_NEAR_ENDPOINT_TOLERANCE=1e-7


def query_road_source(frame):
    _,west,south,east,north=frame
    classes=",".join("'"+x+"'" for x in ROAD_CLASSES)
    return f"""SELECT id, "class" AS road_class, connectors,
                  ST_AsWKB(geometry) AS geometry_wkb
        FROM read_parquet('{OVERTURE}',hive_partitioning=1)
        WHERE subtype='road'
          AND "class" IN ({classes})
          AND bbox.xmin <= {east} AND bbox.xmax >= {west}
          AND bbox.ymin <= {north} AND bbox.ymax >= {south}
        ORDER BY id"""


def _valid_road_segment(row):
    sid, road_class, connectors, raw=row
    if not isinstance(sid,str) or not sid:
        raise ValueError("Missing Overture road segment ID")
    if road_class not in ROAD_CLASSES:
        raise ValueError("Unapproved road class")
    if not isinstance(connectors,list) or len(connectors)<2:
        raise ValueError("Overture road segment must have endpoint connectors")
    geo=wkb.loads(bytes(raw))
    if geo.is_empty or geo.geom_type!="LineString" or not geo.is_valid:
        raise ValueError("Invalid Overture road source line geometry")
    length_m=abs(float(GEOD.geometry_length(geo)))
    if not math.isfinite(length_m) or length_m<=0:
        raise ValueError("Zero/nonfinite source road length")
    links=[]
    for val in connectors:
        if not isinstance(val,dict):raise ValueError("Bad connector map")
        cid=val.get("connector_id")
        at=val.get("at")
        if (not isinstance(cid,str) or not cid or isinstance(at,bool) or
            not isinstance(at,(int,float)) or not math.isfinite(at) or not 0<=at<=1):
            raise ValueError("Bad connector ID or linear reference")
        links.append((float(at),cid))
    links.sort(key=lambda x:(x[0],x[1]))
    if links[0][0]>LENGTH_NEAR_ENDPOINT_TOLERANCE or 1-links[-1][0]>LENGTH_NEAR_ENDPOINT_TOLERANCE:
        raise ValueError("Road endpoint connectors missing")
    if len({cid for at,cid in links})!=len(links):
        raise ValueError("Repeated connector ID within road segment")
    if any(a[0]>=b[0] for a,b in zip(links,links[1:])):
        raise ValueError("Duplicate or nonordered connector fractions")
    return sid,road_class,geo,length_m,links


def graph_from_source_rows(rows):
    """Weighted physical connector graph with explicit source identifiers.

    The geometry-interpolated 'at' coordinates are only approximate connector
    location proxies. Summed edge weights use source geodesic segment length
    times linear-reference fraction; precise routing would segment geometry.
    """
    adj=defaultdict(dict)
    node_positions={}
    seen=set();counts=Counter()
    n_position_disagreement=0
    for row in rows:
        sid,road_class,geo,length_m,links=_valid_road_segment(row)
        if sid in seen: raise ValueError("Duplicate Overture road segment ID")
        seen.add(sid)
        counts[road_class]+=1
        for at,cid in links:
            pt=geo.interpolate(at,normalized=True)
            lon,lat=pt.x,pt.y
            if cid not in node_positions:
                node_positions[cid]=(lon,lat)
            else:
                orig=node_positions[cid]
                _,_,delta=GEOD.inv(orig[0],orig[1],lon,lat)
                if delta>100.:
                    n_position_disagreement+=1
        for (at1,a),(at2,b) in zip(links,links[1:]):
            weight=length_m*(at2-at1)
            if weight<=0: raise ValueError("Degenerate road edge")
            old=adj[a].get(b)
            if old is None or weight<old:
                adj[a][b]=weight
                adj[b][a]=weight
    if len(seen)==0 or len(adj)==0:
        raise ValueError("Empty real road source graph")
    return dict(adj),node_positions,{
        "n_eligible_source_road_segments":len(seen),
        "source_class_counts":dict(sorted(counts.items())),
        "n_source_connector_ids":len(node_positions),
        "n_position_disagreements_gt100m":n_position_disagreement,
        "physical_only_not_vehicle_route":True,
    }


def graph_components(adj):
    component_of={}
    sizes=Counter()
    for node in adj:
        if node in component_of:continue
        idx=len(sizes)
        pending=[node]
        component_of[node]=idx
        count=0
        while pending:
            current=pending.pop()
            count+=1
            for nxt in adj[current]:
                if nxt not in component_of:
                    component_of[nxt]=idx
                    pending.append(nxt)
        sizes[idx]=count
    return component_of,sizes


def haversine_m(a,b):
    lon1,lat1=a;lon2,lat2=b
    y1,y2=math.radians(lat1),math.radians(lat2)
    a_lat=(y2-y1)/2
    a_lon=math.radians(lon2-lon1)/2
    h=math.sin(a_lat)**2+math.cos(y1)*math.cos(y2)*math.sin(a_lon)**2
    return 2*6371008.8*math.asin(min(1.,math.sqrt(max(0.,h))))


def _nearest_connector(anchor,positions):
    if not positions:
        raise ValueError("Cannot snap cell to empty road graph")
    lon,lat=anchor
    cos_lat=math.cos(math.radians(lat))
    nearest=None;best=float("inf")
    for cid,(x,y) in positions.items():
        # Quick local equirectangular approximation used ONLY to rank.
        d=((x-lon)*cos_lat)**2+(y-lat)**2
        if d<best:
            nearest=cid;best=d
    true_m=haversine_m(anchor,positions[nearest])
    return nearest,round(true_m,3)


def shortest_paths_to_targets(adj,source,targets):
    """Positive-weight Dijkstra with early finish for nine cell centroids."""
    if source not in adj:raise ValueError("Source connector outside road graph")
    destinations=set(targets)
    if not destinations.issubset(adj):
        raise ValueError("Destination connector missing")
    out={}
    pq=[(0.,source)];best={source:0.}
    while pq and destinations:
        dist,node=heappop(pq)
        if dist!=best.get(node):continue
        if node in destinations:
            out[node]=dist
            destinations.remove(node)
        for dest,w in adj[node].items():
            if not isinstance(w,(float,int)) or not math.isfinite(w) or w<=0:
                raise ValueError("Invalid road edge")
            alt=dist+w
            if alt<best.get(dest,float("inf")):
                best[dest]=alt;heappush(pq,(alt,dest))
    return out


def network_for_city(name,rows):
    if name not in REGIONS:raise ValueError("Unknown fixed macro source city")
    graph,pos,meta=graph_from_source_rows(rows)
    comp,sizes=graph_components(graph)
    grid=grid_cells(name)
    snap={}
    for cell in grid:
        w,s,e,n=cell["bbox"]
        mid=((w+e)/2,(s+n)/2)
        connector,snap_dist=_nearest_connector(mid,pos)
        snap[cell["cell_id"]]={
            "road_centroid_proxy_node_id":connector,
            "centroid_snap_m":snap_dist,
            "component_id":comp[connector],
            "nearest_connector_within_5km":snap_dist<=SOURCE_SNAP_MAX_M
        }
    # Do not export actual connector IDs: they do not change the conclusion
    # and source geometry/ids may be subject to ODbL share-alike.
    labels=sorted(snap)
    pairs=[]
    for i,ci in enumerate(labels):
        si=snap[ci]
        destination_nodes=[
            snap[cj]["road_centroid_proxy_node_id"] for cj in labels[i+1:]
            if snap[cj]["component_id"]==si["component_id"]
            and snap[cj]["nearest_connector_within_5km"]
            and si["nearest_connector_within_5km"]
        ]
        dist=shortest_paths_to_targets(graph,si["road_centroid_proxy_node_id"],destination_nodes)
        box_a=next(c["bbox"] for c in grid if c["cell_id"]==ci)
        x=(box_a[0]+box_a[2])/2;y=(box_a[1]+box_a[3])/2
        for cj in labels[i+1:]:
            sj=snap[cj]
            box_b=next(c["bbox"] for c in grid if c["cell_id"]==cj)
            v=(box_b[0]+box_b[2])/2,w=(box_b[1]+box_b[3])/2
            straight=haversine_m((x,y),(v,w))
            same=si["component_id"]==sj["component_id"]
            eligible=si["nearest_connector_within_5km"] and sj["nearest_connector_within_5km"]
            path=dist.get(sj["road_centroid_proxy_node_id"]) if same and eligible else None
            # Actual centers sit away from nearest connector, so include
            # straight-line snap distances, an optimistic bound.
            total=(path+si["centroid_snap_m"]+sj["centroid_snap_m"] if path is not None else None)
            circuity=(total/straight if total is not None and straight>0 else None)
            pairs.append({
                "a":ci,"b":cj,
                "within_clipped_road_graph_physically_connected":same and path is not None,
                "source_centroids_snapped_within_5km":eligible,
                "centroid_straight_km":round(straight/1000,6),
                "optimistic_road_connector_path_km":round(total/1000,6) if total is not None else None,
                "optimistic_road_to_euclidean_ratio":round(circuity,6) if circuity is not None else None,
                "road_access_is_not_legal_motor_freight_routes":True,
            })
    if len(pairs)!=36:
        raise ValueError("Expected 9 choose 2 source cell pairs")
    counts=Counter(s["component_id"] for s in snap.values())
    return {
        "country":REGIONS[name][0],**meta,
        "n_distinct_clipped_road_graph_components":len(sizes),
        "largest_component_n_connectors":max(sizes.values()),
        "n_source_cells":9,
        "n_source_cells_snapped_within_5km":sum(s["nearest_connector_within_5km"] for s in snap.values()),
        "n_source_cell_pairs_road_connected":sum(p["within_clipped_road_graph_physically_connected"] for p in pairs),
        "n_source_cell_pairs_with_road_route_km":sum(p["optimistic_road_connector_path_km"] is not None for p in pairs),
        "road_cell_component_sizes":sorted(counts.values(),reverse=True),
        "source_cell_snap_offsets_m":{
            cid:v["centroid_snap_m"] for cid,v in snap.items()
        },
        "cell_pair_road_paths":pairs,
        "window_clipped_unknown_outside_links":True,
        "access_turn_restrictions_not_modeled":True,
    }


def run(names):
    if not names or any(n not in REGIONS for n in names) or len(set(names))!=len(names):
        raise ValueError("Unknown outcome-independent region subset")
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source":OVERTURE,
        "kind":"Overture mapped physical road connector shortest-path PROXY, not traffic or legal route",
        "regions":{},"source_errors":{}
    }
    db=open_network_connection()
    db.execute("INSTALL spatial; LOAD spatial;")
    try:
        for name in names:
            try:
                rows=db.execute(query_road_source(REGIONS[name])).fetchall()
                region=network_for_city(name,rows)
                out["regions"][name]=region
                print(json.dumps({
                    "city":name,"road_segments":region["n_eligible_source_road_segments"],
                    "road_connector_ids":region["n_source_connector_ids"],
                    "graph_components":region["n_distinct_clipped_road_graph_components"],
                    "road_connected_cell_pairs":region["n_source_cell_pairs_road_connected"],
                    "nine_cells_snap_within_5km":region["n_source_cells_snapped_within_5km"],
                    "geodesic_position_disagreements_gt100m":region["n_position_disagreements_gt100m"]
                },ensure_ascii=False),flush=True)
            except Exception as e:
                out["source_errors"][name]=f"{type(e).__name__}: {str(e)[:320]}"
                print(json.dumps({"source_error":name,"detail":out["source_errors"][name]},ensure_ascii=False),flush=True)
    finally:
        db.close()
    out["complete"]=not out["source_errors"] and len(out["regions"])==len(names)
    out["interpretation_warnings"]=[
        "Road graph connectivity is purely physical; no turn restrictions, permitted vehicles or traffic volumes modeled.",
        "Road network cropped at each city frame, so disconnection/detours can reflect source boundary rather than true regional access.",
        "Closest connector to rectangle center is a geometric proxy for botanical cell access, not botanical observation coordinates.",
        "Linear reference-weighted total segment lengths approximate road routes; spatial and restriction QA needed for true freight travel.",
        "A city where all pairs are physically road connected supplies no binary road-connected controls.",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--regions",nargs="+",choices=list(REGIONS),default=["jp_sendai","gb_york"])
    p.add_argument("--output",default="results/audit/overture_road_connector_paths_v28.json")
    a=p.parse_args()
    out=run(a.regions)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"n_regions":len(out["regions"]),
                      "errors":out["source_errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Road-graph source incomplete; do not claim a road-adjusted ecological result")


if __name__=="__main__":
    main()
