#!/usr/bin/env python3
"""v3.2: physical railway route lengths in source-frozen 50km macro grids.

Two source-only strategies (NEITHER chosen by plant outcome):
(1) Nearest explicitly mapped connector INSIDE each 50km cell centroid.
(2) Nearest connector IN THAT CELL belonging to largest graph component.

Use source identifiers and Overture 2026-10 rails with the same 75km
halo as v3.0. Compute paths by weighted Dijkstra; never merge visually
crossing rails lacking a shared connector. Source track connectivity is
not active train service, freight traffic or plant propagule movement.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime,timezone
from heapq import heappop,heappush
from itertools import combinations
import json
import math
from pathlib import Path

from pyproj import Transformer,Geod
from shapely import wkb
from shapely.geometry import Point

from audit_overture_50km_macro_rail_v30 import (
    BLOCKS,PILOT,GRID_CRS,GEOD,RAIL_CLASSES,HALO_M,
    macro_geometry,query_macro_rail,_component_map,open_network_connection
)
from tetsu.rail_graph import build_rail_connector_graph

SNAP_MAX_M=35000.0


def source_network(block,raw_rows):
    """Compute physical mapped rail graph, connector centroids and support."""
    _,polys=macro_geometry(block)
    from_geo=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    to_geo=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    segments=[];point_by_id={};cell_ids=defaultdict(set)
    disagreement=0
    if not raw_rows:raise ValueError("No mapped rail source segments in halo")
    for sid,rail_class,connectors,raw in raw_rows:
        geom=wkb.loads(bytes(raw))
        if geom.is_empty or geom.geom_type!="LineString" or not geom.is_valid:
            raise ValueError("Invalid source rail line geometry")
        length=abs(float(GEOD.geometry_length(geom)))
        if length<=0 or not math.isfinite(length):
            raise ValueError("Invalid rail source length")
        if not isinstance(connectors,list) or len(connectors)<2:
            raise ValueError("Bad source connector list")
        parsed=[]
        for item in connectors:
            if not isinstance(item,dict):raise ValueError("Non-object connector")
            cid=item.get("connector_id");at=item.get("at")
            if not isinstance(cid,str) or not cid or isinstance(at,bool) or not isinstance(at,(int,float)) or not math.isfinite(at) or not 0<=at<=1:
                raise ValueError("Malformed connector ID or at")
            lonlat=geom.interpolate(float(at),normalized=True)
            if cid in point_by_id:
                _,_,meters=GEOD.inv(point_by_id[cid][0],point_by_id[cid][1],lonlat.x,lonlat.y)
                if meters>100:disagreement+=1
            else:
                point_by_id[cid]=(lonlat.x,lonlat.y)
            x,y=from_geo.transform(lonlat.x,lonlat.y)
            pt=Point(x,y)
            for cellid,poly in polys.items():
                if poly.covers(pt):
                    cell_ids[cellid].add(cid)
                    break
            parsed.append({"connector_id":cid,"at":float(at)})
        segments.append({"id":sid,"subtype":"rail","class":rail_class,
                         "length_m":length,"connectors":parsed})
    graph,meta=build_rail_connector_graph(segments,allowed_classes=RAIL_CLASSES)
    comp,ncomp=_component_map(graph)
    sizes=defaultdict(int)
    for node,k in comp.items():sizes[k]+=1
    biggest=max(sizes,key=sizes.get)
    anchors={"nearest_in_cell":{},"nearest_largest_component_in_cell":{}}
    cell_summary={}
    for cid,poly in polys.items():
        pt=poly.centroid
        center=to_geo.transform(pt.x,pt.y)
        nodes=cell_ids[cid]
        base=[]
        for node in nodes:
            if node not in graph:continue
            coords=point_by_id[node]
            _,_,m=GEOD.inv(center[0],center[1],coords[0],coords[1])
            base.append((m,node))
        base.sort()
        options={
            "nearest_in_cell":base,
            "nearest_largest_component_in_cell":[x for x in base if comp[x[1]]==biggest]
        }
        for strategy,choices in options.items():
            if choices and choices[0][0]<=SNAP_MAX_M:
                dist,node=choices[0]
                anchors[strategy][cid]={"node":node,"snap_m":dist,"component":comp[node]}
        cell_summary[cid]={
            "n_source_connector_nodes_in_50km_cell":len(nodes),
            "n_distinct_source_components":len({comp[x] for x in nodes if x in graph}),
            "nearest_in_cell_snap_m":round(base[0][0],3) if base else None,
            "main_component_anchor_eligible":cid in anchors["nearest_largest_component_in_cell"],
        }
    return graph,comp,anchors,cell_summary,{
        "n_source_rail_segments":len(segments),
        "n_source_graph_connectors":len(graph),
        "n_distinct_source_physical_components":ncomp,
        "largest_component_n_nodes":sizes[biggest],
        "n_same_id_positions_disagreeing_over100m":disagreement,
        "same_at_source_links_epsilon":meta["n_same_at_source_links_epsilon"],
        "halo_km":HALO_M/1000,
    }


def path_m(graph,source,targets):
    if source not in graph:raise ValueError("Unmapped source connector")
    remaining=set(targets)
    if any(t not in graph for t in remaining):raise ValueError("Unknown graph target")
    q=[(0.,source)];best={source:0.};out={}
    while q and remaining:
        d,u=heappop(q)
        if d!=best.get(u):continue
        if u in remaining:
            out[u]=d;remaining.remove(u)
        for v,w in graph[u].items():
            if not isinstance(w,(int,float)) or not math.isfinite(w) or w<=0:
                raise ValueError("Malformed physical graph weight")
            nd=d+w
            if nd<best.get(v,math.inf):
                best[v]=nd;heappush(q,(nd,v))
    return out


def route_diagnostics(block,graph,comp,anchors):
    _,polys=macro_geometry(block)
    to_geo=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    centers={cid:to_geo.transform(poly.centroid.x,poly.centroid.y)
             for cid,poly in polys.items()}
    result={}
    for strategy,anchor_map in anchors.items():
        ids=sorted(polys)
        pairs=[]
        dcache={}
        for a in ids:
            aa=anchor_map.get(a)
            if aa:
                wanted=[x["node"] for b,x in anchor_map.items()
                        if x["component"]==aa["component"] and b!=a]
                if wanted:
                    dcache[a]=path_m(graph,aa["node"],wanted)
        for a,b in combinations(ids,2):
            src=anchor_map.get(a);dst=anchor_map.get(b)
            _,_,euclid_m=GEOD.inv(*centers[a],*centers[b])
            valid=bool(src and dst and src["component"]==dst["component"])
            meter=(dcache.get(a,{}).get(dst["node"]) if valid else None)
            length=(meter+src["snap_m"]+dst["snap_m"]) if meter is not None else None
            ratio=length/euclid_m if length is not None and euclid_m>0 else None
            pairs.append({
                "a":a,"b":b,
                "both_cells_have_source_rail_anchor":bool(src and dst),
                "same_clipped_physical_rail_component":valid,
                "rail_connector_path_km":round(meter/1000,6) if meter is not None else None,
                "cell_centroid_snap_augmented_rail_path_km":round(length/1000,6) if length is not None else None,
                "source_cell_centroid_geodesic_km":round(euclid_m/1000,6),
                "rail_to_geodesic_ratio":round(ratio,6) if ratio is not None else None,
            })
        observed=[p for p in pairs if p["rail_to_geodesic_ratio"] is not None]
        ratios=[p["rail_to_geodesic_ratio"] for p in observed]
        from statistics import median
        result[strategy]={
            "n_anchored_cells":len(anchor_map),
            "n_geographical_pairs":120,
            "n_pairs_reachable_in_local_source_graph":len(observed),
            "n_pairs_with_no_clipped_source_route":120-len(observed),
            "n_route_ratios_below_point95_geodesic":sum(x<.95 for x in ratios),
            "min_rail_to_geodesic_ratio":round(min(ratios),6) if ratios else None,
            "median_rail_to_geodesic_ratio":round(median(ratios),6) if ratios else None,
            "max_rail_to_geodesic_ratio":round(max(ratios),6) if ratios else None,
            "n_distinct_rounded_ratio_bins":len(set(round(x,2) for x in ratios)),
            "pair_routes":pairs,
        }
    return result


def analyze_block(block,rows):
    graph,comp,anchors,cells,meta=source_network(block,rows)
    routes=route_diagnostics(block,graph,comp,anchors)
    return {
        "country":BLOCKS[block][0],
        "source":meta,
        "source_cells":cells,
        "routes":routes,
        "not_train_service_or_seed_transfer":True,
    }


def run(blocks):
    if not blocks or any(b not in BLOCKS for b in blocks) or len(set(blocks))!=len(blocks):
        raise ValueError("Must use frozen equal-area block frame")
    data={"utc":datetime.now(timezone.utc).isoformat(),
          "source":"Overture Maps 2026-10 physical mapped RAIL connectors",
          "scope":"Continuous source-only 50km macro physical path length with fixed 75km halo",
          "regions":{},"errors":{}}
    db=open_network_connection()
    db.execute("INSTALL spatial; LOAD spatial")
    try:
        for block in blocks:
            try:
                result=analyze_block(block,db.execute(query_macro_rail(block)).fetchall())
                data["regions"][block]=result
                m=result["routes"]["nearest_in_cell"]
                main=result["routes"]["nearest_largest_component_in_cell"]
                print(json.dumps({
                    "block":block,"segments":result["source"]["n_source_rail_segments"],
                    "sources_anchored":m["n_anchored_cells"],
                    "shortest_pairs":m["n_pairs_reachable_in_local_source_graph"],
                    "median_route_ratio":m["median_rail_to_geodesic_ratio"],
                    "main_component_pairs":main["n_pairs_reachable_in_local_source_graph"],
                    "main_component_median":main["median_rail_to_geodesic_ratio"],
                    "suspicious_subgeodesic_pairs":m["n_route_ratios_below_point95_geodesic"],
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                data["errors"][block]=f"{type(exc).__name__}: {str(exc)[:250]}"
                print(json.dumps({"source_error":block,"detail":data["errors"][block]},ensure_ascii=False),flush=True)
    finally:
        db.close()
    data["complete"]=len(data["regions"])==len(blocks) and not data["errors"]
    data["limitations"]=[
        "The path is physical mapped rail lines, not operating passenger or freight routes, transfer frequencies, direction or turnover.",
        "One source anchor per 2500km2 cell at its center is an artificial regional access proxy, not botanical plant observation location.",
        "The nearest-line and nearest-main-component anchoring strategies must BOTH be shown; do not select based on plant response.",
        "A connector path may travel outside the 75km halo and be omitted from the clipped graph, biasing paths upward.",
        "Segment length times linear-reference fraction is a source geometric approximation, not a precise line-segment geodesic measure.",
        "Equal-area blocks still require land fraction, ecoregion and country-boundary QA and source GBIF event effort.",
        "No country-alien plant records or ecological effects were used in computing these network distances.",
    ]
    return data


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--blocks",nargs="+",choices=list(BLOCKS),
                   default=["jp_tohoku","us_midwest","fr_burgundy"])
    p.add_argument("--output",default="results/audit/overture_50km_continuous_rail_routes_v32.json")
    a=p.parse_args()
    out=run(a.blocks)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],
                      "n_blocks":len(out["regions"]),"errors":out["errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Incomplete physical rail source continuous network distance")
if __name__=="__main__":
    main()
