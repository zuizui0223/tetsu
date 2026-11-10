#!/usr/bin/env python3
"""v3.0: actual Overture physical RAIL topology on 50km equal-area macro cells.

Use the exact outcome-blind 10-of-24 200-km square geographic blocks and
16 EPSG:6933 2500km2 cells from v2.9. Read 2026-10 Overture mapped rails
with 75km HALO around each block to reduce (but NOT eliminate) false
disconnectivity caused by truncating a rail graph at cell boundaries.

Report connectivity of rail-PRESENT vs rail-PRESENT geographic cell pairs,
the indispensable positive/negative-exposure support gate BEFORE any
macro plant-community effect. Not observed train operations, freight,
railway vegetation, seed transport or physical biological evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
from itertools import combinations
import io
import json
import math
from pathlib import Path
import sys

from pyproj import Transformer,Geod
from shapely import wkb
from shapely.geometry import box
from shapely.ops import transform as shapely_transform

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_graph import build_rail_connector_graph
from audit_gbif_equalarea_macroblocks_v29 import (
    BLOCKS,PILOT,GRID_CRS,EACH_CELL_METERS,BLOCK_CELLS_PER_SIDE,
    validate_frame,make_grid
)
from audit_overture_gbif_overlap_v19 import (
    OVERTURE,RAIL_CLASSES,open_network_connection
)

GEOD=Geod(ellps="WGS84")
HALO_M=75000


def macro_geometry(block):
    if block not in BLOCKS:
        raise ValueError("Not a frozen 24-block study region")
    country,lon,lat=BLOCKS[block]
    xy=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    x,y=xy.transform(lon,lat)
    half=EACH_CELL_METERS*BLOCK_CELLS_PER_SIDE/2
    region=box(x-half,y-half,x+half,y+half)
    cells={}
    for k in range(BLOCK_CELLS_PER_SIDE):
        for j in range(BLOCK_CELLS_PER_SIDE):
            cells[f"{block}_x{j}_y{k}"]=box(
                x-half+j*EACH_CELL_METERS,y-half+k*EACH_CELL_METERS,
                x-half+(j+1)*EACH_CELL_METERS,y-half+(k+1)*EACH_CELL_METERS
            )
    if len(cells)!=16 or abs(sum(p.area for p in cells.values())-region.area)>1:
        raise ValueError("Macro equal-area geometry mismatch")
    return region,cells


def bbox_with_halo(block):
    footprint,_=macro_geometry(block)
    expanded=footprint.buffer(HALO_M,cap_style=3,join_style=2)
    inv=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    lonlat=[inv.transform(x,y) for x,y in (
        (expanded.bounds[0],expanded.bounds[1]),
        (expanded.bounds[0],expanded.bounds[3]),
        (expanded.bounds[2],expanded.bounds[1]),
        (expanded.bounds[2],expanded.bounds[3])
    )]
    west=min(x for x,y in lonlat); east=max(x for x,y in lonlat)
    south=min(y for x,y in lonlat); north=max(y for x,y in lonlat)
    if not (-180<west<east<180 and -85<south<north<85):
        raise ValueError("Invalid Overture macro source query frame")
    return west,south,east,north


def query_macro_rail(block):
    west,south,east,north=bbox_with_halo(block)
    classes=",".join("'"+x+"'" for x in RAIL_CLASSES)
    return f"""SELECT id, "class" AS rail_class, connectors,
               ST_AsWKB(geometry) AS geometry_wkb
       FROM read_parquet('{OVERTURE}',hive_partitioning=1)
       WHERE subtype='rail' AND "class" IN ({classes})
         AND bbox.xmin <= {east} AND bbox.xmax >= {west}
         AND bbox.ymin <= {north} AND bbox.ymax >= {south}
       ORDER BY id"""


def _component_map(graph):
    comp={};count=0
    for first in graph:
        if first in comp: continue
        todo=[first];comp[first]=count
        while todo:
            node=todo.pop()
            for nxt in graph[node]:
                if nxt not in comp:
                    comp[nxt]=count;todo.append(nxt)
        count+=1
    return comp,count


def source_rail_macro(block,rows):
    footprint,cell_geoms=macro_geometry(block)
    transform=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    source_segments=[];cell_nodes=defaultdict(set);rail_pass=Counter()
    seen=set();bad=0
    for sid,cls,connectors,raw in rows:
        if not isinstance(sid,str) or not sid or sid in seen:
            raise ValueError("Invalid/duplicated Overture rail source segment")
        seen.add(sid)
        if cls not in RAIL_CLASSES or not isinstance(connectors,list) or len(connectors)<2:
            raise ValueError("Non-eligible rail source class or missing connectors")
        geom=wkb.loads(bytes(raw))
        if geom.is_empty or geom.geom_type!="LineString" or not geom.is_valid:
            raise ValueError("Invalid Overture real source railway geometry")
        length=abs(float(GEOD.geometry_length(geom)))
        if not math.isfinite(length) or length<=0:
            raise ValueError("Missing Overture length source")
        conn=[]
        for c in connectors:
            if not isinstance(c,dict):
                raise ValueError("Malformed connector")
            cid=c.get("connector_id");at=c.get("at")
            if not isinstance(cid,str) or not cid or isinstance(at,bool) or not isinstance(at,(int,float)) or not 0<=at<=1:
                raise ValueError("Invalid connector linear reference")
            conn.append({"connector_id":cid,"at":float(at)})
            pt=geom.interpolate(float(at),normalized=True)
            projected=transform.transform(pt.x,pt.y)
            for cellid,poly in cell_geoms.items():
                if poly.covers(box(projected[0],projected[1],projected[0],projected[1])):
                    cell_nodes[cellid].add(cid)
                    break
        source_segments.append({"id":sid,"class":cls,"subtype":"rail",
                                "length_m":length,"connectors":conn})
        projected_line=shapely_transform(transform.transform,geom)
        if not projected_line.is_empty:
            # 'passage' tests whether an actual mapped centerline intersects
            # the equal-area cell, even if no connector sits inside.
            for cellid,poly in cell_geoms.items():
                if projected_line.intersects(poly) and projected_line.intersection(poly).length>1:
                    rail_pass[cellid]+=1
    if not source_segments:
        raise ValueError("No Overture eligible rail source returned for macro region")
    graph,meta=build_rail_connector_graph(source_segments,allowed_classes=RAIL_CLASSES)
    comp,ncomp=_component_map(graph)
    cover={}
    for cell in make_grid(block):
        cid=cell["cell_id"]
        if cid not in cell_geoms:raise ValueError("GBIF and rail grids differ")
        members={comp[n] for n in cell_nodes[cid]}
        cover[cid]={
            "n_rail_source_lines_crossing_cell":rail_pass[cid],
            "n_source_rail_connectors_in_cell":len(cell_nodes[cid]),
            "n_distinct_mapped_rail_components_with_node_in_cell":len(members),
            "has_mapped_rail_passage":rail_pass[cid]>0,
            "has_rail_connector":bool(members),
            "_components":members
        }
    pairs=[]
    for a,b in combinations(sorted(cover),2):
        sa=cover[a];sb=cover[b]
        both_pass=sa["has_mapped_rail_passage"] and sb["has_mapped_rail_passage"]
        connected=bool(sa["_components"]&sb["_components"])
        pairs.append({
            "a":a,"b":b,
            "both_have_mapped_rail_passage":both_pass,
            "both_have_source_rail_connector":sa["has_rail_connector"] and sb["has_rail_connector"],
            "same_mapped_physical_rail_component_within_75km_halo":connected,
            "source_not_active_service_or_freight":True,
        })
    if len(pairs)!=120:
        raise ValueError("Missing 50km macro source cell pairs")
    exposed=[p for p in pairs if p["both_have_mapped_rail_passage"]]
    with_connectors=[p for p in exposed if p["both_have_source_rail_connector"]]
    pos=sum(p["same_mapped_physical_rail_component_within_75km_halo"] for p in with_connectors)
    neg=len(with_connectors)-pos
    for entry in cover.values():entry.pop("_components")
    return {
        "country":BLOCKS[block][0],
        "n_50km_equal_area_source_cells":16,
        "n_eligible_source_rail_segments_in_block_plus_halo":len(source_segments),
        "n_source_rail_connector_nodes":len(graph),
        "n_mapped_rail_components_in_75km_halo":ncomp,
        "n_equal_area_cells_with_rail_passage":sum(c["has_mapped_rail_passage"] for c in cover.values()),
        "n_equal_area_cells_with_connector":sum(c["has_rail_connector"] for c in cover.values()),
        "n_rail_rail_pairs_with_both_passage":len(exposed),
        "n_rail_rail_pairs_with_both_connector":len(with_connectors),
        "n_rail_rail_pairs_physically_connected":pos,
        "n_rail_rail_pairs_not_connected_within_halo":neg,
        "exposure_positivity_gate_min3_each":pos>=3 and neg>=3,
        "cell_rail_presence":cover,
        "cell_pair_rail_connectivity":pairs,
        "halo_km":HALO_M/1000,
        "no_national_full_graph_or_freight_service_inferred":True,
    }


def run(blocks):
    validate_frame()
    if not blocks or len(set(blocks))!=len(blocks) or any(b not in BLOCKS for b in blocks):
        raise ValueError("Invalid equal-area regional block selection")
    out={"run_utc":datetime.now(timezone.utc).isoformat(),
         "source":OVERTURE,
         "source_question":"Do broad 50km equal-area cell blocks have mapped physical rail-connected and unconnected rail-present cells?",
         "regions":{},"source_errors":{},
         "not_an_ecological_railway_effect":True}
    db=open_network_connection()
    db.execute("INSTALL spatial; LOAD spatial;")
    try:
        for block in blocks:
            try:
                rows=db.execute(query_macro_rail(block)).fetchall()
                o=source_rail_macro(block,rows)
                out["regions"][block]=o
                print(json.dumps({"block":block,"country":o["country"],
                    "source_rail_segments":o["n_eligible_source_rail_segments_in_block_plus_halo"],
                    "railpass_cells":o["n_equal_area_cells_with_rail_passage"],
                    "bothrail_connected":o["n_rail_rail_pairs_physically_connected"],
                    "bothrail_disconnected":o["n_rail_rail_pairs_not_connected_within_halo"],
                    "bothrail_positivity_gate":o["exposure_positivity_gate_min3_each"]},
                    ensure_ascii=False),flush=True)
            except Exception as exc:
                out["source_errors"][block]=f"{type(exc).__name__}: {str(exc)[:280]}"
                print(json.dumps({"source_error":block,"detail":out["source_errors"][block]},ensure_ascii=False),flush=True)
    finally:
        db.close()
    out["complete"]=not out["source_errors"] and len(out["regions"])==len(blocks)
    out["n_blocks_with_bothrail_connected_disconnected_source_support"]=sum(
        v["exposure_positivity_gate_min3_each"] for v in out["regions"].values())
    out["interpretation_limits"]=[
        "Only 75km halo local physical rail connections, no countrywide reachable rail paths and no train operation or freight.",
        "Outside halo connected paths may make mapped disconnected cells falsely disconnected; a negative exposure is only locally unmapped.",
        "No GRIIS plants, observer effort, roads or climate variables are used to select railway samples or infer biological effect.",
        "Cells with rail lines but no connectors are not treated as railway absent; source lines and connector node support both reported.",
        "A minimum3 connected and minimum3 disconnected rail-present pairs is only positivity diagnostic, NOT independent ecological sample size."
    ]
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--blocks",nargs="+",choices=list(BLOCKS),default=list(PILOT))
    a.add_argument("--output",default="results/audit/overture_50km_macro_rail_components_v30.json")
    p=a.parse_args()
    res=run(p.blocks)
    target=Path(p.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":res["complete"],"n_blocks":len(res["regions"]),
        "n_blocks_with_two_rail_treatments":res["n_blocks_with_bothrail_connected_disconnected_source_support"],
        "errors":res["source_errors"]},ensure_ascii=False),flush=True)
    if not res["complete"]:
        raise SystemExit("Actual regional Overture rail exposure source incomplete")


if __name__=="__main__":
    main()
