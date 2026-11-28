"""Exact same-status GERS components supporting observed boundary crossings.

For each crossing, reach = maximum over components of the smaller
component-assigned centerline length on its two sides (km).
This is mapped network continuity, NOT verified traffic or a shortest path.
"""
import argparse, collections, json, math
from pathlib import Path
from pyproj import Transformer,Geod
from build_global_Overture_physical_connector_components_v361 import (
    Disjoint,read_rows,original_connectors,validate_files,SOURCE_RELEASE)
from aggregate_world_rail_by_status_railstatus import category
from aggregate_world_rail_boundary_crossings import crossing_edges
from aggregate_entire_world_GERS_rail_100km_200km_grid_v370 import clip_source_line_to_tiles
from sample_entire_global_GERS_railkm_WGS84_v368b import seg_m

CATS=("active_normal","inactive_normal")
W=200000
CELL_COVER=.8
CROSS_COVER=.9
def key(z):
    return tuple(z["a"]),tuple(z["b"])

def reach(comp,bycell,a,b):
    return max((min(bycell.get((a,r),0.),bycell.get((b,r),0.)) for r in comp),default=0.)

def selftest():
    assert reach({"x","y"},{((0,0),"x"):25,((1,0),"x"):12,
           ((0,0),"y"):8,((1,0),"y"):50},(0,0),(1,0))==12
    assert reach(set(),{},(0,0),(1,0))==0

def run(paths,edge,status,keys):
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or edge["release"]!=SOURCE_RELEASE or not edge["source_verified"]:
        raise ValueError("unverified crossing source")
    if status["source_Overture_release"]!=SOURCE_RELEASE:
        raise ValueError("unverified status source")
    cells={k for ab in keys for k in ab}
    old={c:{key(z):z["crossing_source_segments"] for z in edge["edges"][str(W)][c] if key(z) in keys} for c in CATS}
    known={c:{(z["source_equal_area_grid_ix"],z["source_equal_area_grid_iy"]):z["km"] for z in status["grids_km_by_category_and_width"][c][str(W)] if (z["source_equal_area_grid_ix"],z["source_equal_area_grid_iy"]) in cells} for c in CATS}
    uf={c:Disjoint() for c in CATS};qa=collections.Counter()
    for r in read_rows(paths,("class","rail_flags","connectors")):
        qa["rows"]+=1
        c=category(r)
        if c not in CATS:continue
        qa[c+"_segments"]+=1
        nn=original_connectors(r)
        if len(set(nn))<2:
            qa[c+"_unplaced_segments"]+=1;continue
        for a,b in zip(nn,nn[1:]):
            if a!=b:uf[c].union(a,b)
    if qa["rows"]!=2054855 or any(qa[c+"_segments"]!=status["segments"][c] for c in CATS):
        raise ValueError("source segment census mismatch")
    print(json.dumps({"component_graph_qa":dict(qa),"components":{c:uf[c].components for c in CATS}}),flush=True)
    fwd=Transformer.from_crs(4326,6933,always_xy=True);inv=Transformer.from_crs(6933,4326,always_xy=True);geod=Geod(ellps="WGS84")
    counts={c:collections.Counter() for c in CATS}
    placed_cross={c:collections.Counter() for c in CATS}
    crossing_components={c:collections.defaultdict(set) for c in CATS}
    bycell={c:collections.defaultdict(float) for c in CATS}
    totalcell={c:collections.defaultdict(float) for c in CATS}
    placedcell={c:collections.defaultdict(float) for c in CATS}
    km={c:0. for c in CATS}
    n=0
    for r in read_rows(paths,("class","rail_flags","connectors","geometry")):
        n+=1
        c=category(r)
        if c not in CATS:continue
        geom,L=seg_m(r,geod);km[c]+=L/1000
        nn=original_connectors(r)
        root=uf[c].find(nn[0]) if len(set(nn))>=2 else None
        xx,yy=zip(*[(p[0],p[1]) for p in geom.coords])
        if max(xx)-min(xx)<=180:
            x,y=fwd.transform(xx,yy)
            if all(math.isfinite(v) for v in (*x,*y)):
                for k in crossing_edges(list(zip(x,y)),W)&keys:
                    counts[c][k]+=1
                    if root is not None:
                        placed_cross[c][k]+=1;crossing_components[c][k].add(root)
        try:
            a,_,_=clip_source_line_to_tiles(geom,L,fwd,inv,geod)
        except ValueError:
            if c=="active_normal":raise
            qa[c+"_clip_error"]+=1;continue
        for (i,j),m in a.items():
            k=(i//2,j//2)
            if k not in cells:continue
            totalcell[c][k]+=m/1000
            if root is not None:
                placedcell[c][k]+=m/1000
                bycell[c][k,root]+=m/1000
        if n%500000==0:print(json.dumps({"component_geometry_rows":n}),flush=True)
    if n!=2054855 or any(abs(km[c]-status["total_km"][c])>.01 for c in CATS):
        raise ValueError("full length replay mismatch")
    for c in CATS:
        if dict(counts[c])!=old[c]:raise ValueError("exact original edge replay failed "+c)
        if any(abs(totalcell[c].get(k,0.)-known[c].get(k,0.))>.01 for k in cells):
            raise ValueError("exact original cell lengths differ "+c)
    edges=[];reject=collections.Counter()
    for a,b in sorted(keys):
        z={"a":a,"b":b};ok=True
        for c in CATS:
            k=(a,b);total=counts[c][k]
            # A zero mapped crossing is a defined zero reach. Positive
            # crossings require sufficient connector/length placement.
            if total==0:
                z[c]={"crossings":0,"reach_km":0.}
                continue
            q=placed_cross[c][k]/total
            cv=[placedcell[c][v]/totalcell[c][v] if totalcell[c][v]>0 else 0. for v in (a,b)]
            if q<CROSS_COVER or min(cv)<CELL_COVER:
                reject[c+"_topology_support"]+=1;ok=False;continue
            m=reach(crossing_components[c][k],bycell[c],a,b)
            if m<=0 or m>min(totalcell[c][a],totalcell[c][b])+.001:
                reject[c+"_length_support"]+=1;ok=False;continue
            z[c]={"crossings":total,"reach_km":m,"crossing_placed_fraction":q,"cell_placed_fraction":cv}
        if ok:edges.append(z)
    return {"version":"railstatus-component-reach-v1","source_release":SOURCE_RELEASE,
            "exact_crossing_and_grid_km_replay":True,"n_candidate_edges":len(keys),
            "n_supported_edges":len(edges),"qa":dict(qa),"excluded":dict(reject),"edges":edges,
            "limits":["same-status source connector components, not named railway routes or verified traffic",
             "reach max over components that have a direct within-segment crossing, min of their mapped length in each adjacent cell",
             "component may connect outside cell; not distance along shortest path or guaranteed local continuous habitat",
             "positive crossing requires >=90% of crossing segments and >=80% of both endpoint cell km graph-placed",
             "zero means no mapped direct crossing, not no physical line",
             "status labels same original Overture mapping, no raw connector IDs or coordinates exported"]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true")
    for k in ("rail-dir","census","edges","status","historical","out"):p.add_argument("--"+k,type=Path)
    a=p.parse_args();selftest()
    if a.self_test:return
    read=lambda p:json.loads(p.read_text())
    hist=read(a.historical)
    if hist["version"]!="railstatus-avonet-heldback-2018-22-hwi-pairs-v1" or not hist["exact_2025_full_and_joint_guild_replay"]:
        raise ValueError("wrong frozen HWI input")
    keys={key(z) for z in hist["pairs"]["2025"]}
    if len(keys)!=1645:raise ValueError("wrong frozen geographic frame")
    result=run(validate_files(a.rail_dir,read(a.census)),read(a.edges),read(a.status),keys)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
    print(json.dumps({k:v for k,v in result.items() if k!="edges"},allow_nan=False),flush=True)
if __name__=="__main__":main()
