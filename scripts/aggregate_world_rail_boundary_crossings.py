#!/usr/bin/env python3
"""Source-complete actual railway geometry crossing adjacent equal-area grid boundaries.
Conservative: within-segment interior traversals, not inferred connections at intersections,
not train service. No segment IDs or original coordinates exported.
"""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
from build_global_Overture_physical_connector_components_v361 import validate_files,SOURCE_RELEASE
from sample_entire_global_GERS_railkm_WGS84_v368b import stream_source_rails,seg_m,EXPECTED_CANDIDATE_SEGMENTS,EXPECTED_CANDIDATE_KM
from aggregate_world_rail_by_status_railstatus import category
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START,X_END,Y_END

WIDTHS=(100000,200000)
CATS=("active_normal","inactive_normal","subway","tram_lightrail")
def tile(x,y,w):
    if not (X_START<=x<X_END and Y_START<=y<Y_END):return None
    return (math.floor((x-X_START)/w),math.floor((y-Y_START)/w))

def crossing_edges(points,w):
    """Walk all positive-interior intervals in source order. Do not connect
    through a corner, an exterior interval or a grid-boundary-aligned line.
    An exact crossing at an intermediate vertex of the same polyline counts."""
    prev=None;edges=set()
    for (x0,y0),(x1,y1) in zip(points,points[1:]):
        dx=x1-x0;dy=y1-y0
        if dx==0 and dy==0:continue
        if ((dx==0 and abs((x0-X_START)/w-round((x0-X_START)/w))<1e-10) or
            (dy==0 and abs((y0-Y_START)/w-round((y0-Y_START)/w))<1e-10)):
            prev=None;continue
        ts=[0.,1.]
        for a,b,d,start in ((x0,x1,dx,X_START),(y0,y1,dy,Y_START)):
            if d==0:continue
            low,high=sorted((a,b))
            for i in range(math.floor((low-start)/w)+1,math.ceil((high-start)/w)):
                t=(start+i*w-a)/d
                if 1e-12<t<1-1e-12:ts.append(t)
        ts=sorted(set(ts))
        for a,b in zip(ts,ts[1:]):
            if b-a<=1e-12:continue
            m=(a+b)/2;cur=tile(x0+m*dx,y0+m*dy,w)
            if prev is not None and cur is not None and prev!=cur:
                if abs(cur[0]-prev[0])+abs(cur[1]-prev[1])==1:
                    edges.add(tuple(sorted((prev,cur))))
            prev=cur
    return edges

def run(rail_dir,census):
    from pyproj import Transformer,Geod
    tf=Transformer.from_crs(4326,6933,always_xy=True);geod=Geod(ellps="WGS84")
    paths=validate_files(rail_dir,json.loads(Path(census).read_text()))
    edge={w:{c:Counter() for c in CATS} for w in WIDTHS};qa=Counter();length=Counter()
    for row in stream_source_rails(paths):
        qa["all_rows"]+=1;c=category(row)
        if c is None:continue
        geom,L=seg_m(row,geod);qa[c]+=1;length[c]+=L
        coords=list(geom.coords);xx,yy=zip(*[(p[0],p[1]) for p in coords])
        if max(xx)-min(xx)>180:qa[c+"_dateline_skipped"]+=1;continue
        x,y=tf.transform(xx,yy)
        if not all(math.isfinite(v) for v in (*x,*y)):qa[c+"_nonfinite_skipped"]+=1;continue
        pts=list(zip(x,y))
        for w in WIDTHS:
            if tile(min(x),min(y),w)==tile(max(x),max(y),w) and tile(min(x),min(y),w) is not None:continue
            edge[w][c].update(crossing_edges(pts,w))
        if qa["all_rows"]%250000==0:print(json.dumps({"qa":dict(qa)}),flush=True)
    if qa["all_rows"]!=2054855 or qa["active_normal"]!=EXPECTED_CANDIDATE_SEGMENTS or abs(length["active_normal"]/1000-EXPECTED_CANDIDATE_KM)>.04:
        raise ValueError("Original rail source replay mismatch")
    return {"version":"railstatus-direct-grid-boundary-v1","release":SOURCE_RELEASE,"source_verified":True,
        "qa":dict(qa),"km":{c:length[c]/1000 for c in CATS},
        "edges":{str(w):{c:[{"a":a,"b":b,"crossing_source_segments":n} for (a,b),n in sorted(v.items())] for c,v in groups.items()} for w,groups in edge.items()},
        "scope":"direct within-segment boundary crossing; source-segment counts depend on segmentation and are not train frequencies; endpoint-only crossings may be missed"}
def main():
    p=argparse.ArgumentParser()
    for n in ("rail-dir","census","out"):p.add_argument("--"+n,required=True)
    a=p.parse_args();r=run(a.rail_dir,a.census);d=Path(a.out);d.parent.mkdir(parents=True,exist_ok=True);d.write_text(json.dumps(r,allow_nan=False))
    print(json.dumps({"qa":r["qa"],"km":r["km"],"n_edges":{w:{c:len(v) for c,v in g.items()} for w,g in r["edges"].items()}}),flush=True)
if __name__=="__main__":main()
