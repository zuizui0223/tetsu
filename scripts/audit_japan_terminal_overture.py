"""Source-only independent current-map cross-check of dated Japan terminal groups."""
import argparse
import collections
import hashlib
import json
import pathlib
from shapely import wkb
from shapely.geometry import Point,box
from shapely.ops import transform
from shapely.strtree import STRtree
from pyproj import CRS,Transformer
import audit_japan_terminal_rail_geometry as j
from build_global_Overture_physical_connector_components_v361 import validate_files
from sample_entire_global_GERS_railkm_WGS84_v368b import stream_source_rails
from aggregate_world_rail_by_status_railstatus import category, INACTIVE_STATUS

def wholly_inactive(row):
    return category(row)=="inactive_normal" and any(
        isinstance(r,dict) and set(r.get("values") or ()) & INACTIVE_STATUS
        and (r.get("between") is None or list(r["between"])==[0,1])
        for r in row.get("rail_flags") or ())

def split_current(paths):
    crs=CRS.from_proj4("+proj=aeqd +lat_0=36 +lon_0=137 +datum=WGS84 +units=m +no_defs")
    proj=Transformer.from_crs("EPSG:4326",crs,always_xy=True).transform
    region=box(121,19,156,47)
    qa=collections.Counter()
    yes=[];other=[]
    for row in stream_source_rails(paths):
        qa["all_world_rail_rows"]+=1
        g=wkb.loads(bytes(row["geometry"]))
        if g.is_empty or g.geom_type!="LineString":
            raise ValueError("invalid Overture geometry")
        if not g.intersects(region) or g.bounds[2]-g.bounds[0]>180:
            continue
        clipped=g.intersection(region)
        if clipped.is_empty:
            continue
        cat="explicit_whole_inactive" if wholly_inactive(row) else "potentially_active_or_other"
        qa[cat]+=1
        q=transform(proj,clipped)
        (yes if cat=="explicit_whole_inactive" else other).append(q)
    if qa["all_world_rail_rows"]!=2054855 or not other or not yes:
        raise ValueError("incomplete full-world source pass or regional geometry")
    return yes,other,dict(qa)

def summarize(points,inactive,other):
    t0=STRtree(inactive)
    t1=STRtree(other)
    records=collections.defaultdict(lambda:collections.Counter())
    for key,y,x0,y0,w in points:
        p=Point(x0,y0)
        if len(t1.query(p,predicate="dwithin",distance=500)):
            s="competing_mapped_rail_500m"
        elif len(t0.query(p,predicate="dwithin",distance=80)):
            s="no_competing_rail_explicit_inactive_trace_80m"
        elif len(t0.query(p,predicate="dwithin",distance=500)):
            s="no_competing_rail_inactive_trace_80_to_500m"
        else:
            s="no_mapped_rail_within_500m"
        records[(key,y)][s]+=w/1000
    keys=["competing_mapped_rail_500m","no_competing_rail_explicit_inactive_trace_80m",
          "no_competing_rail_inactive_trace_80_to_500m","no_mapped_rail_within_500m"]
    def agg(items):
        sums={k:sum(v[k] for _,v in items) for k in keys}
        clean=lambda v:sum(v[k] for k in keys[1:])
        return {"components_with_any_screened_points":len(items),
                "km_by_map_class":{k:round(v,3) for k,v in sums.items()},
                "no_competing_map_rail_km":round(sum(sums[k] for k in keys[1:]),3),
                "components_ge_5km_without_competing_map_rail":sum(clean(v)>=5 for _,v in items),
                "components_ge_5km_explicit_inactive_trace":sum(v[keys[1]]>=5 for _,v in items)}
    items=list(records.items())
    return {"all":agg(items),"year_bins":[{"years":[a,b],**agg([((k,y),v) for (k,y),v in items if a<=y<=b])} for a,b in j.BINS]}

def selftest():
    from shapely.geometry import LineString
    p=[("a",1980,0,0,1000),("b",1990,2000,0,1000),("c",2010,5000,0,1000)]
    out=summarize(p,[LineString([(-100,0),(100,0)])],[LineString([(1990,0),(2010,0)])])
    v=out["all"]["km_by_map_class"]
    assert v["no_competing_rail_explicit_inactive_trace_80m"]==1 and v["competing_mapped_rail_500m"]==1 and v["no_mapped_rail_within_500m"]==1
    assert out["all"]["no_competing_map_rail_km"]==2

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    for name in ("rail-dir","census","out"):
        p.add_argument("--"+name,type=pathlib.Path)
    args=p.parse_args()
    selftest()
    if args.self_test:
        return
    rows,geoms,receipt,src=j.read_n05()
    c,attrs=j.classification(rows)
    points=[]
    original=j.audit_geometry(rows,geoms,c,points)
    if (attrs["rows"]!=2653 or len(c)!=576 or len(points)!=4534
        or original["all"]["screened_length_proxy_km"]!=4388.978
        or original["all"]["components_ge_5_screened_km"]!=199):
        raise ValueError("original run 45934992483 failed exact replay")
    if abs(sum(v[4] for v in points)/1000-4388.978)>.00051:
        raise ValueError("midpoint weights do not replay")
    paths=validate_files(args.rail_dir,json.loads(args.census.read_text()))
    inactive,other,qa=split_current(paths)
    result={"version":"japan-n05-Overture-status-crosscheck-v1",
            "exact_japan_replay":True,"n05_source_sha256":receipt["sha256"],
            "overture_census":qa,"original_geometry":original,"map_crosscheck":summarize(points,inactive,other),
            "limits":["2026 map status only; no historical traffic or freight cessation verification",
                      "no mapped trace is not physical absence; age gradient may reflect map coverage",
                      "all non-whole-inactive mapped rail is conservatively competing",
                      "midpoint weighted length, not true abandoned length; terminal groups not independent closures",
                      "source-only aggregate, no coordinates or route IDs, no ecological outcomes"]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps(result,ensure_ascii=False))
if __name__=="__main__":
    main()
