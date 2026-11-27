"""US STB public legal-completion event cohort; no ecological outcomes."""
import argparse
import collections
import datetime as dt
import json
import math
import pathlib
import re
import urllib.parse
import numpy as np
from shapely import wkb
from shapely.geometry import LineString,MultiLineString,box
from shapely.ops import transform,unary_union,linemerge
from shapely.strtree import STRtree
from pyproj import Transformer
from audit_dated_rail_closure_sources import get
from audit_international_dated_rail_sources import epochs
from build_global_Overture_physical_connector_components_v361 import validate_files
from sample_entire_global_GERS_railkm_WGS84_v368b import stream_source_rails
from audit_japan_terminal_overture import wholly_inactive

BASE="https://services3.arcgis.com/6rJKAjBRDRSfjCzV/arcgis/rest/services/"
LAYERS={"abandoned":("Abandoned_rail_lines",530),"railbanked":("Railbanked_rail_lines",171)}
REGION=box(-125,24,-66,50)
PROJ=Transformer.from_crs(4326,5070,always_xy=True)
EQ=Transformer.from_crs(5070,6933,always_xy=True)
FIELDS="FID,Docket,Filed,Approved,Completed,State,ID,Length"
def date(x):
    s=str(x or "").strip()
    if re.fullmatch(r"(18|19|20)\d{2}",s):
        y=int(s)
        return (y,None) if y<=2026 else (None,None)
    for fmt in ("%m/%d/%Y","%m-%d-%Y","%Y-%m-%d","%m/%d/%Y %H:%M:%S","%B %d, %Y","%b %d, %Y"):
        try:
            d=dt.datetime.strptime(s,fmt).date()
            if 1800<=d.year<=2026:return d.year,d
        except ValueError:pass
    return None,None

def fetch_layer(name,expected):
    url=BASE+name+"/FeatureServer/0"
    q=urllib.parse.urlencode({"where":"1=1","outFields":FIELDS,"returnGeometry":"true","outSR":"4326","orderByFields":"FID ASC","f":"json"})
    raw,rc=get(url+"/query?"+q,30_000_000)
    j=json.loads(raw)
    feats=j.get("features")
    if (j.get("error") or j.get("exceededTransferLimit") or not isinstance(feats,list)
            or len(feats)!=expected or j.get("geometryType")!="esriGeometryPolyline"
            or j.get("spatialReference",{}).get("wkid")!=4326):
        raise ValueError(name+" unexpected feature count, truncation or geometry schema")
    rows=[];seen=set();qa=collections.Counter()
    for f in feats:
        a=f.get("attributes",{});fid=a.get("FID")
        if not isinstance(fid,int) or fid in seen:raise ValueError("invalid or duplicate FID")
        seen.add(fid)
        docket=re.sub(r"\s+","",str(a.get("Docket") or "")).upper()
        y,d=date(a.get("Completed"));ap,ad=date(a.get("Approved"));fl,fd=date(a.get("Filed"))
        qa["completion_year_valid" if y is not None else "completion_year_unknown"]+=1
        qa["completion_full_date" if d else "completion_not_full_date"]+=1
        if not docket:qa["missing_docket"]+=1
        if y is not None and ((ap is not None and ap>y) or (fl is not None and fl>y) or (ad is not None and d is not None and ad>d) or (fd is not None and d is not None and fd>d)):
            qa["date_order_conflict"]+=1
            y=None
        paths=f.get("geometry",{}).get("paths",[]) if f.get("geometry") else []
        g=None
        try:
            ll=[LineString(v) for v in paths if len(v)>=2 and all(len(p)==2 and all(math.isfinite(float(q)) for q in p) for p in v)]
            if len(ll)==len(paths) and ll:
                g=unary_union(ll)
                if g.is_empty or not g.is_valid or not REGION.covers(g):
                    g=None
        except (ValueError,TypeError):g=None
        qa["valid_conus_geometry" if g is not None else "geometry_excluded"]+=1
        rows.append((docket,y,transform(PROJ.transform,g) if g is not None else None))
    return rows,{"receipt":rc,"n_features":len(rows),"quality":dict(qa)}

def group_events(data):
    by=collections.defaultdict(list);types=collections.defaultdict(set)
    for kind,rr in data.items():
        for docket,y,g in rr:
            if docket:
                by[(kind,docket)].append((y,g));types[docket].add(kind)
    events=[];qa=collections.Counter()
    for (kind,docket),rr in sorted(by.items()):
        if len(types[docket])!=1:
            qa[kind+"_cross_pathway_docket"]+=1;continue
        yy={y for y,g in rr}
        if None in yy or len(yy)!=1:
            qa[kind+"_unknown_or_conflicting_completion_year"]+=1;continue
        if any(g is None for y,g in rr):
            qa[kind+"_incomplete_geometry"]+=1;continue
        y=next(iter(yy))
        geom=unary_union([g for yy,g in rr])
        if geom.geom_type=="MultiLineString":geom=linemerge(geom)
        if geom.geom_type not in ("LineString","MultiLineString") or geom.length<=0:
            qa[kind+"_invalid_union"]+=1;continue
        events.append((kind,docket,y,geom))
        qa[kind+"_complete_docket"]+=1
    return events,dict(qa)

def current_map(paths):
    active=[];inactive=[];qa=collections.Counter()
    for r in stream_source_rails(paths):
        qa["world_rows"]+=1
        g=wkb.loads(bytes(r["geometry"]))
        if g.is_empty or g.geom_type!="LineString":raise ValueError("invalid source line")
        if g.bounds[2]-g.bounds[0]>180 or not g.intersects(REGION):continue
        g=g.intersection(REGION)
        if g.is_empty:continue
        q=transform(PROJ.transform,g)
        if wholly_inactive(r):
            inactive.append(q);qa["explicit_whole_inactive"]+=1
        else:
            active.append(q);qa["potential_competing_rail"]+=1
    if qa["world_rows"]!=2054855 or not active or not inactive:raise ValueError("incomplete source geometry scan")
    return active,inactive,dict(qa)

def metric(events,active=None,inactive=None):
    at=STRtree(active) if active is not None else None
    it=STRtree(inactive) if inactive is not None else None
    out=[]
    for kind,key,y,g in events:
        lines=list(g.geoms) if g.geom_type=="MultiLineString" else [g]
        km=clean=trace=0.;anchor=None
        for line in lines:
            L=line.length;n=max(1,math.ceil(L/1000));km+=L/1000
            for i in range(n):
                p=line.interpolate((i+.5)*L/n)
                if at is not None and len(at.query(p,predicate="dwithin",distance=500)):continue
                clean+=L/n/1000
                if it is not None and len(it.query(p,predicate="dwithin",distance=80)):trace+=L/n/1000
                if anchor is None:
                    a,b=EQ.transform(p.x,p.y)
                    anchor=(math.floor(a/200000),math.floor(b/200000))
        out.append({"kind":kind,"year":y,"km":km,"screened_km":clean,"trace_km":trace,"anchor":anchor})
    def agg(rr):
        good=[a for a in rr if a["screened_km"]>=5]
        time=[a for a in good if epochs(a["year"])]
        return {"n_dockets":len(rr),"union_km":round(sum(a["km"] for a in rr),3),
            "screened_km_proxy":round(sum(a["screened_km"] for a in rr),3),
            "explicit_inactive80_km_proxy":round(sum(a["trace_km"] for a in rr),3),
            "n_ge5km":len(good),"n_temporal_ge5km":len(time),"n_temporal_anchor_blocks":len({a["anchor"] for a in time}),
            "temporal_completion_year_counts":dict(sorted(collections.Counter(a["year"] for a in time).items()))}
    return {kind:agg([a for a in out if a["kind"]==kind]) for kind in LAYERS}

def selftest():
    assert date("11/02/2008")[0]==2008 and date("2008")[0]==2008 and date("unknown")[0] is None
    assert date("13/01/2012")[0] is None and date("12/31/99")[0] is None
    g=LineString([(0,0),(6000,0)])
    events,qa=group_events({"abandoned":[("A",2008,g),("A",2008,g),("B",2009,g),("C",2008,g),("C",2009,g)],"railbanked":[("B",2009,g)]})
    assert len(events)==1 and events[0][3].length==6000
    assert qa["abandoned_unknown_or_conflicting_completion_year"]==1
    assert qa["abandoned_cross_pathway_docket"]==1

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true")
    for k in ("rail-dir","census","out"):p.add_argument("--"+k,type=pathlib.Path)
    a=p.parse_args();selftest()
    if a.self_test:return
    data={};receipt={}
    for kind,(name,n) in LAYERS.items():
        data[kind],receipt[kind]=fetch_layer(name,n)
    events,qa=group_events(data)
    raw=metric(events)
    print(json.dumps({"stb_source_quality":receipt,"group_quality":qa,"before_modern_screen":raw}),flush=True)
    paths=validate_files(a.rail_dir,json.loads(a.census.read_text()))
    active,inactive,mapqa=current_map(paths)
    result={"version":"stb-legal-completion-docket-census-v1","source":receipt,"group_quality":qa,
        "before_modern_screen":raw,"map_qa":mapqa,"after_modern_screen":metric(events,active,inactive),
        "limits":["Completed string is administrative completion of legal abandonment or railbanking, NOT final train date",
        "docket deduplicated within each pathway, cross-pathway/unknown/conflicting years excluded; different dockets may still overlap",
        "strict four-digit years and calendar dates only, unresolved strings not inferred; exact date semantic verification needed",
        "CONUS bounding geometry; approximate applicant alignments; EPSG5070 midpoint distances approximate",
        ">=5km sampled proxy after all non-whole-inactive rail within 500m excluded; modern map absence is not physical absence",
        "1985/1990/1995/2000–2022 same four-epoch temporal availability as Japanese audit",
        "anchor blocks descriptive, not geographic matching or independent closure count",
        "no raw coordinates/dockets/individual records exported; no ecological outcome; source STB Railroad Map Depot"]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
    print(json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":main()
