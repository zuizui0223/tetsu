"""Independent dated railway history source audit: Spain and US STB.
No ecological outcomes, no source threshold changes, no raw geometry export.
"""
import argparse
import collections
import hashlib
import io
import json
import math
import pathlib
import re
import urllib.parse
import zipfile
import shapefile
from pyproj import CRS,Transformer
from shapely.geometry import shape
from shapely.ops import transform,unary_union,linemerge
from shapely.strtree import STRtree
from audit_dated_rail_closure_sources import get,safe_source

DV="https://dataverse.csuc.cat"
PID="doi:10.34810/data917"
STB="https://www.arcgis.com/sharing/rest/content/items/59c5662600854756a7e6f18bca1a0f44"
YEARS={1985,1990,1995,*range(2000,2023)}
FILES={"IberianGauge_Lines.zip":2548,"NarrowGauge_Lines.zip":1681,"HighSpeed_lines.zip":50}
REQ={"opening","closure","reopening","line_code","length","f_id"}
CR=CRS.from_epsg(3035)
EQ=Transformer.from_crs(CR,6933,always_xy=True)

def year(x):
    try:
        y=float(x)
        return int(y) if math.isfinite(y) and y==int(y) and 1848<=y<=2023 else None
    except (ValueError,TypeError,OverflowError):
        return None

def absent(x):
    return x is None or str(x).strip() in ("","0","0.0")

def token(x):
    if absent(x):return "empty_or_zero"
    if year(x) is not None:return "valid_year"
    return "other_value"

def epochs(y):
    p=[x for x in YEARS if x<y and x-5 in YEARS]
    q=[x for x in YEARS if x>y and x+5 in YEARS]
    return (max(p)-5,max(p),min(q),min(q)+5) if p and q else None

def read_spain_zip(data,fn,expected):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        dbfs=[n for n in z.namelist() if n.lower().endswith(".dbf") and not n.startswith("__MACOSX/")]
        if len(dbfs)!=1:raise ValueError(fn+" unexpected DBF count")
        stem=dbfs[0][:-4]
        names={n.lower():n for n in z.namelist()}
        parts={}
        for ext in ("dbf","shp","shx","prj"):
            name=names.get((stem+"."+ext).lower())
            if not name or z.getinfo(name).file_size>50_000_000:raise ValueError(fn+" missing or oversized "+ext)
            parts[ext]=z.read(name)
        enc="utf-8"
        cpg=names.get((stem+".cpg").lower())
        if cpg:
            label=z.read(cpg).decode("ascii").strip().upper()
            if label in ("UTF-8","UTF8","65001"):enc="utf-8"
            elif label in ("1252","CP1252","WINDOWS-1252"):enc="cp1252"
            elif label in ("ISO-8859-1","ISO8859-1"):enc="latin-1"
            else:raise ValueError("unrecognized CPG: "+label)
        rd=shapefile.Reader(dbf=io.BytesIO(parts["dbf"]),shp=io.BytesIO(parts["shp"]),shx=io.BytesIO(parts["shx"]),encoding=enc)
        fs=[f[0] for f in rd.fields[1:]]
        print(json.dumps({"spain_layer_schema":{"file":fn,"fields":fs,"rows":len(rd),"encoding":enc,"crs":parts["prj"].decode("ascii","replace")}}),flush=True)
        if not REQ.issubset({f.lower() for f in fs}) or len(rd)!=expected:
            raise ValueError(fn+" published row/schema census mismatch")
        tf=Transformer.from_crs(CRS.from_wkt(parts["prj"].decode("ascii")),CR,always_xy=True)
        rows=[];seen=set()
        for rec in rd.iterShapeRecords():
            r={k.lower():v for k,v in rec.record.as_dict().items()}
            if r["f_id"] in seen:raise ValueError(fn+" duplicate feature ID")
            seen.add(r["f_id"])
            g=shape(rec.shape.__geo_interface__)
            if g.is_empty or not g.is_valid or g.geom_type not in ("LineString","MultiLineString"):
                raise ValueError(fn+" invalid line")
            g=transform(tf.transform,g)
            if not math.isfinite(g.length) or g.length<=0:raise ValueError("invalid projected length")
            rows.append((fn,r,g))
        return rows,{"file":fn,"fields":fs,"rows":len(rows),"encoding":enc,
                     "member_sha256":{k:hashlib.sha256(v).hexdigest() for k,v in parts.items()}}

def inspect_spain(rows):
    qa=collections.Counter();candidate=collections.defaultdict(list);competing=[];history=[]
    for fn,r,g in rows:
        a,b,c=(year(r[k]) for k in ("opening","closure","reopening"))
        qa["closure_"+token(r["closure"])]+=1
        qa["reopening_"+token(r["reopening"])]+=1
        # Conservative: only finite valid closure with valid prior opening,
        # explicitly empty/zero reopening, and a named line code. Unknown
        # values are competing, not an assumed closure. Verify README
        # sentinel conventions before ecological use.
        if a is not None and b is not None and a<=b and absent(r["reopening"]) and str(r["line_code"] or "").strip():
            candidate[(fn,str(r["line_code"]).strip(),b)].append(g)
            history.append((b,g))
            qa["terminal_candidate_rows"]+=1
        else:
            competing.append(g)
            qa["other_or_reopened_rows"]+=1
    if not competing or not candidate:raise ValueError("no viable dated candidates or competing routes")
    tree=STRtree(competing)
    ht=STRtree([g for y,g in history])
    stats=[]
    for (fn,k,y),gg in sorted(candidate.items()):
        g=unary_union(gg)
        if g.geom_type=="MultiLineString":g=linemerge(g)
        lines=list(g.geoms) if g.geom_type=="MultiLineString" else [g]
        total=clean=0.;points=[]
        for line in lines:
            if line.geom_type!="LineString":raise ValueError("unexpected merged geometry")
            L=line.length;total+=L
            n=max(1,math.ceil(L/1000))
            for i in range(n):
                p=line.interpolate((i+.5)*L/n)
                if len(tree.query(p,predicate="dwithin",distance=500)):continue
                if any(history[int(h)][0]>y for h in ht.query(p,predicate="dwithin",distance=500)):continue
                clean+=L/n;points.append(p)
        if points:
            x0,y0=EQ.transform(points[0].x,points[0].y)
            block=(math.floor(x0/200000),math.floor(y0/200000))
        else:block=None
        stats.append({"year":y,"gauge":fn,"total":total/1000,"clean":clean/1000,"block":block,"epoch":epochs(y)})
    def aggregate(x):
        good=[a for a in x if a["clean"]>=5]
        temp=[a for a in good if a["epoch"]]
        return {"route_year_groups":len(x),"terminal_union_km":round(sum(a["total"] for a in x),3),
                "screened_km_proxy":round(sum(a["clean"] for a in x),3),
                "groups_ge5km":len(good),"temporal_ge5km_groups":len(temp),
                "temporal_200km_anchor_blocks":len({a["block"] for a in temp}),
                "temporal_epoch_patterns":dict(collections.Counter(str(a["epoch"]) for a in temp))}
    return {"attribute_counts":dict(qa),"all":aggregate(stats),
            "by_era":[{"years":[a,b],**aggregate([s for s in stats if a<=s["year"]<=b])} for a,b in [(1848,1984),(1985,1990),(1991,2000),(2001,2010),(2011,2015),(2016,2023)]],
            "limits":["route-code x gauge x closure-year groups not independent closure events",
                      "blank/zero reopening only a conservative source candidate; confirm README convention",
                      "within-source passenger/rail history and 500m geometry screen, not modern Overture cross-check",
                      "point-weighted lengths approximate; no raw geography or labels exported"]}

def spain():
    url=DV+"/api/datasets/:persistentId/?persistentId="+urllib.parse.quote(PID,safe=":")
    raw,rec=get(url,6_000_000);j=json.loads(raw)
    if j.get("status")!="OK":raise ValueError("Dataverse metadata error")
    ver=j["data"]["latestVersion"]
    if ver.get("versionNumber")!=1 or ver.get("versionMinorNumber")!=0:raise ValueError("not dataset V1.0")
    files={f["dataFile"]["filename"]:f for f in ver["files"]}
    if not set(FILES)|{"Readme.txt"} <= set(files):raise ValueError("missing published files")
    print(json.dumps({"spain_dataset_metadata":{"receipt":rec,"file_names":sorted(files),"version":[1,0]}}),flush=True)
    def fetch(fn):
        f=files[fn];df=f["dataFile"]
        if f.get("restricted") or df.get("restricted"):raise ValueError("restricted source: "+fn)
        raw,rc=get(DV+"/api/access/datafile/"+str(df["id"]),12_000_000)
        checksum=df.get("checksum",{})
        alg=str(checksum.get("type","")).lower().replace("-","")
        expected=checksum.get("value")
        if not expected or alg not in ("md5","sha1","sha256","sha512") or hashlib.new(alg,raw).hexdigest().lower()!=str(expected).lower():
            raise ValueError("Dataverse file checksum mismatch or absent: "+fn)
        return raw,{**rc,"file_id":df["id"],"published_checksum":checksum}
    readme,rm=fetch("Readme.txt")
    # Store checksum and source-document presence; do not infer absent=never-reopened
    # without a separate documented semantic audit.
    text=readme.decode("utf-8-sig",errors="replace")
    print(json.dumps({"readme_schema_mentions":[s.strip()[:150] for s in text.splitlines() if re.search(r"reopen|closure|closed|cierre|reapertura|no data|missing|null",s,re.I)][:15]}),flush=True)
    rows=[];receipts={}
    for fn,n in FILES.items():
        raw,rc=fetch(fn);v,a=read_spain_zip(raw,fn,n);rows.extend(v);receipts[fn]={**rc,**a}
    return {"status":"source geometry audit","metadata":rec,"readme":rm,"files":receipts,
            "attribution":"Esteban-Oliver & Martí-Henneberg (2023), Railways of Spain GIS (1848–2023), V1, doi:10.34810/data917, CC BY 4.0",
            "audit":inspect_spain(rows)}

def stb():
    meta,rec=get(STB+"?f=json",2_000_000);m=json.loads(meta)
    if m.get("error") or m.get("access")!="public" or m.get("type")!="Web Map":
        raise ValueError("public official map metadata unavailable")
    raw,rc=get(STB+"/data?f=json",4_000_000);j=json.loads(raw)
    layers=j.get("operationalLayers",[])
    selected=[]
    for l in layers:
        title=str(l.get("title",""))
        if not re.search(r"abandon|railbank",title,re.I):continue
        url=l.get("url")
        if not url:
            selected.append({"title":title,"status":"no public service URL"})
            continue
        u=urllib.parse.urlsplit(url)
        if u.scheme!="https" or not u.hostname or not u.hostname.endswith(".arcgis.com") or not re.search(r"/FeatureServer/\d+$",u.path):
            selected.append({"title":title,"status":"unexpected service endpoint; not followed"});continue
        try:
            md,mr=get(url+"?f=json",3_000_000);d=json.loads(md)
            if "error" in d:raise ValueError("service metadata access error")
            fields=[{"name":f["name"],"alias":f.get("alias"),"type":f.get("type")} for f in d.get("fields",[])]
            count,cr=get(url+"/query?"+urllib.parse.urlencode({"where":"1=1","returnCountOnly":"true","f":"json"}),1_000_000)
            c=json.loads(count)
            if not isinstance(c.get("count"),int):raise ValueError("count error")
            selected.append({"title":title,"status":"public layer schema and count only","service":mr,"count_receipt":cr,"n_features":c["count"],"fields":fields,"geometry_type":d.get("geometryType"),"copyright":d.get("copyrightText")})
        except Exception as e:
            selected.append({"title":title,"status":"blocked","error":type(e).__name__,"reason":str(e)[:160]})
    return {"status":"map metadata audit","map":rec,"data":rc,"layers":selected,
            "limits":"Public official map; approximate applicant-supplied locations, legal abandonment/railbanking not actual traffic cessation. Do not infer dates from ambiguous fields. If underlying data restricted request from STB; no bypass."}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",default="results/railstatus/dated_international_sources.json")
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    assert epochs(1990) is None and epochs(2006)==(2000,2005,2007,2012)
    assert year("2000")==2000 and year("9999") is None and not absent("9999")
    if a.self_test:return
    result={"version":"international-dated-rail-source-audit-v1","spain":safe_source(spain),"us_stb":safe_source(stb),
            "japan_temporal_result":"20 groups, 11 blocks, 58 sites; insufficient frozen minimum 25; no rasters"}
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
    print(json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":main()
