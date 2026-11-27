"""Reconcile published paper counts against checksum-verified Spain V1 archives.
Source inventory only, not a relaxation of the earlier analysis gate.
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
from shapely.geometry import shape
from shapely.ops import transform
from pyproj import CRS,Transformer
import audit_international_dated_rail_sources as s
from audit_dated_rail_closure_sources import get,safe_source

def val(x):
    return str(x if x is not None else "").strip()

def audit_zip(data,fn,paper,external_cpg=None):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        dbfs=[n for n in z.namelist() if n.lower().endswith(".dbf") and not n.startswith("__MACOSX/")]
        if len(dbfs)!=1:raise ValueError("nonunique DBF")
        stem=dbfs[0][:-4];names={n.lower():n for n in z.namelist()}
        parts={e:z.read(names[(stem+"."+e).lower()]) for e in ("dbf","shp","shx","prj")}
        internal=names.get((stem+".cpg").lower())
        labels=[b.decode("ascii").strip().upper() for b in ([z.read(internal)] if internal else [])+([external_cpg] if external_cpg else [])]
        def encoding(x):
            if x in ("UTF-8","UTF8","65001"):return "utf-8"
            if x in ("1252","CP1252","WINDOWS-1252"):return "cp1252"
            if x in ("ISO-8859-1","ISO8859-1"):return "latin-1"
            raise ValueError("unknown CPG "+x)
        encs={encoding(x) for x in labels}
        if len(encs)>1:raise ValueError("conflicting CPG")
        enc=next(iter(encs)) if encs else "utf-8"
        r=shapefile.Reader(dbf=io.BytesIO(parts["dbf"]),shp=io.BytesIO(parts["shp"]),shx=io.BytesIO(parts["shx"]),encoding=enc)
        fields=[f[0].lower() for f in r.fields[1:]]
        if not s.REQ.issubset(fields):raise ValueError("missing required attributes")
        n=len(r)
        fids=collections.Counter();linecodes=set();sectionkeys=set();events=set();gh=collections.Counter();rowhash=collections.Counter()
        status=collections.Counter();date_tokens=collections.defaultdict(collections.Counter)
        closure_years=collections.Counter();valid_rows=[];geom_total=0.
        tf=Transformer.from_crs(CRS.from_wkt(parts["prj"].decode("ascii")),s.CR,always_xy=True)
        for rec in r.iterShapeRecords():
            a={k.lower():v for k,v in rec.record.as_dict().items()}
            fid=val(a["f_id"]);code=val(a["line_code"]);fids[fid]+=1
            if code:linecodes.add(code)
            sectionkeys.add((code,val(a.get("section"))))
            y0,y1,y2=[s.year(a[k]) for k in ("opening","closure","reopening")]
            for k,y in zip(("opening","closure","reopening"),(y0,y1,y2)):
                date_tokens[k]["valid_year" if y is not None else "empty" if not val(a[k]) else "zero" if s.absent(a[k]) else "other_nonempty"]+=1
            if y0 is None:status["invalid_opening"]+=1
            if y0 is not None and y1 is not None and y1<y0:status["closure_before_opening"]+=1
            if y2 is not None and (y1 is None or y2<y1):status["reopening_without_valid_prior_closure"]+=1
            if not code:status["missing_line_code"]+=1
            if y1 is not None:closure_years[str((y1//10)*10)]+=1
            if y0 is not None and y1 is not None and y0<=y1 and s.absent(a["reopening"]) and code:
                events.add((code,y1));status["provisional_terminal_rows"]+=1
            try:
                g=shape(rec.shape.__geo_interface__)
                if g.is_empty or not g.is_valid or g.geom_type not in ("LineString","MultiLineString"):
                    raise ValueError("invalid line")
                h=hashlib.sha256(g.normalize().wkb).hexdigest();gh[h]+=1
                semantic=tuple((k,val(a[k])) for k in sorted(a) if k!="f_id")
                rh=hashlib.sha256(repr((h,semantic)).encode()).hexdigest();rowhash[rh]+=1
                g=transform(tf.transform,g)
                if not math.isfinite(g.length) or g.length<=0:raise ValueError("nonpositive geometry")
                geom_total+=g.length/1000
                valid_rows.append((fn,a,g))
            except (ValueError,TypeError,AttributeError):
                status["invalid_geometry"]+=1
        if len(valid_rows)+status["invalid_geometry"]!=n:raise ValueError("record census mismatch")
        out={"paper_track_sections":paper,"v1_dbf_rows":n,"difference":n-paper,
             "unique_nonempty_fid":len(fids)-int("" in fids),"missing_fid_rows":fids[""],
             "duplicate_fid_excess_rows":sum(max(0,v-1) for k,v in fids.items() if k),
             "distinct_line_codes":len(linecodes),"distinct_line_section_pairs":len(sectionkeys),
             "unique_normalized_geometries":len(gh),"duplicate_normalized_geometry_excess":sum(v-1 for v in gh.values()),
             "unique_geometry_and_all_non_fid_attributes":len(rowhash),"exact_non_fid_duplicate_excess":sum(v-1 for v in rowhash.values()),
             "provisional_line_code_closure_year_pairs":len(events),"row_checks":dict(status),
             "date_value_classes":{k:dict(v) for k,v in date_tokens.items()},
             "closure_decade_counts":dict(sorted(closure_years.items())),
             "projected_row_length_sum_km":round(geom_total,3),"encoding":enc,"fields":fields,
             "member_sha256":{k:hashlib.sha256(b).hexdigest() for k,b in parts.items()}}
        return out

def main_audit():
    url=s.DV+"/api/datasets/:persistentId/?persistentId="+urllib.parse.quote(s.PID,safe=":")
    raw,receipt=get(url,6_000_000);j=json.loads(raw)
    if j.get("status")!="OK":raise ValueError("metadata unavailable")
    v=j["data"]["latestVersion"]
    if (v["versionNumber"],v["versionMinorNumber"])!=(1,0):raise ValueError("release changed")
    ff={f["dataFile"]["filename"]:f for f in v["files"]}
    def fetch(fn):
        f=ff[fn];df=f["dataFile"]
        if f.get("restricted") or df.get("restricted"):raise ValueError("restricted data")
        b,rc=get(s.DV+"/api/access/datafile/"+str(df["id"]),12_000_000)
        c=df.get("checksum",{});alg=str(c.get("type","")).lower().replace("-","")
        if alg not in ("md5","sha1","sha256","sha512") or hashlib.new(alg,b).hexdigest().lower()!=str(c.get("value","")).lower():
            raise ValueError("published checksum mismatch")
        return b,{**rc,"file_id":df["id"],"published_checksum":c}
    readme,rm=fetch("Readme.txt")
    lines=readme.decode("utf-8-sig",errors="replace").splitlines()
    # Record which semantic concepts appear, not full copyrighted README.
    concepts={k:sum(bool(re.search(pattern,x,re.I)) for x in lines) for k,pattern in {
        "reopening":r"reopen|reapertura","second_closure":r"second time|second closure|segunda",
        "zero_or_null":r"\bzero\b|\bnull\b|missing|blank|empty|sin datos",
        "validation":r"validat|robust|consisten"}.items()}
    res={}
    for fn,paper in s.FILES.items():
        try:
            data,rc=fetch(fn)
            cpgname=next((n for n in ff if n.lower()==(fn[:-4]+".cpg").lower()),None)
            cpg,cpgr=fetch(cpgname) if cpgname else (None,None)
            a=audit_zip(data,fn,paper,cpg)
            res[fn]={"status":"audited","receipt":rc,"cpg_receipt":cpgr,"audit":a}
        except Exception as e:res[fn]={"status":"blocked","type":type(e).__name__,"reason":str(e)[:250]}
    return {"status":"V1 release reconciliation only","metadata":receipt,"readme":rm,"readme_concept_counts":concepts,
        "layers":res,"decision":"The paper reports track sections, the DBF may count feature representations. Do not assume identical units or change the old gate. Distinct geometries/IDs and date consistency can guide a new explicitly V1-specific cohort, but discrepancy is not resolved without author/source explanation.",
        "limits":["no biological or satellite outcomes","no route names, feature IDs, geometries or raw data exported",
                  "candidate counts only, missing reopening semantics not inferred","paper https://doi.org/10.1017/S021261092400003X; data https://doi.org/10.34810/data917, CC BY 4.0"]}

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",default="results/railstatus/spain_v1_reconciliation.json")
    a=p.parse_args()
    r={"version":"spain-v1-source-reconciliation-v1","result":safe_source(main_audit)}
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False))
    print(json.dumps(r,ensure_ascii=False,allow_nan=False),flush=True)
