#!/usr/bin/env python3
"""Full NHM ADEPT source retrieval, exact-name match and real 2025 observed
plant trait coverage. Source preflight only; no ecological model or imputation.
No raw traits, species IDs, or observer IDs published.
"""
import argparse,csv,gzip,hashlib,io,json,math,re,time,unicodedata
import urllib.request,urllib.parse,urllib.error
from collections import Counter,defaultdict
from pathlib import Path
from audit_original_NHM_global_functional_trait_coverage_v372 import (
    SOURCE_RESOURCE,EXPECTED_TOTAL,COLUMNS,source_value)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    source_client,open_original_source,taxonomy_source,source_coarse_cell,
    OBS_ETAG,TAXA_ETAG,OBS_BYTES,TAXA_BYTES,YEAR,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS)
ROOT="https://data.nhm.ac.uk"
CSV=ROOT+"/dataset/9d05fda3-e22c-4919-bc29-76dd09cfa903/resource/"+SOURCE_RESOURCE+"/download/combined.csv"
API=ROOT+"/api/3/action/datastore_search"
DUMP=ROOT+"/datastore/dump/"+SOURCE_RESOURCE+"?format=csv"
SEED_LENGTH="seed max_ length [mm]"
HEIGHT="plant max_ height [m]"
FIELDS=tuple(dict.fromkeys(COLUMNS+("seed max_ width [mm]","seed min_ width [mm]","clonality","perennial organ")))
CAT=("life form","habit","dispersal mode","clonality")
NUM=(HEIGHT,SEED_LENGTH,"seed max_ diameter [mm]","seed max_ width [mm]")
# Exact source tokens only, not substrings. Unclassified values remain unknown.
ANNUAL={"annual","annual herb","annual herbs"}
PERENNIAL={"perennial","perennial herb","perennial herbs"}
WOODY={"tree","trees","shrub","shrubs","woody"}
W=100000

def norm(v):
    return " ".join(unicodedata.normalize("NFKC",str(v or "")).strip().casefold().split())
def number(v):
    if v is None or isinstance(v,(list,dict,bool)):return None
    s=str(v).strip()
    if not re.fullmatch(r"[+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",s):return None
    x=float(s)
    return x if math.isfinite(x) and x>0 else None
def scalar_label(v):
    # Do not print descriptions or arbitrary long raw text from publisher.
    if isinstance(v,(list,dict)):return None
    s=norm(v)
    return s if source_value(s) and len(s)<=40 and re.fullmatch(r"[a-z][a-z ,;/|\-]*",s) else None
def life_class(r):
    labels={scalar_label(r.get(k)) for k in ("life form","habit")}
    hits=[name for name,values in (("annual",ANNUAL),("perennial",PERENNIAL),("woody",WOODY)) if labels&values]
    return hits[0] if len(hits)==1 else None

def request(url,cap=50_000_000,tries=3):
    last=None
    for i in range(tries):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"tetsu-NHM-exact-source-replay/1.0","Accept":"*/*"})
            with urllib.request.urlopen(req,timeout=120) as r:
                b=r.read(cap+1)
            if len(b)>cap:raise ValueError("Source exceeds cap")
            return b
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            last=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in (408,429,500,502,503,504):raise
            if i+1<tries:time.sleep(2**i)
    raise RuntimeError("Publisher source unavailable") from last

def validate_rows(rows):
    if len(rows)!=EXPECTED_TOTAL:raise ValueError("Original row count mismatch")
    if any(not set(FIELDS)<=set(r) for r in rows):raise ValueError("Original columns missing")
    if any(not norm(r["taxon"]) for r in rows):raise ValueError("Missing source taxon")
    # Detect a repeated API page or truncated/replicated export using content,
    # even when source IDs are not numerical. Identical rows cannot be treated
    # as independent records without source-level documentation.
    fingerprints=[hashlib.sha256(json.dumps({k:r[k] for k in FIELDS},sort_keys=True,ensure_ascii=False).encode()).hexdigest() for r in rows]
    if len(set(fingerprints))!=len(rows):raise ValueError("Repeated projected source records")
    return hashlib.sha256("\n".join(sorted(fingerprints)).encode()).hexdigest()

def bulk(url):
    b=request(url)
    if b.startswith(b"\x1f\x8b"):b=gzip.decompress(b)
    rows=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"),newline="")))
    manifest=validate_rows(rows)
    return rows,{"method":"publisher CSV","url":url,"source_bytes":len(b),"source_SHA256":hashlib.sha256(b).hexdigest(),"normalized_projected_manifest":manifest}

def paginated():
    rows=[];pages=[];ids=set();offset=0;id_types=Counter()
    while offset<EXPECTED_TOTAL:
        params={"resource_id":SOURCE_RESOURCE,"limit":1000,"offset":offset,"sort":"_id asc",
                "fields":",".join(("_id",)+FIELDS),"include_total":"true"}
        raw=request(API+"?"+urllib.parse.urlencode(params),cap=25_000_000)
        z=json.loads(raw)
        if z.get("success") is not True:raise ValueError("DataStore request failed")
        q=z["result"];rr=q["records"]
        if q.get("total")!=EXPECTED_TOTAL or not isinstance(rr,list) or not rr or len(rr)>1000:raise ValueError("DataStore page or total changed")
        for r in rr:
            if "_id" not in r or r["_id"] is None or str(r["_id"])=="":raise ValueError("Missing stable source ID")
            rid=str(r["_id"])
            if rid in ids:raise ValueError("Source ID repeated across offset pages")
            ids.add(rid);id_types[type(r["_id"]).__name__]+=1
            if not set(FIELDS)<=set(r):raise ValueError("Selected DataStore columns missing")
            rows.append({k:r[k] for k in FIELDS})
        offset+=len(rr);pages.append(hashlib.sha256(raw).hexdigest())
        if offset%10000<1000 or offset==EXPECTED_TOTAL:print(json.dumps({"NHM_rows":offset}),flush=True)
    manifest=validate_rows(rows)
    if len(ids)!=EXPECTED_TOTAL:raise ValueError("Incomplete stable source IDs")
    return rows,{"method":"publisher stable-ID sorted offset pagination","pages":len(pages),"id_types":dict(id_types),
                 "page_manifest_SHA256":hashlib.sha256("\n".join(pages).encode()).hexdigest(),"normalized_projected_manifest":manifest}

def retrieve():
    errors=[]
    for label,fn in (("direct CSV",lambda:bulk(CSV)),("DataStore CSV dump",lambda:bulk(DUMP)),("stable offset API",paginated)):
        try:
            rows,qa=fn()
            return rows,dict(qa,earlier_route_failures=errors)
        except Exception as e:
            # No response body or taxon/ID contents in logs.
            msg=str(e)[:250]
            errors.append({"route":label,"exception":type(e).__name__,"message":msg})
            print(json.dumps({"source_route_failed":errors[-1]}),flush=True)
    raise RuntimeError("No full publisher source passed completeness: "+json.dumps(errors))

def active_names(client,canonical,expected_sha):
    rows,raw,stream=open_original_source(client,"taxa.csv.gz",TAXA_ETAG,TAXA_BYTES)
    names=defaultdict(list)
    try:
        for r in rows:
            tid=int(r["taxon_id"])
            if canonical.get(tid)==tid and r["rank"].strip().lower()=="species":
                names[norm(r["name"])].append(tid)
    finally:stream.close()
    if raw.bytes!=TAXA_BYTES or raw.hash.hexdigest()!=expected_sha:raise ValueError("Frozen taxon source SHA mismatch")
    if sum(map(len,names.values()))!=sum(k==v for k,v in canonical.items()):raise ValueError("Incomplete canonical species names")
    return names

def match(rows,names):
    by=defaultdict(list);qa=Counter();labels={k:Counter() for k in CAT}
    for r in rows:
        by[norm(r["taxon"])].append(r)
        for k in FIELDS:
            if source_value(r.get(k)):qa["source_available_"+k]+=1
        for k in CAT:
            val=scalar_label(r.get(k))
            if val is not None:labels[k][val]+=1
        for k in NUM:
            if number(r.get(k)) is not None:qa["source_positive_scalar_"+k]+=1
    joined={}
    for name,rr in by.items():
        if len(rr)!=1:
            qa["source_ambiguous_normalized_names"]+=1;continue
        ids=names.get(name,[])
        if len(ids)!=1:
            qa["unmatched_or_ambiguous_iNat_names"]+=1;continue
        r=rr[0];tid=ids[0]
        if tid in joined:raise ValueError("Species joined twice")
        row={"matched":True,"life_class":life_class(r)}
        row.update({k:number(r.get(k)) for k in NUM})
        row["dispersal_annotated"]=source_value(r.get("dispersal mode"))
        row["life_annotated"]=source_value(r.get("life form")) or source_value(r.get("habit"))
        joined[tid]=row
    qa["distinct_normalized_source_names"]=len(by)
    qa["exact_unique_active_iNat_matches"]=len(joined)
    qa["exact_matches_classified_lifeform"]=sum(r["life_class"] is not None for r in joined.values())
    for k in NUM:qa["exact_matches_positive_scalar_"+k]=sum(r[k] is not None for r in joined.values())
    qa["exact_matches_dispersal_annotated"]=sum(r["dispersal_annotated"] for r in joined.values())
    # Only frequent short category labels are exported to guide a later
    # controlled coding dictionary. No species-by-trait table is exported.
    return joined,dict(qa),{k:{a:n for a,n in v.most_common(30) if n>=10} for k,v in labels.items()}

def scan(client,canonical,plant):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    data=defaultdict(lambda:[0,set(),set()]);qa=Counter()
    try:
        fields=set(rows.fieldnames or [])
        if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=fields:raise ValueError("Source schema drift")
        if {"geoprivacy","coordinates_obscured","taxon_geoprivacy"}&fields:raise ValueError("Source geoprivacy drift")
        for r in rows:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"source_scan":dict(qa)}),flush=True)
            if not r["observed_on"].startswith(YEAR+"-"):continue
            if r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:tid=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if tid is None or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k is None or source_coarse_cell(lon,lat,tf,200000) is None:continue
            qa["plants"]+=1
            b=data[k];b[0]+=1;b[1].add(tid)
            if r["observer_id"]:b[2].add(int(r["observer_id"]))
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or qa["rows"]!=288721938 or qa["plants"]!=7827936 or raw.hash.hexdigest()!=plant["source_summary"]["source_iNat_original_2026_10_27_observations_gzip_SHA256_full"]:
        raise ValueError("Exact original observation replay failed")
    old={ (r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]):
        (r["source_year_2025_photo_linked_plants_observation_count"],r["source_2025_plant_species_richness"],r["source_2025_distinct_observers_with_photo_linked_plant_records"])
        for r in plant["source_grids"][str(W)]["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]}
    if old!={k:(v[0],len(v[1]),len(v[2])) for k,v in data.items()}:raise ValueError("Exact v369 grid replay failed")
    return data,dict(qa)

def counts(species,joined):
    j=[joined[t] for t in species if t in joined]
    known=[r["life_class"] for r in j if r["life_class"] is not None]
    out={"species":len(species),"matched":len(j),"life_annotated":sum(r["life_annotated"] for r in j),
         "life_classified":len(known),"annual":known.count("annual"),"perennial":known.count("perennial"),"woody":known.count("woody"),
         "dispersal_annotated":sum(r["dispersal_annotated"] for r in j)}
    for k in NUM:
        v=[r[k] for r in j if r[k] is not None]
        out["known_"+k]=len(v)
        out["mean_log_"+k]=sum(map(math.log,v))/len(v) if v else None
    return out

def main():
    p=argparse.ArgumentParser()
    for n in ("plant","audit","grid"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();plant=json.loads(a.plant.read_text())
    rows,source=retrieve();client=source_client()
    canonical,_,_,_,sha=taxonomy_source(client)
    if sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Original taxonomy mismatch")
    names=active_names(client,canonical,sha);joined,qa,labels=match(rows,names)
    audit={"version":"railstatus-NHM-exact-trait-preflight-v1","DOI":"10.5519/p3dm31kc","source_resource":SOURCE_RESOURCE,
         "publisher_rows":len(rows),"source":source,"taxonomy_sha256":sha,"match_qa":qa,"frequent_short_categories":labels,
         "scope":"exact unique normalized full scientific name, not synonym resolution; strict positive scalar numerical values; conservative unvalidated life-form categories"}
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False))
    print(json.dumps(audit,allow_nan=False),flush=True)
    # Avoid expensive scan if not enough unique species can be linked.
    if len(joined)<1000:
        print("Insufficient exact-name source support; no observation scan",flush=True);return
    del rows,names
    data,scanqa=scan(client,canonical,plant)
    all_species=set().union(*(v[1] for v in data.values()))
    overall=counts(all_species,joined)
    cells=[]
    for (i,j),(n,s,o) in sorted(data.items()):
        if n<50 or len(o)<5 or len(s)<20:continue
        cells.append({"key":[i,j],"observations":n,"observers":len(o),**counts(s,joined)})
    gates={}
    for k in ("life_classified","dispersal_annotated",*("known_"+k for k in NUM)):
        gates[k]={"cells_with_at_least_20_known_and_half_coverage":sum(r[k]>=20 and r[k]/r["species"]>=.5 for r in cells),
                  "cells_with_at_least_20_known_and_quarter_coverage":sum(r[k]>=20 and r[k]/r["species"]>=.25 for r in cells)}
    result={"version":"railstatus-NHM-observed-trait-coverage-v1","source":audit,"observation_replay":scanqa,
         "overall_observed_species_coverage":overall,"n_eligible_cells":len(cells),"coverage_gates":gates,"cells":cells,
         "limits":["no synonym matching; exact name is not validated WFO ID equivalence","positive scalar parser deliberately excludes multivalue text","unclassified means missing, not absence","trait annotated dispersal is not validated dispersal guild","coarse whole-landscape observations, not trackside flora","no railway association estimated"]}
    a.grid.parent.mkdir(parents=True,exist_ok=True);a.grid.write_text(json.dumps(result,allow_nan=False))
    print(json.dumps({k:v for k,v in result.items() if k not in ("cells","source")},allow_nan=False),flush=True)
if __name__=="__main__":main()
