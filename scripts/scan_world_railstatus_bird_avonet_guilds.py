#!/usr/bin/env python3
"""AVONET BirdLife traits and observer-disjoint bird guild beta.
Reuse fixed 2023-25 source pair geography, one exact iNat scan.
No raw occurrence, observer or species IDs are exported.
"""
import argparse,csv,hashlib,io,json,math,time,zipfile
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from audit_world_NHM_rail_trait_match import norm,request
from scan_world_railstatus_cross_taxa import taxonomy,source_client,open_original_source,source_coarse_cell,OBS_ETAG,TAXA_ETAG,OBS_BYTES,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import TAXA_BYTES
from scan_world_railstatus_observer_split_beta import bucket,add,pair_rows,overlap,split_counts,eligible,observer_side,SALT
from scan_world_railstatus_observer_annual_beta import pair_key
YEARS=("2023","2024","2025")
W=200000
AXES={"habitat":("dense","open"),"migration":("resident","migrant"),"HWI":("low_HWI","high_HWI")}
GUILDS=tuple(g for v in AXES.values() for g in v)
MIN_GUILD_HALF=(10,3,5)
OFFICIAL="https://figshare.com/ndownloader/files/34480865?private_link=b990722d72a26b5bfead"
DERIVATIVE_COMMIT="34baf9b9241d3c2d55b92287412294cc128b0f72"
DERIVATIVE="https://raw.githubusercontent.com/BrentPease1/alan/"+DERIVATIVE_COMMIT+"/data/traits/avonet.csv"
REQUIRED=("Species1","Habitat.Density","Migration","Hand-Wing.Index")

def csv_rows(raw):
    # Both comma and tab separated publisher tables are supported.
    txt=raw.decode("utf-8-sig")
    sep="\t" if txt.splitlines()[0].count("\t")>txt.splitlines()[0].count(",") else ","
    return list(csv.DictReader(io.StringIO(txt),delimiter=sep))

def validate(rows):
    if not 10500<=len(rows)<=11200 or not all(set(REQUIRED)<=set(r) for r in rows):
        raise ValueError("Not full BirdLife AVONET table or required columns absent")
    names=[norm(r["Species1"]) for r in rows]
    if not all(names) or len(set(names))!=len(names):raise ValueError("Duplicate or missing BirdLife names")
    return rows

def workbook(raw):
    import openpyxl
    w=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    try:
        for sheet in w.worksheets:
            if "birdlife" not in sheet.title.casefold():continue
            it=sheet.values;header=next(it);names=[str(x) if x is not None else "" for x in header]
            if set(REQUIRED)<=set(names):
                yield sheet.title,[dict(zip(names,row)) for row in it if any(x is not None for x in row)]
    finally:w.close()
def original_rows(raw):
    z=zipfile.ZipFile(io.BytesIO(raw))
    if "xl/workbook.xml" in z.namelist():
        candidates=list(workbook(raw))
    else:
        candidates=[]
        for m in z.infolist():
            n=m.filename.casefold()
            if m.file_size>150_000_000:continue
            if "birdlife" in n and n.endswith((".csv",".txt",".tsv")):
                rr=csv_rows(z.read(m))
                if rr and set(REQUIRED)<=set(rr[0]):candidates.append((m.filename,rr))
            elif n.endswith(".xlsx") and ("avonet" in n or "supplement" in n):
                candidates.extend((m.filename+":"+s,r) for s,r in workbook(z.read(m)))
    good=[]
    for name,rows in candidates:
        try:good.append((name,validate(rows)))
        except ValueError:continue
    if len(good)!=1:raise ValueError("Cannot identify unique full original BirdLife table")
    return good[0]
def retrieve():
    errors=[]
    try:
        raw=request(OFFICIAL,cap=150_000_000);name,rows=original_rows(raw)
        return rows,{"route":"original Figshare archive","source_url":OFFICIAL,"source_member":name,
                     "source_SHA256":hashlib.sha256(raw).hexdigest(),"source_bytes":len(raw)}
    except Exception as e:
        errors.append({"route":"original","exception":type(e).__name__,"message":str(e)[:180]})
        print(json.dumps({"AVONET_source_failure":errors[-1]}),flush=True)
    # Explicitly documented derivative only, not represented as original
    # publisher bytes. Its Git blob identity is frozen and checked.
    raw=request(DERIVATIVE,cap=10_000_000)
    blob=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
    if blob!="e210b459a740800551a6b4fa2f79c48f4ab98018":raise ValueError("Derivative Git blob mismatch")
    rows=validate(csv_rows(raw))
    return rows,{"route":"pinned published analysis derivative; not original archive bytes","source_url":DERIVATIVE,
                 "source_SHA256":hashlib.sha256(raw).hexdigest(),"git_blob_sha1":blob,"source_bytes":len(raw),"earlier_failures":errors}
def num(x):
    try:v=float(x)
    except (ValueError,TypeError):return None
    return v if math.isfinite(v) else None
def names(client,canonical,sha):
    rr,raw,stream=open_original_source(client,"taxa.csv.gz",TAXA_ETAG,TAXA_BYTES)
    n=defaultdict(list)
    try:
        for r in rr:
            tid=int(r["taxon_id"])
            if canonical.get(tid)==("birds",tid) and r["rank"].lower()=="species":
                n[norm(r["name"])].append(tid)
    finally:stream.close()
    if raw.bytes!=TAXA_BYTES or raw.hash.hexdigest()!=sha:raise ValueError("Taxonomy source SHA mismatch")
    return n
def match(rows,n):
    joined={};qa=Counter()
    for r in rows:
        ids=n.get(norm(r["Species1"]),[])
        if len(ids)!=1:continue
        h,m,v=map(num,(r["Habitat.Density"],r["Migration"],r["Hand-Wing.Index"]))
        joined[ids[0]]=(h if h in (1,2,3) else None,m if m in (1,2,3) else None,v if v is not None and 0<=v<=100 else None)
    known=[v[2] for v in joined.values() if v[2] is not None]
    if len(joined)<8000 or len(known)<7000:raise ValueError("Insufficient AVONET exact-name support")
    median=float(np.median(known))
    out={}
    for t,(h,m,v) in joined.items():
        g=[]
        if h in (1,3):g.append("dense" if h==1 else "open")
        if m in (1,3):g.append("resident" if m==1 else "migrant")
        if v is not None:g.append("high_HWI" if v>=median else "low_HWI")
        out[t]=tuple(g);qa.update(g)
    if min(qa[g] for g in GUILDS)<100:raise ValueError("Insufficient source trait class")
    return out,{"matched":len(joined),"HWI_median_matched_external_pool":median,"guild_species_counts":dict(qa),
                "habitat_known":sum(v[0] is not None for v in joined.values()),"migration_known":sum(v[1] is not None for v in joined.values()),"HWI_known":len(known)}
def scan(client,canonical,trait,cells):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    qa=Counter()
    full={y:defaultdict(bucket) for y in YEARS}
    guild={y:{g:defaultdict(bucket) for g in GUILDS} for y in YEARS}
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"scan":dict(qa)}),flush=True)
            y=r["observed_on"][:4]
            if y not in YEARS or not r["observed_on"].startswith(y+"-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0]!="birds" or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k not in cells:continue
            o=int(r["observer_id"]) if r["observer_id"] else None;t=hit[1]
            add(full[y][k],t,o);qa[y+"_bird_records"]+=1
            for g in trait.get(t,()):add(guild[y][g][k],t,o)
    finally:stream.close()
    if qa["rows"]!=288721938 or raw.bytes!=OBS_BYTES:raise ValueError("Incomplete original source")
    return full,guild,dict(qa),raw.hash.hexdigest()
def group_pair(a,b):
    ac,bc=split_counts(a),split_counts(b)
    if not all(eligible(v,MIN_GUILD_HALF) for v in (*ac,*bc)):return None
    vals={f"{i}{j}":overlap(a[4][i],b[4][j]) for i in (0,1) for j in (0,1)}
    return {"split_counts_a":ac,"split_counts_b":bc,"split_similarity":vals}
def main():
    p=argparse.ArgumentParser()
    for n in ("split","annual","out","audit"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();src=json.loads(a.split.read_text());annual=json.loads(a.annual.read_text())
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or any(src[k]!=annual[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or src["split_seed"]!=SALT.decode():raise ValueError("Wrong frozen source pair artifacts")
    prior={y:{pair_key(p):p for p in annual["pairs"][y]["birds"]} for y in YEARS[:2]}
    prior["2025"]={pair_key(p):p for p in src["pairs"][str(W)]["birds"] if "split_similarity" in p}
    common=set.intersection(*(set(v) for v in prior.values()))
    cells={k for pair in common for k in pair}
    rows,source=retrieve();client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("Taxonomy SHA mismatch")
    traits,qa=match(rows,names(client,canonical,sha));audit={"source":source,"match":qa,"source_rows":len(rows),"paper_DOI":"10.1111/ele.13898"}
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False))
    print(json.dumps({"bird_trait_source_audit":audit}),flush=True)
    full,guild,scanqa,osha=scan(client,canonical,traits,cells)
    if osha!=src["obs_sha256"]:raise ValueError("Observation SHA mismatch")
    result={}
    for y in YEARS:
        replay={pair_key(p):p for p in pair_rows(full[y]) if pair_key(p) in common}
        if set(replay)!=common:raise ValueError("Annual full pair replay support mismatch")
        out=[]
        for k in sorted(common):
            # JSON normalizes tuples to lists, matching prior serialized artifact.
            if json.loads(json.dumps(replay[k]))!=prior[y][k]:raise ValueError("Exact annual full/split bird similarity replay failed")
            x,z=k
            rr={"a":x,"b":z,"guilds":{}}
            for g in GUILDS:
                A=guild[y][g].get(x);B=guild[y][g].get(z)
                if A is None or B is None:continue
                m=group_pair(A,B)
                if m is not None:rr["guilds"][g]=m
            out.append(rr)
        result[y]=out
        print(json.dumps({"year":y,"all_year_common_candidates":len(out),"both_guild_supported":{axis:sum(all(g in x["guilds"] for g in gs) for x in out) for axis,gs in AXES.items()}}),flush=True)
    obj={"version":"railstatus-avonet-bird-guild-observer-beta-v1","source":audit,"source_qa":scanqa,
         "obs_sha256":osha,"taxa_sha256":sha,"width":W,"split_seed":SALT.decode(),
         "n_all_year_common_candidates":len(common),"exact_full_pair_replay":True,"guild_half_gate":MIN_GUILD_HALF,"pairs":result,
         "limits":["traits are species averages and primary habitat classifications, not actual individual movements","open habitat includes aquatic/open habitat types unless excluded by a future separate test","no eBird occurrence data","fixed 2025 geography and same original platform","only coarse aggregates"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
