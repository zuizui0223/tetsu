#!/usr/bin/env python3
"""One frozen original-source scan: earlier-year butterfly edge affinity.
Only coarse cell summaries, no individual occurrence/taxon/observer IDs exported."""
import argparse,json,math,hashlib,time
from pathlib import Path
from collections import Counter,defaultdict
from scan_world_railstatus_disturbance_traits import (
    taxonomy,source_client,open_original_source,source_coarse_cell,OBS_ETAG,OBS_BYTES,MAX_UNCERT_M,
    PROGRESS_EVERY_SOURCE_ROWS,taxonomy_metadata,load_lep,family_refs,summarize,
    WIDTHS,LEP_SHA,extract,bucket,add,full_counts,split_counts,SALT,MIN_FULL,MIN_SPLIT)
YEARS=("2023","2024","2025")
G="insects"
TRAIT="edge_affinity"

def validate_close(a,b):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise ValueError("Trait summary schema changed")
        for k in a:validate_close(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        if not isinstance(b,(list,tuple)) or len(a)!=len(b):raise ValueError("Trait count schema changed")
        for x,y in zip(a,b):validate_close(x,y)
    elif isinstance(a,(int,float)) and isinstance(b,(int,float)):
        if not math.isclose(a,b,abs_tol=1e-10,rel_tol=1e-10):raise ValueError("Exact trait replay mismatch")
    elif a!=b:raise ValueError("Trait summary replay mismatch")

def full_gate(f,h):
    return all(x>=y for x,y in zip(f,MIN_FULL)) and all(all(x>=y for x,y in zip(z,MIN_SPLIT)) for z in h)

def one_summary(v,t,traits,fam,refs):
    f=full_counts(v);h=tuple(tuple(x) for x in split_counts(v))
    tt=[{"counts":c,"traits":{TRAIT:summarize(s,traits,fam,refs,G)[TRAIT]}} for s,c in zip([t[1],*t[4]],[full_counts(t),*split_counts(t)])]
    if any(any(x>y for x,y in zip(q["counts"],p)) for q,p in zip(tt,(f,*h))):raise ValueError("Target counts exceed insects")
    return {"full":f,"halves":h,"target_full":tt[0],"target_halves":tt[1:]}

def scan(client,canonical,allowed,candidates):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    fields=set(rr.fieldnames or [])
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=fields:raise ValueError("Original schema changed")
    if {"geoprivacy","coordinates_obscured","taxon_geoprivacy"}&fields:raise ValueError("Geoprivacy schema changed")
    data={w:{y:defaultdict(bucket) for y in YEARS} for w in WIDTHS}
    target={w:{y:defaultdict(bucket) for y in YEARS} for w in WIDTHS}
    qa=Counter();t0=time.monotonic()
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"edge_temporal_scan":dict(qa),"seconds":round(time.monotonic()-t0,1)}),flush=True)
            y=r["observed_on"][:4]
            if y not in YEARS or not r["observed_on"].startswith(y+"-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0]!=G or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            t=hit[1];k1=source_coarse_cell(lon,lat,tf,100000);k2=source_coarse_cell(lon,lat,tf,200000)
            if k1 is None or k2 is None:continue
            qa[y+"_global_insects"]+=1;o=int(r["observer_id"]) if r["observer_id"] else None
            for w,k in ((100000,k1),(200000,k2)):
                if k not in candidates[w]:continue
                add(data[w][y][k],t,o)
                if t in allowed:add(target[w][y][k],t,o)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or qa["rows"]!=288721938:raise ValueError("Incomplete original census")
    return data,target,dict(qa),raw.hash.hexdigest()

def main():
    p=argparse.ArgumentParser()
    for k in ("traits","split","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();old=json.loads(a.traits.read_text());src=json.loads(a.split.read_text())
    if old["version"]!="railstatus-disturbance-trait-cells-v1" or not old["exact_full_and_split_replay"] or old["source"]["source_SHA256"]["LepTraits"]!=LEP_SHA or any(old[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong prior trait data")
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or src["split_seed"]!=SALT.decode() or src["full_gate_observations_observers_species"]!=list(MIN_FULL) or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT):raise ValueError("Wrong prior split data")
    cand={w:extract(src["pairs"][str(w)][G]) for w in WIDTHS}
    oldmap={w:{tuple(r["key"]):r for r in old["grids"][str(w)][G]} for w in WIDTHS}
    for w in WIDTHS:
        if {k:(tuple(r["full"]),tuple(tuple(z) for z in r["halves"])) for k,r in oldmap[w].items()}!=cand[w]:raise ValueError("Prior cell counts changed")
    client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("Frozen taxonomy changed")
    names,allowed,fam,roots=taxonomy_metadata(client,canonical,sha)
    traits=load_lep(names[G],allowed[G]);refs=family_refs(traits,fam,G)
    pool=old["source"]["pool_counts"][G]
    if len(allowed[G])!=pool["eligible_taxa"] or len(traits)!=pool["matched_taxa"] or sum(r.get(TRAIT) is not None for r in traits.values())!=pool["known"][TRAIT] or roots!=old["source"]["target_roots"]:raise ValueError("Frozen butterfly matching changed")
    del names
    data,target,qa,osha=scan(client,canonical,allowed[G],cand)
    if osha!=src["obs_sha256"] or qa["2025_global_insects"]!=old["qa"][G]:raise ValueError("Original 2025 replay failed")
    result={}
    for w in WIDTHS:
        result[str(w)]={}
        if set(data[w]["2025"])!=set(cand[w]):raise ValueError("2025 candidate support changed")
        for y in YEARS:
            cells=[]
            for k,v in sorted(data[w][y].items()):
                f=full_counts(v);h=tuple(tuple(z) for z in split_counts(v))
                if y=="2025" and (f,h)!=cand[w][k]:raise ValueError("2025 full/half census mismatch")
                if not full_gate(f,h):continue
                q=one_summary(v,target[w][y].get(k,bucket()),traits,fam,refs)
                if y=="2025":
                    r=oldmap[w][k]
                    for x,z in zip([q["target_full"],*q["target_halves"]],[r["target_full"],*r["target_halves"]]):
                        validate_close(x["counts"],z["counts"]);validate_close(x["traits"][TRAIT],z["traits"][TRAIT])
                cells.append({"key":list(k),**q})
            result[str(w)][y]=cells
            print(json.dumps({"edge_annual_source":{"width":w,"year":y,"candidate_cells":len(cand[w]),"observed_cells":len(data[w][y]),"split_eligible_cells":len(cells)}}),flush=True)
    out={"version":"railstatus-butterfly-edge-annual-cells-v1","source_sha256":{"traits":hashlib.sha256(a.traits.read_bytes()).hexdigest(),"split":hashlib.sha256(a.split.read_bytes()).hexdigest()},"obs_sha256":osha,"taxa_sha256":sha,"LepTraits_SHA256":LEP_SHA,"split_seed":SALT.decode(),"exact_2025_full_half_trait_replay":True,"years":list(YEARS),"qa":qa,"grids":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
