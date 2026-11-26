#!/usr/bin/env python3
"""One original iNaturalist census scan; fixed observer/record sampling.
Only coarse cell summaries leave runner; no identifiers or coordinates exported.
"""
import argparse,hashlib,json,math,time
from array import array
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scan_world_railstatus_cross_taxa import (
    taxonomy,source_client,open_original_source,source_coarse_cell,
    OBS_ETAG,TAXA_ETAG,OBS_BYTES,YEAR,MAX_UNCERT_M,
    SOURCE_REGION_M,PROGRESS_EVERY_SOURCE_ROWS)
GROUPS=("plants","birds","insects")
K=20
MIN_RECORDS=50
REPS=256
SEED="tetsu-fixed-observer-richness-v1"

def counts(z):
    return tuple(int(z[k]) for k in ("observations","observers","species"))

def candidates(prior):
    out={}
    for w in SOURCE_REGION_M:
        gmap={}
        for g in GROUPS:
            q={}
            for r in prior["grids"][str(w)][g]:
                k=(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"])
                n=counts(r)
                if k in q or min(n)<0 or max(n[1:])>n[0]:raise ValueError("Bad prior grid census")
                if n[0]>=MIN_RECORDS and n[1]>=K:q[k]=n
            gmap[g]=q
        # Freeze all pairwise plant-comparator effort-supported cells before reading raw records.
        keep={g:set(gmap["plants"])&set(gmap[g]) for g in GROUPS[1:]}
        keep["plants"]=set.union(*keep.values())
        out[w]={g:{k:gmap[g][k] for k in keep[g]} for g in GROUPS}
    return out

def bucket():
    return array("q"),array("q")

def scan(client,canonical,cand):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    fs=set(rows.fieldnames or [])
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=fs:raise ValueError("Source schema drift")
    if {"geoprivacy","coordinates_obscured","taxon_geoprivacy"}&fs:raise ValueError("Geoprivacy schema changed")
    data={w:{g:defaultdict(bucket) for g in GROUPS} for w in SOURCE_REGION_M}
    qa=Counter();t0=time.monotonic()
    try:
        for r in rows:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:
                print(json.dumps({"scan":dict(qa),"seconds":round(time.monotonic()-t0,1)}),flush=True)
            if not r["observed_on"].startswith(YEAR+"-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:
                hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (ValueError,TypeError):continue
            if hit is None or hit[0] not in GROUPS or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            g,s=hit
            small=source_coarse_cell(lon,lat,tf,100000);large=source_coarse_cell(lon,lat,tf,200000)
            if small is None or large is None:continue
            o=int(r["observer_id"]) if r["observer_id"] else -1
            if o<0 and r["observer_id"]:raise ValueError("Negative source observer ID")
            qa[g+"_eligible_global"]+=1
            for w,k in ((100000,small),(200000,large)):
                if k in cand[w][g]:
                    a,b=data[w][g][k];a.append(o);b.append(s);qa[f"{w}_{g}_candidate_records"]+=1
    finally:stream.close()
    if qa["rows"]!=288721938 or raw.bytes!=OBS_BYTES:raise ValueError("Incomplete frozen source")
    return data,dict(qa),raw.hash.hexdigest()

def summarize(observer,species,expected,seed,reps=REPS,k=K):
    o=np.asarray(observer,dtype=np.int64);s=np.asarray(species,dtype=np.int64)
    if o.ndim!=1 or s.shape!=o.shape or np.any(o< -1) or np.any(s<=0):raise ValueError("Invalid candidate records")
    known=o>=0
    actual=(len(o),len(np.unique(o[known])),len(np.unique(s)))
    if actual!=tuple(expected):raise ValueError(f"Exact candidate full census mismatch {actual} != {expected}")
    o=o[known];s=s[known]
    ix=np.lexsort((s,o));o=o[ix];s=s[ix]
    _,start,n=np.unique(o,return_index=True,return_counts=True)
    if len(start)<k or len(s)<k:raise ValueError("Fixed effort gate changed")
    rng=np.random.default_rng(seed)
    ob=np.empty(reps);rec=np.empty(reps)
    for j in range(reps):
        u=rng.choice(len(start),size=k,replace=False)
        # One uniform eligible source record from each selected distinct observer.
        take=start[u]+(rng.random(k)*n[u]).astype(np.int64)
        ob[j]=len(np.unique(s[take]))
        rec[j]=len(np.unique(s[rng.choice(len(s),size=k,replace=False)]))
    def stats(x):
        v=np.log1p(x)
        return {"mean_species":float(x.mean()),"mean_log1p_species":float(v.mean()),
                "mc_se_log1p":float(v.std(ddof=1)/math.sqrt(reps))}
    return {"full":list(actual),"known_records":len(s),"observer_balanced":stats(ob),"record_balanced":stats(rec)}

def main():
    p=argparse.ArgumentParser()
    for k in ("prior","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();prior=json.loads(a.prior.read_text())
    if prior["version"]!="railstatus-cross-taxa-v1" or prior["obs_etag"]!=OBS_ETAG or prior["taxa_etag"]!=TAXA_ETAG or prior["year"]!=YEAR:
        raise ValueError("Wrong frozen census")
    cand=candidates(prior)
    print(json.dumps({"candidate_counts":{str(w):{g:len(x) for g,x in gs.items()} for w,gs in cand.items()}}),flush=True)
    client=source_client();canonical,tcounts,tsha=taxonomy(client)
    if tsha!=prior["taxa_sha256"] or tcounts!=prior["active_species_counts"]:raise ValueError("Frozen taxonomy mismatch")
    data,qa,osha=scan(client,canonical,cand)
    if osha!=prior["obs_sha256"] or any(qa[g+"_eligible_global"]!=prior["qa"][g] for g in GROUPS):
        raise ValueError("Original source replay failed")
    grids={}
    for w in SOURCE_REGION_M:
        grids[str(w)]={}
        for g in GROUPS:
            if set(data[w][g])!=set(cand[w][g]):raise ValueError("Candidate cell census changed")
            z=[]
            for k in sorted(cand[w][g]):
                o,s=data[w][g].pop(k)
                seed=int.from_bytes(hashlib.sha256(f"{SEED}:{w}:{g}:{k[0]}:{k[1]}".encode()).digest()[:8],"big")
                q=summarize(o,s,cand[w][g][k],seed)
                z.append({"key":list(k),**q})
            grids[str(w)][g]=z
            print(json.dumps({"sampled_grid":{"width_m":w,"group":g,"n_cells":len(z),
                "median_mc_se":float(np.median([r["observer_balanced"]["mc_se_log1p"] for r in z])) if z else None}}),flush=True)
    out={"version":"railstatus-fixed-observer-richness-v1","year":YEAR,"source_sha256":{"prior":hashlib.sha256(a.prior.read_bytes()).hexdigest(),
          "observations":osha,"taxa":tsha},"qa":qa,"K":K,"min_full_records":MIN_RECORDS,"reps":REPS,"seed":SEED,
          "exact_candidate_full_counts_replayed":True,"grids":grids,
          "limits":["observed species accumulation, not true richness","fixed 20 different observers one random record per observer; not field visits","no individual identifiers or coordinates exported"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
