#!/usr/bin/env python3
"""HWI contrasts within fixed AVONET habitat x migration strata.
Use same exact original observation source and replay all prior HWI guilds.
No occurrence, taxon, or observer identities exported.
"""
import argparse,json,hashlib,math
from collections import Counter,defaultdict
from pathlib import Path
from scan_world_railstatus_bird_avonet_guilds import (
    validate,csv_rows,DERIVATIVE,request,names,match,norm,num,group_pair,
    YEARS,W,MIN_GUILD_HALF,pair_key,taxonomy,source_client,open_original_source,
    source_coarse_cell,OBS_ETAG,OBS_BYTES,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS,
    bucket,add,pair_rows,SALT)
SOURCE_SHA="a8f0e0c889a05485f306029433cd45068836324e2bc6cf5b1ed8e34417f3d5fe"
STRATA=tuple((h,m) for h in (1,2,3) for m in (1,2,3))
AXES={f"h{h}_m{m}":(f"h{h}_m{m}_low",f"h{h}_m{m}_high") for h,m in STRATA}
HWI=("low_HWI","high_HWI")
GUILDS=HWI+tuple(g for v in AXES.values() for g in v)

def conditioned_match(rows,n):
    old,qa=match(rows,n)
    if qa["matched"]!=9872 or qa["HWI_median_matched_external_pool"]!=21.0:
        raise ValueError("Exact prior AVONET taxonomy/HWI support changed")
    out={};counts=Counter()
    for r in rows:
        ids=n.get(norm(r["Species1"]),[])
        if len(ids)!=1:continue
        t=ids[0];h,m,v=map(num,(r["Habitat.Density"],r["Migration"],r["Hand-Wing.Index"]))
        g=tuple(x for x in old[t] if x in HWI)
        if v is not None and 0<=v<=100:
            expected=("high_HWI" if v>=21 else "low_HWI",)
            if g!=expected:raise ValueError("HWI class replay mismatch")
            if h in (1,2,3) and m in (1,2,3):
                g=g+(AXES[f"h{int(h)}_m{int(m)}"][int(v>=21)],)
        out[t]=g;counts.update(g)
    if any(counts[g]!=qa["guild_species_counts"][g] for g in HWI):
        raise ValueError("Original HWI denominators changed")
    return out,{"original_match":qa,"joint_species_counts":dict(counts)}

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

def main():
    p=argparse.ArgumentParser()
    for n in ("split","annual","prior_guild","out","audit"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();read=lambda x:json.loads(x.read_text())
    src,annual,prior=map(read,(a.split,a.annual,a.prior_guild))
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or prior["version"]!="railstatus-avonet-bird-guild-observer-beta-v1" or not prior["exact_full_pair_replay"] or prior["guild_half_gate"]!=list(MIN_GUILD_HALF):
        raise ValueError("Wrong original artifacts")
    for k in ("obs_sha256","taxa_sha256","split_seed"):
        if len({src[k],annual[k],prior[k]})!=1:raise ValueError("Source identity mismatch")
    if src["split_seed"]!=SALT.decode() or prior["source"]["source"]["source_SHA256"]!=SOURCE_SHA:
        raise ValueError("Frozen source changed")
    previous={y:{pair_key(p):p for p in prior["pairs"][y]} for y in YEARS}
    full_previous={y:{pair_key(p):p for p in annual["pairs"][y]["birds"]} for y in YEARS[:2]}
    full_previous["2025"]={pair_key(p):p for p in src["pairs"][str(W)]["birds"] if "split_similarity" in p}
    common=set.intersection(*(set(v) for v in full_previous.values()))
    if len(common)!=1645 or any(set(previous[y])!=common for y in YEARS):
        raise ValueError("Frozen original geography changed")
    cells={k for pair in common for k in pair}
    raw=request(DERIVATIVE,cap=10_000_000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError("Frozen AVONET derivative SHA changed")
    rows=validate(csv_rows(raw));client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("Original taxonomy changed")
    trait,qa=conditioned_match(rows,names(client,canonical,sha))
    audit={"source_SHA256":SOURCE_SHA,"source_rows":len(rows),"match":qa,"trait_source":prior["source"]["source"],
           "HWI_threshold":21.0,"joint_strata":AXES}
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False))
    print(json.dumps({"conditional_hwi_source_audit":audit}),flush=True)
    full,guild,scanqa,osha=scan(client,canonical,trait,cells)
    if osha!=src["obs_sha256"]:raise ValueError("Observation SHA mismatch")
    result={}
    for y in YEARS:
        replay={pair_key(p):p for p in pair_rows(full[y]) if pair_key(p) in common}
        if set(replay)!=common:raise ValueError("Original full pair support changed")
        out=[]
        for k in sorted(common):
            if json.loads(json.dumps(replay[k]))!=full_previous[y][k]:
                raise ValueError("Original full pair outcomes changed")
            x,z=k;rr={"a":x,"b":z,"guilds":{}}
            for g in GUILDS:
                A=guild[y][g].get(x);B=guild[y][g].get(z)
                m=group_pair(A,B) if A is not None and B is not None else None
                if g in HWI and json.loads(json.dumps(m))!=previous[y][k]["guilds"].get(g):
                    raise ValueError("Original HWI guild pair replay failed")
                if m is not None:rr["guilds"][g]=m
            out.append(rr)
        result[y]=out
        print(json.dumps({"year":y,"pairs":len(out),"stratum_support":{ax:sum(all(g in p["guilds"] for g in gs) for p in out) for ax,gs in AXES.items()}}),flush=True)
    obj={"version":"railstatus-avonet-conditional-hwi-pairs-v1","source":audit,"source_qa":scanqa,
         "obs_sha256":osha,"taxa_sha256":sha,"width":W,"split_seed":SALT.decode(),
         "n_all_year_common_candidates":len(common),"guild_half_gate":MIN_GUILD_HALF,
         "exact_full_pair_replay":True,"exact_hwi_pair_replay":True,"pairs":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
