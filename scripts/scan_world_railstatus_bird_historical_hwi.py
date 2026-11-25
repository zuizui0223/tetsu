#!/usr/bin/env python3
"""One heldback 2018-22 pooled temporal test, same frozen geography and HWI.
2025 records replay the exact verified original full and joint-trait aggregates.
No raw taxon, occurrence or observer identifiers exported.
"""
import argparse,json,hashlib,math
from collections import Counter,defaultdict
from pathlib import Path
from scan_world_railstatus_bird_conditional_hwi import (
    validate,csv_rows,DERIVATIVE,request,names,conditioned_match,group_pair,
    W,MIN_GUILD_HALF,pair_key,taxonomy,source_client,open_original_source,
    source_coarse_cell,OBS_ETAG,OBS_BYTES,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS,
    bucket,add,pair_rows,SALT,SOURCE_SHA,AXES,HWI,GUILDS)
PERIOD="2018_2022"
PERIODS=(PERIOD,"2025")
DISCOVERY=("2023","2024","2025")

def period(date):
    if not isinstance(date,str) or len(date)<10 or date[4]!="-" or date[7]!="-":return None
    y=date[:4]
    if y in ("2018","2019","2020","2021","2022"):return PERIOD
    if y=="2025":return "2025"
    return None

def scan(client,canonical,trait,cells):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    qa=Counter()
    full={y:defaultdict(bucket) for y in PERIODS}
    guild={y:{g:defaultdict(bucket) for g in GUILDS} for y in PERIODS}
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"scan":dict(qa)}),flush=True)
            y=period(r["observed_on"])
            if y is None or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
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
    for n in ("split","annual","conditional","out","audit"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();read=lambda x:json.loads(x.read_text())
    src,annual,prior=map(read,(a.split,a.annual,a.conditional))
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or prior["version"]!="railstatus-avonet-conditional-hwi-pairs-v1" or not prior["exact_full_pair_replay"] or not prior["exact_hwi_pair_replay"] or prior["guild_half_gate"]!=list(MIN_GUILD_HALF):
        raise ValueError("Wrong original artifacts")
    for k in ("obs_sha256","taxa_sha256","split_seed"):
        if len({src[k],annual[k],prior[k]})!=1:raise ValueError("Frozen source identity mismatch")
    if src["split_seed"]!=SALT.decode() or prior["source"]["source_SHA256"]!=SOURCE_SHA:raise ValueError("Frozen trait source changed")
    full_previous={y:{pair_key(p):p for p in annual["pairs"][y]["birds"]} for y in DISCOVERY[:2]}
    full_previous["2025"]={pair_key(p):p for p in src["pairs"][str(W)]["birds"] if "split_similarity" in p}
    common=set.intersection(*(set(v) for v in full_previous.values()))
    previous={y:{pair_key(p):p for p in prior["pairs"][y]} for y in DISCOVERY}
    if len(common)!=1645 or any(set(previous[y])!=common for y in DISCOVERY):raise ValueError("Frozen original candidate geography changed")
    cells={k for pair in common for k in pair}
    raw=request(DERIVATIVE,cap=10_000_000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError("Frozen AVONET derivative SHA changed")
    rows=validate(csv_rows(raw));client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("Original taxonomy changed")
    trait,qa=conditioned_match(rows,names(client,canonical,sha))
    if qa!=prior["source"]["match"]:raise ValueError("Exact trait membership count replay failed")
    audit={"source_SHA256":SOURCE_SHA,"source_rows":len(rows),"match":qa,"HWI_threshold":21.0,
           "heldback_years":[2018,2019,2020,2021,2022],"fixed_candidate_pairs":len(common)}
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False))
    print(json.dumps({"historical_hwi_source_audit":audit}),flush=True)
    full,guild,scanqa,osha=scan(client,canonical,trait,cells)
    if osha!=src["obs_sha256"]:raise ValueError("Observation SHA mismatch")
    result={}
    for y in PERIODS:
        replay={pair_key(p):p for p in pair_rows(full[y]) if pair_key(p) in common and "split_similarity" in p}
        if y=="2025" and set(replay)!=common:raise ValueError("2025 full pair support changed")
        out=[]
        for k in sorted(replay):
            if y=="2025" and json.loads(json.dumps(replay[k]))!=full_previous[y][k]:raise ValueError("2025 full pair replay failed")
            x,z=k;rr={"a":x,"b":z,"guilds":{}}
            for g in GUILDS:
                A=guild[y][g].get(x);B=guild[y][g].get(z)
                m=group_pair(A,B) if A is not None and B is not None else None
                if y=="2025" and json.loads(json.dumps(m))!=previous[y][k]["guilds"].get(g):
                    raise ValueError("2025 exact conditional guild replay failed")
                if m is not None:rr["guilds"][g]=m
            out.append({"full":replay[k],**rr})
        result[y]=out
        print(json.dumps({"period":y,"full_pair_support":len(out),"HWI_pair_support":sum(all(g in p["guilds"] for g in HWI) for p in out),
             "joint_stratum_support":{ax:sum(all(g in p["guilds"] for g in gs) for p in out) for ax,gs in AXES.items()}}),flush=True)
    obj={"version":"railstatus-avonet-heldback-2018-22-hwi-pairs-v1","source":audit,"source_qa":scanqa,
         "obs_sha256":osha,"taxa_sha256":sha,"width":W,"split_seed":SALT.decode(),
         "n_frozen_candidates":len(common),"guild_half_gate":MIN_GUILD_HALF,
         "exact_2025_full_and_joint_guild_replay":True,"pairs":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
