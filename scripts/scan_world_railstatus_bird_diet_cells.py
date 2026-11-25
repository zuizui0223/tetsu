#!/usr/bin/env python3
"""Frozen AVONET trophic guild counts on existing disjoint plant/bird geography.
One exact source scan; only coarse per-cell counts, not IDs, are exported.
"""
import argparse,hashlib,json,math,time
from pathlib import Path
from collections import Counter,defaultdict
from scan_world_railstatus_bird_avonet_guilds import (
    validate,csv_rows,DERIVATIVE,request,names,match,norm,taxonomy,source_client,
    open_original_source,source_coarse_cell,OBS_ETAG,OBS_BYTES,MAX_UNCERT_M,
    PROGRESS_EVERY_SOURCE_ROWS,bucket,add,SALT)
from scan_world_railstatus_observer_split_beta import full_counts,split_counts,eligible,MIN_FULL,MIN_SPLIT
from analyze_world_railstatus_plant_bird_disjoint import extract,YEARS,GROUPS
from scan_world_railstatus_bird_conditional_hwi import SOURCE_SHA
GUILDS=("fruit_seed","invertivore","frugivore","granivore")
GATE=(10,3,5)
CATEGORIES={"frugivore","granivore","nectarivore","invertivore","vertivore","scavenger","omnivore","aquaticpredator","herbivoreaquatic","aquaticherbivore","herbivoreterrestrial","terrestrialherbivore"}

def category(x):
    return "".join(c for c in str(x or "").casefold() if c.isalpha())

def traits(rows,n):
    previous,qa=match(rows,n)
    if qa["matched"]!=9872 or qa["HWI_median_matched_external_pool"]!=21.:raise ValueError("Frozen AVONET match changed")
    if "Trophic.Niche" not in rows[0]:raise ValueError("Published trophic niche column missing")
    counts=Counter(category(r["Trophic.Niche"]) for r in rows)
    print(json.dumps({"AVONET_trophic_categories":dict(counts)}),flush=True)
    if set(counts)-CATEGORIES-{""}:raise ValueError("Unrecognized trophic niche code; inspect original metadata, do not infer")
    out={};tot=Counter()
    for r in rows:
        ids=n.get(norm(r["Species1"]),[])
        if len(ids)!=1:continue
        t=ids[0];cat=category(r["Trophic.Niche"])
        if cat in ("frugivore","granivore"):g=("fruit_seed",cat)
        elif cat=="invertivore":g=("invertivore",)
        else:g=()
        out[t]=g;tot.update(g)
    if len(out)!=len(previous) or min(tot[k] for k in GUILDS)<100:raise ValueError("Trophic guild match insufficient")
    return out,{"source_categories":dict(counts),"matched":len(out),"matched_guild_species":dict(tot),"previous_HWI_match":qa}

def candidates(src,annual):
    out={}
    for w in (100000,200000):
        for y in YEARS if w==200000 else ("2025",):
            g={k:extract(src["pairs"][str(w)][k] if y=="2025" else annual["pairs"][y][k]) for k in GROUPS}
            out[w,y]={k:g["birds"][k] for k in set(g["birds"])&set(g["plants"])}
    return out

def scan(client,canonical,trait,cand):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=set(rr.fieldnames or []):raise ValueError("Source schema changed")
    full={k:defaultdict(bucket) for k in cand};guild={k:{g:defaultdict(bucket) for g in GUILDS} for k in cand}
    qa=Counter()
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"scan":dict(qa)}),flush=True)
            y=r["observed_on"][:4]
            if y not in YEARS or not r["observed_on"].startswith(y+"-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0]!="birds" or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            t=hit[1];o=int(r["observer_id"]) if r["observer_id"] else None
            for w in (100000,200000) if y=="2025" else (200000,):
                k=source_coarse_cell(lon,lat,tf,w)
                if k not in cand[w,y]:continue
                add(full[w,y][k],t,o)
                qa[f"{w}_{y}_records"]+=1
                for g in trait.get(t,()):add(guild[w,y][g][k],t,o)
    finally:stream.close()
    if qa["rows"]!=288721938 or raw.bytes!=OBS_BYTES:raise ValueError("Incomplete exact original source")
    return full,guild,dict(qa),raw.hash.hexdigest()

def main():
    p=argparse.ArgumentParser()
    for k in ("split","annual","out","audit"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src=json.loads(a.split.read_text());annual=json.loads(a.annual.read_text())
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or any(src[k]!=annual[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or src["split_seed"]!=SALT.decode():raise ValueError("Frozen split sources changed")
    cand=candidates(src,annual)
    raw=request(DERIVATIVE,cap=10_000_000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError("Frozen AVONET derivative changed")
    rows=validate(csv_rows(raw));client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("Frozen taxonomy changed")
    trait,tqa=traits(rows,names(client,canonical,sha))
    audit={"source_SHA256":SOURCE_SHA,"source_url":DERIVATIVE,"source_rows":len(rows),"trait_match":tqa,"guild_gate_per_half":GATE}
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False))
    print(json.dumps({"diet_source_audit":audit}),flush=True)
    full,guild,qa,osha=scan(client,canonical,trait,cand)
    if osha!=src["obs_sha256"]:raise ValueError("Observation SHA changed")
    out={}
    for (w,y),expected in cand.items():
        rr=[]
        if set(full[w,y])!=set(expected):raise ValueError("Frozen cell support not exactly replayed")
        for k,(f,h) in sorted(expected.items()):
            z=full[w,y][k]
            if full_counts(z)!=f or tuple(split_counts(z))!=h:raise ValueError("Full or split source counts changed")
            result={"key":list(k),"full":list(f),"halves":h,"guilds":{}}
            for g in GUILDS:
                if k not in guild[w,y][g]:continue
                hs=split_counts(guild[w,y][g][k])
                if all(eligible(v,GATE) for v in hs):result["guilds"][g]=hs
            rr.append(result)
        out.setdefault(str(w),{})[y]=rr
        print(json.dumps({"width":w,"year":y,"exact_replay_cells":len(rr),"supported_guild_cells":{g:sum(g in r["guilds"] for r in rr) for g in GUILDS},"primary_joint":sum(all(g in r["guilds"] for g in ("fruit_seed","invertivore")) for r in rr)}),flush=True)
    obj={"version":"railstatus-avonet-trophic-disjoint-cell-counts-v1","source":audit,"source_qa":qa,"obs_sha256":osha,"taxa_sha256":sha,"split_seed":SALT.decode(),"exact_full_and_split_cell_replay":True,"grids":out}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
