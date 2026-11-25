#!/usr/bin/env python3
"""v416 condition HWI ranks on source taxonomy and ecological trait strata."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import scan_world_railstatus_bird_avonet_guilds as base
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_ancestry_ids
SOURCE_SHA="a8f0e0c889a05485f306029433cd45068836324e2bc6cf5b1ed8e34417f3d5fe"
AXES0=dict(base.AXES)
EXTRA_AXES={
 "family_relative":("family_rel_low","family_rel_high"),
 "family_absolute":("family_abs_low","family_abs_high"),
 "family_ecology_relative":("ecology_rel_low","ecology_rel_high"),
 "family_ecology_absolute":("ecology_abs_low","ecology_abs_high")}
AXES={**AXES0,**EXTRA_AXES}
MIN_STRATUM=4
def names_family(client,canonical,sha):
    rows,raw,stream=base.open_original_source(client,"taxa.csv.gz",base.TAXA_ETAG,base.TAXA_BYTES)
    names=defaultdict(list);families=set();anc={}
    try:
        for r in rows:
            t=int(r["taxon_id"]);rank=r["rank"].strip().lower()
            if rank=="family":families.add(t)
            if canonical.get(t)==("birds",t) and rank=="species":
                names[base.norm(r["name"])].append(t)
                anc[t]=tuple(map(int,source_ancestry_ids(r["ancestry"])))
    finally:stream.close()
    if raw.bytes!=base.TAXA_BYTES or raw.hash.hexdigest()!=sha:raise ValueError("Original taxonomy SHA mismatch")
    family={t:f for t,a in anc.items() if (f:=next((k for k in reversed(a) if k in families),None)) is not None}
    return names,family
def rank_groups(records,ecology):
    """Disjoint symmetric rank pairs in each source stratum, excluding ties."""
    groups=defaultdict(list)
    for t,(f,h,m,v) in records.items():
        if f is None or v is None:continue
        if ecology and (h not in (1,2,3) or m not in (1,2,3)):continue
        k=(f,h,m) if ecology else (f,)
        groups[k].append((v,t))
    low=set();high=set();kept=0;ties=0;nn=[]
    for k,ls in sorted(groups.items()):
        if len(ls)<MIN_STRATUM:continue
        ls.sort();n=0
        for i in range(len(ls)//2):
            x,y=ls[i],ls[-1-i]
            if y[0]<=x[0]:
                ties+=1;continue
            low.add(x[1]);high.add(y[1]);n+=1
        if n:kept+=1;nn.append(n)
    if low&high or len(low)!=len(high):raise ValueError("Invalid matched source labels")
    return low,high,{"eligible_strata":len(groups),"paired_strata":kept,"paired_species_per_label":len(low),
        "skipped_equal_HWI_pairs":ties,"median_species_pairs_per_stratum":float(np.median(nn)) if nn else None}
def conditioned(rows,names,family):
    orig,qa=base.match(rows,names)
    records={}
    for r in rows:
        ids=names.get(base.norm(r["Species1"]),[])
        if len(ids)!=1:continue
        t=ids[0];h,m,v=map(base.num,(r["Habitat.Density"],r["Migration"],r["Hand-Wing.Index"]))
        records[t]=(family.get(t),h,m,v if v is not None and 0<=v<=100 else None)
    if set(records)!=set(orig):raise ValueError("Original trait match denominator changed")
    out={t:list(v) for t,v in orig.items()};d={}
    median=qa["HWI_median_matched_external_pool"]
    for prefix,eco in (("family",False),("ecology",True)):
        lo,hi,z=rank_groups(records,eco)
        if len(lo)<100:raise ValueError("Insufficient source taxonomic pairs")
        for t in lo|hi:
            out[t].append(prefix+"_rel_"+("low" if t in lo else "high"))
            out[t].append(prefix+"_abs_"+("high" if records[t][3]>=median else "low"))
        z.update(absolute_low=sum(records[t][3]<median for t in lo|hi),absolute_high=sum(records[t][3]>=median for t in lo|hi))
        d[prefix]=z
    return {t:tuple(v) for t,v in out.items()},qa,d
def main():
    p=argparse.ArgumentParser()
    for n in ("split","annual","prior_guild","out","audit"):p.add_argument("--"+n.replace("_","-"),dest=n,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    src,annual,old=map(read,(a.split,a.annual,a.prior_guild))
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or old["version"]!="railstatus-avonet-bird-guild-observer-beta-v1" or not old["exact_full_pair_replay"] or any(src[k]!=annual[k] or src[k]!=old[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or src["split_seed"]!=base.SALT.decode():raise ValueError("Wrong original pair sources")
    # Reuse identical public derivative bytes, not a silently different original archive.
    import hashlib
    raw=base.request(base.DERIVATIVE,cap=10_000_000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA or old["source"]["source"]["source_SHA256"]!=SOURCE_SHA:raise ValueError("Prior exact AVONET derivative changed")
    rows=base.validate(base.csv_rows(raw));client=base.source_client();canonical,_,sha=base.taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("Taxa SHA mismatch")
    names,family=names_family(client,canonical,sha)
    traits,qa,cond=conditioned(rows,names,family)
    if qa!=old["source"]["match"]:raise ValueError("Original trait support mismatch")
    audit={"source":old["source"],"conditional_source_labels":cond,"matched_with_source_family":sum(t in family for t in traits),"min_stratum":MIN_STRATUM}
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False))
    print(json.dumps({"v416_source_QA":audit}),flush=True)
    prior={y:{base.pair_key(p):p for p in annual["pairs"][y]["birds"]} for y in base.YEARS[:2]}
    prior["2025"]={base.pair_key(p):p for p in src["pairs"][str(base.W)]["birds"] if "split_similarity" in p}
    common=set.intersection(*(set(v) for v in prior.values()))
    if len(common)!=1645 or len(common)!=old["n_all_year_common_candidates"]:raise ValueError("Prior candidate geography mismatch")
    og={y:{base.pair_key(p):p for p in old["pairs"][y]} for y in base.YEARS}
    if any(set(og[y])!=common for y in base.YEARS):raise ValueError("Prior guild geography changed")
    # Existing scan uses a module-level guild list; extend only this job.
    base.GUILDS=tuple(g for gs in AXES.values() for g in gs)
    cells={k for pair in common for k in pair}
    full,guild,scanqa,osha=base.scan(client,canonical,traits,cells)
    if osha!=src["obs_sha256"] or scanqa!=old["source_qa"]:raise ValueError("Exact source replay mismatch")
    result={}
    for y in base.YEARS:
        replay={base.pair_key(p):p for p in base.pair_rows(full[y]) if base.pair_key(p) in common}
        if set(replay)!=common:raise ValueError("Full pair support mismatch")
        out=[]
        for k in sorted(common):
            if json.loads(json.dumps(replay[k]))!=prior[y][k]:raise ValueError("Exact old full/split similarity mismatch")
            x,z=k;rr={"a":x,"b":z,"guilds":{}}
            for g in base.GUILDS:
                A=guild[y][g].get(x);B=guild[y][g].get(z)
                if A is None or B is None:continue
                m=base.group_pair(A,B)
                if m is not None:rr["guilds"][g]=m
            old_g={g:v for g,v in rr["guilds"].items() if g in tuple(g for gs in AXES0.values() for g in gs)}
            if json.loads(json.dumps(old_g))!=og[y][k]["guilds"]:raise ValueError("Original six guilds did not exactly replay")
            out.append(rr)
        result[y]=out
        print(json.dumps({"v416_year":y,"n":len(out),"both_guild_supported":{axis:sum(all(g in p["guilds"] for g in gs) for p in out) for axis,gs in AXES.items()}}),flush=True)
    obj={"version":"v416-conditional-HWI-pair-guilds","source":old["source"],"conditional_audit":audit,"source_qa":scanqa,
         "obs_sha256":osha,"taxa_sha256":sha,"width":base.W,"split_seed":base.SALT.decode(),"n_all_year_common_candidates":len(common),
         "exact_full_pair_replay":True,"exact_original_six_guild_replay":True,"guild_half_gate":base.MIN_GUILD_HALF,"pairs":result,
         "limits":["within source family/ecology rank groups balance source species counts only, not observed local composition",
                   "relative HWI differs from global HWI category; absolute groups reuse exactly same source taxa",
                   "family from original iNat hierarchy, not full phylogenetic distances",
                   "same observation platform, fixed 2025 geography and 2026 rail map; no raw identifiers exported"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
