#!/usr/bin/env python3
"""Trait-specific actual crossing test: exact original split-pair replay.
All four observer-endpoint halves must support both directional guilds.
Only coarse aggregates exported, no occurrence, observer or species IDs.
"""
import argparse,json,math,hashlib
from collections import defaultdict,Counter
from pathlib import Path
from scan_world_railstatus_disturbance_traits import (
    taxonomy,source_client,taxonomy_metadata,load_life,load_lep,load_birds,
    LIFE_SHA,LEP_SHA,SOURCE_SHA,WIDTHS,extract,OBS_ETAG,OBS_BYTES,
    open_original_source,source_coarse_cell,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS,
    bucket,add,full_counts,split_counts,SALT,MIN_FULL,MIN_SPLIT)
from scan_world_railstatus_observer_split_beta import overlap,eligible
AXES={"plants":("annual_nonwoody",),"insects":("edge_affinity","disturbance_affinity")}
GROUPS=tuple(AXES)
GATE=(10,3,5)

def key(p):return tuple(p["a"]),tuple(p["b"])
def counts(z):return tuple(tuple(x) for x in split_counts(z))
def categories(t,g,traits):
    r=traits[g].get(t,{})
    return [(axis,int(r[axis])) for axis in AXES[g] if r.get(axis) in (0.,1.)]
def pack_pair(aa,bb,dena,denb):
    # aa and bb are the two directional guilds (0 and 1).
    ac=[counts(x) for x in aa];bc=[counts(x) for x in bb]
    if not all(eligible(h,GATE) for x in ac+bc for h in x):return None
    coverage=[]
    for h in (0,1):
        for vv,den in ((ac,dena),(bc,denb)):
            n=vv[0][h][2]+vv[1][h][2]
            if not 0<n<=den[h]:raise ValueError("Guild denominator mismatch")
            coverage.append(n/den[h])
    return {"min_half_coverage":min(coverage),
            "guilds":[{"split_counts_a":ac[g],"split_counts_b":bc[g],
              "split_similarity":{f"{i}{j}":overlap(aa[g][4][i],bb[g][4][j]) for i in (0,1) for j in (0,1)}} for g in (0,1)]}

def main():
    p=argparse.ArgumentParser()
    for n in ("taxa","split","trait_cells","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();prior,src,old=map(lambda p:json.loads(p.read_text()),(a.taxa,a.split,a.trait_cells))
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_seed"]!=SALT.decode():raise ValueError("Wrong split")
    if old["version"]!="railstatus-disturbance-trait-cells-v1" or not old["exact_full_and_split_replay"] or any(len({x[k] for x in (prior,src,old)})!=1 for k in ("obs_sha256","taxa_sha256")):raise ValueError("Wrong trait cell artifact")
    if old["source"]["source_SHA256"]["LifeForm"]!=LIFE_SHA or old["source"]["source_SHA256"]["LepTraits"]!=LEP_SHA:raise ValueError("Wrong trait source")
    cand={w:{g:extract(src["pairs"][str(w)][g]) for g in GROUPS} for w in WIDTHS}
    tc={w:{g:{tuple(r["key"]):r for r in old["grids"][str(w)][g]} for g in GROUPS} for w in WIDTHS}
    for w in WIDTHS:
        for g in GROUPS:
            if {k:(tuple(v["full"]),tuple(tuple(q) for q in v["halves"])) for k,v in tc[w][g].items()}!=cand[w][g]:raise ValueError("Trait cell support changed")
    client=source_client();canonical,_,tsha=taxonomy(client)
    if tsha!=src["taxa_sha256"]:raise ValueError("Taxonomy hash changed")
    names,allowed,fam,roots=taxonomy_metadata(client,canonical,tsha)
    traits={"plants":load_life(names["plants"],allowed["plants"]),"insects":load_lep(names["insects"],allowed["insects"])}
    axes={g:{t:categories(t,g,traits) for t in traits[g]} for g in GROUPS}
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    full={w:{g:defaultdict(bucket) for g in GROUPS} for w in WIDTHS}
    target={w:{g:defaultdict(bucket) for g in GROUPS} for w in WIDTHS}
    guild={w:{g:{ax:{i:defaultdict(bucket) for i in (0,1)} for ax in AXES[g]} for g in GROUPS} for w in WIDTHS}
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES);qa=Counter()
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=set(rr.fieldnames or []):raise ValueError("Observation schema changed")
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"guild_pair_scan":dict(qa)}),flush=True)
            if not r["observed_on"].startswith("2025-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0] not in GROUPS or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            g,t=hit;k1=source_coarse_cell(lon,lat,tf,100000);k2=source_coarse_cell(lon,lat,tf,200000)
            if k1 is None or k2 is None:continue
            qa[g]+=1;o=int(r["observer_id"]) if r["observer_id"] else None
            for w,k in ((100000,k1),(200000,k2)):
                if k not in cand[w][g]:continue
                add(full[w][g][k],t,o)
                if t not in allowed[g]:continue
                add(target[w][g][k],t,o)
                for ax,i in axes[g].get(t,()):add(guild[w][g][ax][i][k],t,o)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or qa["rows"]!=288721938 or raw.hash.hexdigest()!=src["obs_sha256"] or any(qa[g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Exact original source replay failed")
    out={}
    for w in WIDTHS:
        out[str(w)]={}
        for g in GROUPS:
            if set(full[w][g])!=set(cand[w][g]):raise ValueError("Full endpoint support changed")
            for k,v in full[w][g].items():
                if (full_counts(v),counts(v))!=cand[w][g][k]:raise ValueError("Full/half count mismatch")
                t=target[w][g].get(k,bucket())
                if tuple(tc[w][g][k]["target_full"]["counts"])!=full_counts(t) or tuple(tuple(x["counts"]) for x in tc[w][g][k]["target_halves"])!=counts(t):raise ValueError("Target taxon census changed")
                for ax in AXES[g]:
                    gg=[guild[w][g][ax][i].get(k,bucket()) for i in (0,1)]
                    for h in (0,1):
                        z=tc[w][g][k]["target_halves"][h]["traits"][ax]
                        n=sum(counts(x)[h][2] for x in gg)
                        if n!=z["known"]:raise ValueError("Exact directional known count changed")
            pairs=[]
            for p in src["pairs"][str(w)][g]:
                if "split_similarity" not in p:continue
                x,y=key(p);v,u=full[w][g][x],full[w][g][y]
                if overlap(v[1],u[1])!=p["species_full"] or any(overlap(v[4][i],u[4][j])!=p["split_similarity"][f"{i}{j}"] for i in (0,1) for j in (0,1)):raise ValueError("Exact source pair similarity changed")
                z={"a":x,"b":y,"axes":{}}
                for ax in AXES[g]:
                    aa=[guild[w][g][ax][i].get(x,bucket()) for i in (0,1)]
                    bb=[guild[w][g][ax][i].get(y,bucket()) for i in (0,1)]
                    da=[t["traits"][ax]["denominator"] for t in tc[w][g][x]["target_halves"]];db=[t["traits"][ax]["denominator"] for t in tc[w][g][y]["target_halves"]]
                    q=pack_pair(aa,bb,da,db)
                    if q is not None:z["axes"][ax]=q
                if z["axes"]:pairs.append(z)
            out[str(w)][g]=pairs
            print(json.dumps({"guild_pair_source":{"width":w,"taxon":g,"n_pairs":len(pairs),"support":{ax:{str(c):sum(ax in p["axes"] and p["axes"][ax]["min_half_coverage"]>=c for p in pairs) for c in (.5,.25)} for ax in AXES[g]}}}),flush=True)
    result={"version":"railstatus-disturbance-guild-pairs-v1","obs_sha256":src["obs_sha256"],"taxa_sha256":tsha,"split_seed":SALT.decode(),
        "exact_full_target_guild_count_and_pair_replay":True,"trait_cells_sha256":hashlib.sha256(a.trait_cells.read_bytes()).hexdigest(),
        "guild_gate":GATE,"axes":AXES,"pairs":out}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
