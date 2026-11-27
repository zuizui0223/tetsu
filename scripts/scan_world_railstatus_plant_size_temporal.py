#!/usr/bin/env python3
"""Annual height validation on frozen 200km plant geography."""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter,defaultdict
from scan_world_railstatus_verified_plant_size import (
    taxonomy,source_client,taxonomy_metadata,load_life,fetch_fig,parse_scaled,add_size,
    LIFE_SHA,SIZE_SHA,AXES,OLD,check_old,extract,bucket,full_counts,split_counts,
    summarize,family_refs,SALT)
from scan_world_railstatus_cross_taxa import (
    open_original_source,source_coarse_cell,OBS_ETAG,OBS_BYTES,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS)
from analyze_world_railstatus_plant_bird_disjoint import YEARS
W=200000

def read(p):return json.loads(p.read_text())
def frozen_candidates(src,annual):
    c={"2025":extract(src["pairs"][str(W)]["plants"])}
    for y in YEARS[:2]:
        c[y]={k:v for k,v in extract(annual["pairs"][y]["plants"]).items() if k in c["2025"]}
    return c

def scan(client,canonical,allowed,c):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    full={y:defaultdict(bucket) for y in YEARS};target={y:defaultdict(bucket) for y in YEARS};qa=Counter()
    from scan_world_railstatus_observer_split_beta import add
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=set(rr.fieldnames or []):raise ValueError("Observation schema changed")
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"height_annual_scan":dict(qa)}),flush=True)
            y=r["observed_on"][:4]
            if y not in YEARS or not r["observed_on"].startswith(y+"-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0]!="plants" or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k not in c[y] or source_coarse_cell(lon,lat,tf,100000) is None:continue
            t=hit[1];o=int(r["observer_id"]) if r["observer_id"] else None
            add(full[y][k],t,o)
            if t in allowed:add(target[y][k],t,o)
            qa[y+"_plant_records"]+=1
    finally:stream.close()
    if qa["rows"]!=288721938 or raw.bytes!=OBS_BYTES:raise ValueError("Incomplete original source")
    return full,target,dict(qa),raw.hash.hexdigest()

def compare(new,old):
    for x,y in zip([new["target_full"],*new["target_halves"]],[old["target_full"],*old["target_halves"]]):
        if x["counts"]!=y["counts"]:raise ValueError("Prior target counts differ")
        for f in (*OLD,*AXES):
            for k,v in x["traits"][f].items():
                q=y["traits"][f][k]
                if v is None or q is None:
                    if v!=q:raise ValueError("Prior trait missingness changed")
                elif abs(v-q)>1e-8:raise ValueError("Prior trait moments changed")

def main():
    p=argparse.ArgumentParser()
    for k in ("split","annual","size","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src,ann,old=map(read,(a.split,a.annual,a.size))
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or ann["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or ann["width_m"]!=W or ann["qa"]["rows"]!=288721938:raise ValueError("Wrong annual source")
    if old["version"]!="railstatus-verified-scaled-plant-size-cells-v1" or not old["exact_full_split_target_and_life_replay"] or old["source_SHA256"]!={"LifeForm":LIFE_SHA,"size_archive":SIZE_SHA} or any(len({v[k] for v in (src,ann,old)})!=1 for k in ("obs_sha256","taxa_sha256","split_seed")) or src["split_seed"]!=SALT.decode():raise ValueError("Wrong height source")
    c=frozen_candidates(src,ann)
    om={tuple(r["key"]):r for r in old["grids"][str(W)]["plants"]}
    if {k:(tuple(v["full"]),tuple(tuple(x) for x in v["halves"])) for k,v in om.items()}!=c["2025"]:raise ValueError("Frozen 2025 cell support drift")
    print(json.dumps({"annual_height_candidates":{y:len(v) for y,v in c.items()}}),flush=True)
    client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("Taxonomy changed")
    names,allowed,fam,_=taxonomy_metadata(client,canonical,sha)
    life=load_life(names["plants"],allowed["plants"])
    rows,audit=parse_scaled(fetch_fig(25600650,45643398,SIZE_SHA));audit.update(add_size(life,rows,names["plants"],allowed["plants"]));audit.update(taxonomy_sha256=sha,LifeForm_sha256=LIFE_SHA)
    if audit!=old["source"]:raise ValueError("Prior size source not reproduced")
    ref=family_refs(life,fam,"plants");del names,rows
    full,target,qa,osha=scan(client,canonical,allowed["plants"],c)
    if osha!=src["obs_sha256"]:raise ValueError("Original observation hash changed")
    grids={}
    for y in YEARS:
        if set(full[y])!=set(c[y]):raise ValueError("Annual candidate cells changed")
        rr=[]
        for k,v in sorted(full[y].items()):
            f=full_counts(v);h=tuple(tuple(x) for x in split_counts(v))
            if (f,h)!=c[y][k]:raise ValueError("Original annual count replay failed")
            t=target[y].get(k,bucket());ss=[t[1],*t[4]];cnt=[full_counts(t),*split_counts(t)]
            qq=[{"counts":q,"traits":summarize(s,life,fam,ref,"plants")} for s,q in zip(ss,cnt)]
            z={"key":k,"full":f,"halves":h,"target_full":qq[0],"target_halves":qq[1:]}
            if y=="2025":compare(z,om[k])
            rr.append(z)
        grids[y]=rr
        print(json.dumps({"annual_height_cells":y,"n":len(rr)}),flush=True)
    out={"version":"railstatus-annual-verified-plant-size-v1","width":W,"source":audit,
        "obs_sha256":osha,"taxa_sha256":sha,"split_seed":SALT.decode(),"exact_annual_counts_and_2025_moments_replay":True,
        "source_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","annual","size")},"qa":qa,"grids":grids}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
