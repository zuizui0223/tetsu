#!/usr/bin/env python3
"""One frozen original-source scan for insect counts in 2023–25.
Restrict to prior 2025 split-supported insect x plant-bird diet cell geography.
Exact 2025 full and half counts replay; no taxon, observer or occurrence IDs exported.
"""
import argparse,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict
from scan_world_railstatus_cross_taxa import (
    taxonomy,source_client,open_original_source,source_coarse_cell,OBS_ETAG,
    OBS_BYTES,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS)
from scan_world_railstatus_observer_split_beta import (
    bucket,add,full_counts,split_counts,eligible,MIN_FULL,MIN_SPLIT,SALT)
from analyze_world_railstatus_plant_bird_disjoint import extract,YEARS
W=200000

def candidate_cells(src,diet):
    old=extract(src["pairs"][str(W)]["insects"])
    dc={y:{tuple(r["key"]) for r in diet["grids"][str(W)][y]} for y in YEARS}
    fixed=set(old)&dc["2025"]
    if not fixed:raise ValueError("Empty frozen insect candidate geography")
    return {y:fixed&dc[y] for y in YEARS},{k:old[k] for k in fixed}

def scan(client,canonical,cand):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=set(rr.fieldnames or []):raise ValueError("Source schema changed")
    data={y:defaultdict(bucket) for y in YEARS};qa=Counter()
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"scan":dict(qa)}),flush=True)
            y=r["observed_on"][:4]
            if y not in YEARS or not r["observed_on"].startswith(y+"-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0]!="insects" or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k not in cand[y]:continue
            add(data[y][k],hit[1],int(r["observer_id"]) if r["observer_id"] else None);qa[y+"_records"]+=1
    finally:stream.close()
    if qa["rows"]!=288721938 or raw.bytes!=OBS_BYTES:raise ValueError("Incomplete frozen source")
    return data,dict(qa),raw.hash.hexdigest()

def main():
    p=argparse.ArgumentParser()
    for k in ("split","diet","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src=json.loads(a.split.read_text());diet=json.loads(a.diet.read_text())
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_seed"]!=SALT.decode() or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or src["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Wrong split")
    if diet["version"]!="railstatus-avonet-trophic-disjoint-cell-counts-v1" or not diet["exact_full_and_split_cell_replay"] or any(diet[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong diet source")
    cand,old=candidate_cells(src,diet)
    print(json.dumps({"fixed_candidates":{y:len(v) for y,v in cand.items()}}),flush=True)
    client=source_client();canonical,_,tsha=taxonomy(client)
    if tsha!=src["taxa_sha256"]:raise ValueError("Frozen taxonomy changed")
    data,qa,osha=scan(client,canonical,cand)
    if osha!=src["obs_sha256"]:raise ValueError("Observation hash changed")
    if set(data["2025"])!=set(old):raise ValueError("2025 insect source support changed")
    result={}
    for y in YEARS:
        rr=[]
        for k,v in sorted(data[y].items()):
            f=full_counts(v);h=tuple(tuple(x) for x in split_counts(v))
            if y=="2025" and (f,h)!=old[k]:raise ValueError("2025 exact insect full/half replay failed")
            if eligible(f,MIN_FULL) and all(eligible(x,MIN_SPLIT) for x in h):
                rr.append({"key":list(k),"full":f,"halves":h})
        result[y]=rr
        print(json.dumps({"annual_insect_year":y,"candidate_cells":len(cand[y]),"observed_cells":len(data[y]),"supported_cells":len(rr)}),flush=True)
    if len(result["2025"])!=len(old):raise ValueError("2025 split eligibility changed")
    out={"version":"railstatus-annual-insect-disjoint-cells-v1","width_m":W,"years":list(YEARS),"obs_sha256":osha,"taxa_sha256":tsha,"split_seed":SALT.decode(),
         "exact_2025_full_and_half_replay":True,"source_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","diet")},"source_qa":qa,"cells":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
