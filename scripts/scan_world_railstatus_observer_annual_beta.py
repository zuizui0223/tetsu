#!/usr/bin/env python3
"""Frozen 2023/2024 retest on the 2025 observer-split bird/plant geography.
Only 200 km coarse pair summaries are exported; no individual IDs.
"""
import argparse,json,math,time
from collections import Counter,defaultdict
from pathlib import Path
from scan_world_railstatus_cross_taxa import (
    taxonomy,source_client,open_original_source,source_coarse_cell,
    OBS_ETAG,TAXA_ETAG,OBS_BYTES,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS)
from scan_world_railstatus_observer_split_beta import bucket,add,pair_rows,MIN_FULL,MIN_SPLIT,SALT
YEARS=("2023","2024")
GROUPS=("plants","birds")
W=200000

def pair_key(p):
    return tuple(p["a"]),tuple(p["b"])

def candidate_support(source):
    pairs={g:{pair_key(p) for p in source["pairs"][str(W)][g] if "split_similarity" in p} for g in GROUPS}
    cells={g:{k for pair in pairs[g] for k in pair} for g in GROUPS}
    return pairs,cells

def scan(client,canonical,cells):
    from pyproj import Transformer
    trans=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    fields=set(rows.fieldnames or [])
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=fields:raise ValueError("Source schema drift")
    if {"geoprivacy","coordinates_obscured","taxon_geoprivacy"}&fields:raise ValueError("Geoprivacy schema changed")
    data={y:{g:defaultdict(bucket) for g in GROUPS} for y in YEARS}
    qa=Counter();t0=time.monotonic()
    try:
        for r in rows:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"qa":dict(qa),"seconds":round(time.monotonic()-t0,1)}),flush=True)
            date=r["observed_on"];y=date[:4]
            if y not in YEARS or not date.startswith(y+"-"):continue
            if r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]))
            except (TypeError,ValueError):continue
            if hit is None or hit[0] not in GROUPS:continue
            g,species=hit
            try:acc=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if not math.isfinite(acc) or not 0<=acc<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,trans,W)
            if k is None:continue
            qa[y+"_"+g+"_all_grid"]+=1
            if k not in cells[g]:continue
            o=int(r["observer_id"]) if r["observer_id"] else None
            qa[y+"_"+g+"_fixed_cells"]+=1
            if o is None:qa[y+"_"+g+"_missing_observer"]+=1
            add(data[y][g][k],species,o)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or qa["rows"]!=288721938:raise ValueError("Exact original source incomplete")
    return data,dict(qa),raw.hash.hexdigest()

def main():
    p=argparse.ArgumentParser()
    for n in ("split","prior","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();src=json.loads(a.split.read_text());prior=json.loads(a.prior.read_text())
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"]:raise ValueError("Wrong source split")
    if src["split_seed"]!=SALT.decode() or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or src["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Frozen split protocol changed")
    if prior["version"]!="railstatus-cross-taxa-v1" or prior["obs_etag"]!=OBS_ETAG or prior["taxa_etag"]!=TAXA_ETAG:raise ValueError("Wrong prior source")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Prior source mismatch")
    pairs,cells=candidate_support(src)
    client=source_client();canonical,counts,tsha=taxonomy(client)
    if tsha!=src["taxa_sha256"] or counts!=prior["active_species_counts"]:raise ValueError("Taxonomy mismatch")
    data,qa,osha=scan(client,canonical,cells)
    if osha!=src["obs_sha256"]:raise ValueError("Observation SHA mismatch")
    result={}
    for y in YEARS:
        result[y]={}
        for g in GROUPS:
            if qa[y+"_"+g+"_fixed_cells"]<=0:raise ValueError("No annual original observations")
            rr=[r for r in pair_rows(data[y][g]) if pair_key(r) in pairs[g] and "split_similarity" in r]
            result[y][g]=rr
            print(json.dumps({"year":y,"group":g,"fixed_2025_candidate_pairs":len(pairs[g]),"annual_eligible_pairs":len(rr)}),flush=True)
    out={"version":"railstatus-fixed-2025-observer-split-annual-v1","obs_sha256":osha,"taxa_sha256":tsha,
         "qa":qa,"width_m":W,"split_seed":SALT.decode(),"years":YEARS,"groups":GROUPS,
         "candidate_pairs":{g:len(pairs[g]) for g in GROUPS},"pairs":result,
         "limits":["only 2025 split-supported candidate pair geography","same platform, observer and species may repeat over years","2026 railway status not historical operations","no raw IDs or coordinates"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({k:v for k,v in out.items() if k!="pairs"}),flush=True)
if __name__=="__main__":main()
