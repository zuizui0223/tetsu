#!/usr/bin/env python3
"""Single exact replay of v369 frozen 2025 photo-linked iNat data for four taxa.
Publish only grid counts, never observation IDs, taxon IDs or coordinates.
"""
import argparse, json, math, time
from collections import Counter, defaultdict
from pathlib import Path
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    source_client, open_original_source, source_ancestry_ids, source_coarse_cell,
    OBS_ETAG, TAXA_ETAG, OBS_BYTES, TAXA_BYTES, YEAR, MAX_UNCERT_M,
    SOURCE_REGION_M, PROGRESS_EVERY_SOURCE_ROWS)
ROOTS={"plants":(47126,"plantae"),"birds":(3,"aves"),
       "insects":(47158,"insecta"),"fungi":(47170,"fungi")}
INFRA={"subspecies","variety","form","subvariety"}

def taxonomy(client):
    rows, raw, stream = open_original_source(client,"taxa.csv.gz",TAXA_ETAG,TAXA_BYTES)
    taxa={}
    try:
        if not {"taxon_id","ancestry","rank","active","name"}<=set(rows.fieldnames or []):
            raise ValueError("Taxonomy columns drifted")
        for row in rows:
            tid=int(row["taxon_id"])
            if tid in taxa: raise ValueError("Duplicate taxon")
            taxa[tid]=(row["rank"].strip().lower(),
                       row["active"].strip().lower() in ("t","true","1"),
                       tuple(map(int,source_ancestry_ids(row["ancestry"]))),
                       row["name"].strip().lower() if tid in {r[0] for r in ROOTS.values()} else "")
    finally: stream.close()
    if raw.bytes!=TAXA_BYTES: raise ValueError("Incomplete taxonomy")
    for group,(root,name) in ROOTS.items():
        if root not in taxa or taxa[root][3]!=name:
            raise ValueError("Unexpected taxonomic root: "+group)
    canonical={}
    roots={v[0]:k for k,v in ROOTS.items()}
    for tid,(rank,active,anc,_) in taxa.items():
        if rank!="species" or not active: continue
        groups=[roots[x] for x in anc if x in roots]
        if len(groups)>1: raise ValueError("Overlapping roots")
        if groups: canonical[tid]=(groups[0],tid)
    species_counts=Counter(g for g,s in canonical.values())
    if any(species_counts[g]<1000 for g in ROOTS): raise ValueError("Taxonomic coverage implausibly low")
    for tid,(rank,active,anc,_) in taxa.items():
        if not active or rank not in INFRA: continue
        hit=next((canonical[x] for x in reversed(anc)
                  if x in canonical and taxa[x][0]=="species"),None)
        if hit is not None: canonical[tid]=hit
    return canonical,dict(species_counts),raw.hash.hexdigest()

def bucket():
    return [0,set(),set()]

def add(cells,group,key,species,observer):
    entry=cells[group][key]
    entry[0]+=1
    entry[1].add(species)
    if observer is not None: entry[2].add(observer)

def serialize(cells):
    out={}
    for width,groups in cells.items():
        out[str(width)]={}
        for group,cmap in groups.items():
            out[str(width)][group]=[
                {"source_equal_area_grid_ix":i,"source_equal_area_grid_iy":j,
                 "observations":v[0],"species":len(v[1]),"observers":len(v[2])}
                for (i,j),v in sorted(cmap.items())]
    return out

def scan(client,canonical):
    from pyproj import Transformer
    trans=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    fields=set(rows.fieldnames or [])
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=fields:
        raise ValueError("Observations schema drift")
    if {"geoprivacy","coordinates_obscured","taxon_geoprivacy"}&fields:
        raise ValueError("Geoprivacy schema changed")
    cells={m:{g:defaultdict(bucket) for g in ROOTS} for m in SOURCE_REGION_M}
    qa=Counter()
    t0=time.monotonic()
    try:
        for row in rows:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:
                print(json.dumps({"qa":dict(qa),"seconds":round(time.monotonic()-t0,1)}),flush=True)
            if not row["observed_on"].startswith(YEAR+"-"): continue
            if row["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research": continue
            try: hit=canonical.get(int(row["taxon_id"]))
            except (TypeError,ValueError): continue
            if hit is None: continue
            group,species=hit
            try:
                acc=float(row["positional_accuracy"]); lon=float(row["longitude"]); lat=float(row["latitude"])
            except (TypeError,ValueError): continue
            if not math.isfinite(acc) or not 0<=acc<=MAX_UNCERT_M: continue
            small=source_coarse_cell(lon,lat,trans,100000)
            large=source_coarse_cell(lon,lat,trans,200000)
            if small is None or large is None: continue
            observer=int(row["observer_id"]) if row["observer_id"] else None
            qa[group]+=1
            add(cells[100000],group,small,species,observer)
            add(cells[200000],group,large,species,observer)
    finally: stream.close()
    if raw.bytes!=OBS_BYTES or qa["rows"]<50000000: raise ValueError("Incomplete frozen original observations")
    out=serialize(cells)
    for g in ROOTS:
        if qa[g]<=0 or any(sum(z["observations"] for z in out[str(m)][g])!=qa[g] for m in SOURCE_REGION_M):
            raise ValueError("Cross-resolution source count mismatch "+g)
    return out,dict(qa),raw.hash.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    client=source_client()
    canonical,counts,tsha=taxonomy(client)
    grids,qa,osha=scan(client,canonical)
    result={"version":"railstatus-cross-taxa-v1","year":YEAR,
            "obs_etag":OBS_ETAG,"taxa_etag":TAXA_ETAG,"obs_sha256":osha,"taxa_sha256":tsha,
            "roots":ROOTS,"active_species_counts":counts,"qa":qa,"grids":grids,
            "scope":"2025 public photo-linked research-grade observations; same v369 filters; whole landscapes, not presence/absence"}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,allow_nan=False),encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k!="grids"}),flush=True)
if __name__=="__main__":main()
