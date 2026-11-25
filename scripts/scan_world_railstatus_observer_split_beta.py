#!/usr/bin/env python3
"""Exact source replay: adjacent landscape species similarity for four taxa
and a globally disjoint observer split. No individual observer/taxon IDs or
coordinates leave the runner. Whole-landscape observation data only.
"""
import argparse,hashlib,json,math,time
from collections import Counter,defaultdict
from functools import lru_cache
from pathlib import Path
from scan_world_railstatus_cross_taxa import (
    taxonomy,ROOTS,source_client,open_original_source,source_coarse_cell,
    OBS_ETAG,TAXA_ETAG,OBS_BYTES,YEAR,MAX_UNCERT_M,SOURCE_REGION_M,PROGRESS_EVERY_SOURCE_ROWS)

GRID="source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"
PAIRS="source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"
SALT=b"tetsu-railstatus-observer-split-v1"
MIN_FULL=(50,5,20)
MIN_SPLIT=(25,5,10)

@lru_cache(maxsize=262144)
def observer_side(observer):
    # Same observer always has the same side across cells, taxa and scales.
    return hashlib.blake2b(str(observer).encode(),key=SALT,digest_size=8).digest()[0]&1

def bucket():
    return [0,set(),set(),[0,0],[set(),set()]]

def add(entry,species,observer):
    entry[0]+=1;entry[1].add(species)
    if observer is None:return
    entry[2].add(observer)
    side=observer_side(observer);entry[3][side]+=1;entry[4][side].add(species)

def full_counts(entry):
    return (entry[0],len(entry[2]),len(entry[1]))

def split_counts(entry):
    b=[0,0]
    for o in entry[2]:b[observer_side(o)]+=1
    return [(entry[3][i],b[i],len(entry[4][i])) for i in range(2)]

def eligible(counts,gate):
    return all(x>=y for x,y in zip(counts,gate))

def overlap(a,b):
    sh=len(a&b);union=len(a)+len(b)-sh
    return {"shared":sh,"union":union,"jaccard":sh/union if union else None,
            "simpson_similarity":sh/min(len(a),len(b)) if min(len(a),len(b)) else None}

def pair_rows(cells):
    valid={k:v for k,v in cells.items() if eligible(full_counts(v),MIN_FULL)}
    counts={k:split_counts(v) for k,v in valid.items()}
    out=[]
    for (i,j),a in sorted(valid.items()):
        for bkey in ((i+1,j),(i,j+1)):
            b=valid.get(bkey)
            if b is None:continue
            ac,bc=counts[i,j],counts[bkey]
            r={"a":[i,j],"b":list(bkey),"full_a":full_counts(a),"full_b":full_counts(b),
                "species_full":overlap(a[1],b[1]),"observer_full":overlap(a[2],b[2])}
            # Publish split statistics only if both observer halves pass fixed
            # gates in both endpoints. All four combinations have same support.
            if all(eligible(c,MIN_SPLIT) for c in (*ac,*bc)):
                r["split_counts_a"]=ac;r["split_counts_b"]=bc
                r["split_similarity"]={f"{x}{y}":overlap(a[4][x],b[4][y]) for x in range(2) for y in range(2)}
                r["same_side_observer_overlap"]={str(s):overlap(
                    {o for o in a[2] if observer_side(o)==s},
                    {o for o in b[2] if observer_side(o)==s}) for s in range(2)}
            out.append(r)
    return out

def scan(client,canonical):
    from pyproj import Transformer
    trans=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    fields=set(rows.fieldnames or [])
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=fields:raise ValueError("Source schema drift")
    if {"geoprivacy","coordinates_obscured","taxon_geoprivacy"}&fields:raise ValueError("Geoprivacy schema changed")
    cells={w:{g:defaultdict(bucket) for g in ROOTS} for w in SOURCE_REGION_M}
    qa=Counter();t0=time.monotonic()
    try:
        for r in rows:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"qa":dict(qa),"seconds":round(time.monotonic()-t0,1)}),flush=True)
            if not r["observed_on"].startswith(YEAR+"-"):continue
            if r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]))
            except (TypeError,ValueError):continue
            if hit is None:continue
            g,species=hit
            try:acc=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if not math.isfinite(acc) or not 0<=acc<=MAX_UNCERT_M:continue
            small=source_coarse_cell(lon,lat,trans,100000);large=source_coarse_cell(lon,lat,trans,200000)
            if small is None or large is None:continue
            o=int(r["observer_id"]) if r["observer_id"] else None
            qa[g]+=1
            if o is None:qa[g+"_missing_observer"]+=1
            add(cells[100000][g][small],species,o);add(cells[200000][g][large],species,o)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or qa["rows"]!=288721938:raise ValueError("Exact source row/byte mismatch")
    return cells,dict(qa),raw.hash.hexdigest()

def validate_counts(cmap,prior):
    new={k:full_counts(v) for k,v in cmap.items()}
    old={(z["source_equal_area_grid_ix"],z["source_equal_area_grid_iy"]):
         (z["observations"],z["observers"],z["species"]) for z in prior}
    if new!=old:raise ValueError("Exact original four-taxon grid replay mismatch")

def validate_plant_pairs(pairs,prior):
    new={(tuple(r["a"]),tuple(r["b"])):(r["species_full"]["shared"],r["species_full"]["union"]) for r in pairs}
    old={(tuple(r["source_cell_A_original_grid_ix_iy"]),tuple(r["source_cell_B_original_grid_ix_iy"])):
         (r["source_shared_plants_species_n"],r["source_union_plants_species_n"]) for r in prior}
    if new!=old:raise ValueError("Exact original plant Jaccard pair replay mismatch")

def main():
    p=argparse.ArgumentParser()
    for n in ("prior","plant","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();prior=json.loads(a.prior.read_text());plant=json.loads(a.plant.read_text())
    if prior["version"]!="railstatus-cross-taxa-v1" or prior["obs_etag"]!=OBS_ETAG or prior["taxa_etag"]!=TAXA_ETAG:raise ValueError("Wrong prior")
    client=source_client();canonical,counts,tsha=taxonomy(client)
    if tsha!=prior["taxa_sha256"] or tsha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Taxonomy replay mismatch")
    if counts!=prior["active_species_counts"]:raise ValueError("Taxonomy counts changed")
    cells,qa,osha=scan(client,canonical)
    if osha!=prior["obs_sha256"] or osha!=plant["source_summary"]["source_iNat_original_2026_10_27_observations_gzip_SHA256_full"]:raise ValueError("Observations SHA changed")
    if any(qa[g]!=prior["qa"][g] for g in ROOTS):raise ValueError("Source group count mismatch")
    result={}
    for w in SOURCE_REGION_M:
        result[str(w)]={}
        for g in ROOTS:
            cmap=cells[w][g];validate_counts(cmap,prior["grids"][str(w)][g])
            prs=pair_rows(cmap)
            if g=="plants":validate_plant_pairs(prs,plant["source_grids"][str(w)][PAIRS])
            result[str(w)][g]=prs
            print(json.dumps({"width_m":w,"taxon":g,"n_pairs":len(prs),"n_split_eligible":sum("split_similarity" in r for r in prs)}),flush=True)
    out={"version":"railstatus-four-taxon-observer-split-beta-v1","qa":qa,"obs_sha256":osha,"taxa_sha256":tsha,
         "full_gate_observations_observers_species":MIN_FULL,"split_gate_observations_observers_species":MIN_SPLIT,
         "split_seed":SALT.decode(),"exact_grid_replay":True,"exact_plant_beta_replay":True,"pairs":result,
         "limits":["same observer may observe both cells in full/same-side comparisons","cross-side comparisons use disjoint observer IDs but not independent access patterns or detection","split support may favor well-sampled regions","whole landscapes, not trackside flora","no individual IDs or exact coordinates exported"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({k:v for k,v in out.items() if k!="pairs"}),flush=True)
if __name__=="__main__":main()
