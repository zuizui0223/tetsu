#!/usr/bin/env python3
"""v405: fixed 2018–22 references, annual 2023–25 temporal retest on common geography."""
import argparse,json,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from pyproj import Transformer
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
 source_client,open_original_source,taxonomy_source,source_coarse_cell,OBS_ETAG,OBS_BYTES,MAX_UNCERT_M,OBS_QUALITY)
from analyze_global_multiyear_naturalized_v397 import entry
from analyze_global_rail_documented_naturalized_v379 import grid_counts,grid_key
from analyze_global_crossfit_climate_niche_v399 import (
 REFERENCE_YEARS,encode_refs,W,MIN_KNOWN,MIN_COVERAGE,OUTCOMES,BOOT,TARGET,group_indices,slopes)
from analyze_global_niche_shared_taxa_v401 import traits,metrics,NICHE_FIELDS
from analyze_global_niche_common_support_v400 import prepare,MODES,high,RESP,COVERAGE_FIELDS,INTERACTION_FIELDS,COLD
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import source_climate_lookup,block_weights
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
from analyze_global_turnover_climate_block_v387 import load_extra_climate
YEARS=("2023","2024","2025")
SEED=40520261125
def scan_years(client,canonical,plant):
    ref=defaultdict(set);tests={y:defaultdict(entry) for y in YEARS};test=tests['2025'];counts=Counter();total=0
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    try:
        required={"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}
        if not required<=set(rows.fieldnames or []):raise ValueError("Original source schema drift")
        for r in rows:
            total+=1
            if total%10000000==0:print(json.dumps({"v399_original_rows":total,"eligible_records":counts}),flush=True)
            y=str(r["observed_on"])[:4]
            if y not in (*REFERENCE_YEARS,*YEARS) or not str(r["observed_on"]).startswith(y+"-"):continue
            if str(r["quality_grade"]).strip().lower().replace(" ","_").replace("-","_")!=OBS_QUALITY:continue
            try:
                tid=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"])
                lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if tid is None or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k is None:continue
            counts[y]+=1
            if y in YEARS:
                a=tests[y][k];a[0]+=1;a[1].add(tid)
                if r["observer_id"]:a[2].add(int(r["observer_id"]))
            else:ref[tid].add(k)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or raw.hash.hexdigest()!=plant["source_summary"]["source_iNat_original_2026_10_27_observations_gzip_SHA256_full"] or total!=288721938 or counts["2025"]!=7827936:
        raise ValueError("Original snapshot full replay mismatch")
    old={grid_key(r):(r["source_year_2025_photo_linked_plants_observation_count"],r["source_2025_plant_species_richness"],r["source_2025_distinct_observers_with_photo_linked_plant_records"]) for r in grid_counts(plant,W)}
    if old!={k:(a[0],len(a[1]),len(a[2])) for k,a in test.items()}:raise ValueError("2025 original exact cell replay failed")
    return ref,tests,{"all_source_rows":total,"eligible_by_year":dict(counts),"source_observation_SHA256":raw.hash.hexdigest(),
                    "exact_2025_all_cell_effort_richness_replay":True}

def annual_row(r,entry,known,mode):
    n,spp,observers=entry
    if n<50 or len(observers)<5 or len(spp)<20:return None
    k=len(known);coverage=k/len(spp)
    if k<MIN_KNOWN or coverage<MIN_COVERAGE:return None
    a=dict(r,observations=n,observers=len(observers),species=len(spp),
        log_observations=math.log1p(n),log_observers=math.log1p(len(observers)),log_species=math.log1p(len(spp)),
        known_niche_species=k,trait_coverage=coverage)
    j=MODES.index(mode)
    return prepare(dict(a,**metrics(np.stack([v[j][1] for v in known],axis=1),a)))
def rebuild(records,tests,encoded):
    bg=defaultdict(list)
    for r in records:
        a=prepare(r["variants"][0]["shared_rarefied"])
        bg[(a["continent"],a["block"])].append(r)
    out={};cont_cache={}
    for (cont,b),rs in sorted(bg.items()):
        caches=[{},cont_cache.setdefault(cont,{})]
        for rec in rs:
            k=tuple(rec["key"]);old=[prepare(v["shared_rarefied"]) for v in rec["variants"]]
            annual={}
            for y in YEARS:
                a=tests[y].get(k)
                if a is None:annual[y]=None;continue
                known=[]
                for tid in sorted(a[1]):
                    vv=[]
                    for j,m in enumerate(MODES):
                        if tid not in caches[j]:caches[j][tid]=traits(tid,old[j],m,encoded)
                        vv.append(caches[j][tid])
                    if all(v is not None for v in vv):known.append(vv)
                rr=[annual_row(old[j],a,known,m) for j,m in enumerate(MODES)]
                if (rr[0] is None)!=(rr[1] is None):raise ValueError("Identical annual shared-taxon gates disagree")
                annual[y]=rr if rr[0] is not None else None
                if y=="2025":
                    if annual[y] is None:raise ValueError("Original 2025 shared support disappeared")
                    for j in range(2):
                        for f in ("known_niche_species","trait_coverage","species","observations","observers",*NICHE_FIELDS):
                            if abs(rr[j][f]-old[j][f])>1e-9:raise ValueError("Original v401 2025 rarefied replay mismatch: "+f)
            out[k]=annual
    return out
def stat(X,g,w,fields):
    c=slopes(X,g,w,fields,(TARGET,COLD))
    if c is None:return None
    a,b=np.asarray(c[TARGET]),np.asarray(c[COLD])
    return np.stack((a,b,a+b,a-b))
STAT=("rail_at_15C","rail_cold","rail_at_5C","rail_at_25C")
def paired(rows,fields,seed):
    g,b=group_indices(rows["2025"],"block")
    n=len(rows["2025"]);out={"n_cells":n,"n_blocks":len(b),"cells_by_continent":dict(Counter(r["continent"] for r in rows["2025"]))}
    if n<100 or len(b)<20 or n-len(g)-len(fields)<30:return dict(out,status="insufficient common-year support")
    X={y:np.asarray([[r[k] for k in fields+RESP] for r in rows[y]],float) for y in YEARS}
    for y in YEARS:
        if [r["key"] for r in rows[y]]!=[r["key"] for r in rows["2025"]]:raise ValueError("Year keys changed")
    w=block_weights(rows["2025"])
    def est(w):
        v={y:stat(X[y],g,w,fields) for y in YEARS}
        if any(x is None for x in v.values()):return None
        v["2023_2024_mean"]=(v["2023"]+v["2024"])/2.
        v["2025_minus_earlier_mean"]=v["2025"]-v["2023_2024_mean"]
        return v
    c=est(w)
    if c is None:return dict(out,status="insufficient residual exposure variation")
    bc=defaultdict(list)
    for k in sorted(b):bc[k[0]].append(k)
    rng=np.random.default_rng(seed);draws={y:[] for y in c}
    for _ in range(BOOT):
        bw=np.zeros(n)
        for bs in bc.values():
            m=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for k,v in zip(bs,m):bw[b[k]]=v/len(b[k])
        z=est(bw)
        if z is not None:
            for y,v in z.items():draws[y].append(v)
    out["status"]="estimated on identical three-year geography"
    out["coefficients"]={y:{s:dict(zip(RESP,map(float,a[i]))) for i,s in enumerate(STAT)} for y,a in c.items()}
    out["block_percentile_95_descriptive"]={y:{s:{t:{
        "valid":len(v),"low":float(np.quantile(np.asarray(v)[:,i,j],.025)) if len(v)>=.95*BOOT else None,
        "high":float(np.quantile(np.asarray(v)[:,i,j],.975)) if len(v)>=.95*BOOT else None}
        for j,t in enumerate(RESP)} for i,s in enumerate(STAT)} for y,v in draws.items()}
    return out
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","multiyear","plants","ne","climate","out","aggregates"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    src,orig,multi,plant=map(read,(a.source,a.niche,a.multiyear,a.plants));q=src["source_QA"]
    if src["version"]!="v401-original-shared-taxon-rarefaction-coarse" or orig["version"]!="v399-crossfit-niche-coarse-aggregate" or multi["version"]!="v397-original-2023-2025-coarse-aggregates":raise ValueError("Wrong source artifacts")
    if q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original source replay")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong map or climate")
    records=[r for r in src["cells"] if r["variants"]]
    if len(records)!=1323 or len({tuple(r["key"]) for r in records})!=1323:raise ValueError("Original v401 denominator changed")
    land80={tuple(r["key"]) for p in orig["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    client=source_client();canonical,_,_,_,sha=taxonomy_source(client)
    if sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Original taxa source changed")
    refs,tests,newqa=scan_years(client,canonical,plant)
    if any(newqa["eligible_by_year"][y]!=n for y,n in q["eligible_by_year"].items()):raise ValueError("v399 reference yearly counts changed")
    if {y:newqa["eligible_by_year"][y] for y in YEARS}!=multi["source_replay"]["annual_eligible_original_plant_records"]:raise ValueError("v397 original annual totals changed")
    keys=set(tests["2025"])
    for ks in refs.values():keys.update(ks)
    land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);geo,_=source_climate_lookup(land,raster,scale,keys,W)
    finally:raster.close();mem.close()
    extra,_,_=load_extra_climate(a.climate,set(geo),W)
    encoded=encode_refs(refs,geo,extra);del refs
    oq=orig["source_QA"]
    if len(encoded[0])!=oq["reference_taxa_with_at_least_20_climate_cells"] or len(geo)!=oq["climate_source_cells"] or len(extra)!=oq["climate_complete_cells"]:raise ValueError("Original reference support changed")
    data=rebuild(records,tests,encoded)
    a.aggregates.parent.mkdir(parents=True,exist_ok=True)
    a.aggregates.write_text(json.dumps({"version":"v405-original-annual-niche-coarse","source_QA":newqa,
       "cells":[{"key":k,"annual":v} for k,v in sorted(data.items())]},allow_nan=False))
    out=[]
    for gate in (.5,.8):
        for h in (False,True):
            ks=sorted(k for k,v in data.items() if (gate==.5 or k in land80) and (not h or high(v["2025"][0])))
            if len(ks)!={(.5,False):1323,(.5,True):742,(.8,False):1116,(.8,True):618}[gate,h]:raise ValueError("Original 2025 support changed")
            def ok(v,y):return v[y] is not None and (not h or high(v[y][0]))
            p={"land":gate,"support":"high_effort" if h else "all","original_2025_cells":len(ks),
               "annual_supported":{y:sum(ok(data[k],y) for k in ks) for y in YEARS},"fits":[]}
            common=[k for k in ks if all(ok(data[k],y) for y in YEARS)]
            p["common_three_year_cells"]=len(common)
            for j,m in enumerate(MODES):
                rr={y:[data[k][y][j] for k in common] for y in YEARS}
                f={"reference":m,"models":{}}
                for label,fields in (("competing",INTERACTION_FIELDS),("coverage_sensitivity",COVERAGE_FIELDS)):
                    f["models"][label]=paired(rr,fields,SEED+int(gate*10)+h)
                p["fits"].append(f)
            out.append(p);print(json.dumps({"v405_panel":p},allow_nan=False),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps({"version":"v405-original-annual-niche-replication","source_QA":newqa,"panels":out,
      "limits":["posthoc temporal sensitivity to v404, not independent pre-registered experiment",
       "same fixed 2018-22 realized reference and 2026 mapped railway, no establishment or opening history",
       "2025-selected common geography and annual species sets, not independent ecosystems",
       "mean rainfall affinity is not physiological drought tolerance or aridity/PET",
       "fixed reference and rarefaction draws condition spatial bootstrap; no individual data exported"]},allow_nan=False))
if __name__=="__main__":main()
