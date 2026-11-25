#!/usr/bin/env python3
"""v412 exact same-month and same-observer/month annual record balancing."""
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
 REFERENCE_YEARS,encode_refs,W,MIN_KNOWN,MIN_COVERAGE,TARGET)
from analyze_global_niche_shared_taxa_v401 import traits,metrics,NICHE_FIELDS
from analyze_global_niche_common_support_v400 import prepare,MODES,high,RESP,COVERAGE_FIELDS,INTERACTION_FIELDS,COLD
from analyze_global_niche_temporal_v405 import annual_row,YEARS,paired
from analyze_global_niche_annual_change_v410 import make_change,fit
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import source_climate_lookup
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
from analyze_global_turnover_climate_block_v387 import load_extra_climate
DRAWS=8
SEED=41220261125
def scan_strata(client,canonical,plant,keys):
    ref=defaultdict(set);tests={y:defaultdict(entry) for y in YEARS};test=tests['2025'];counts=Counter();total=0
    strata={y:defaultdict(lambda:defaultdict(list)) for y in YEARS}
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
                if r["observer_id"]:
                    oid=int(r["observer_id"]);a[2].add(oid)
                    if k in keys:
                        try:month=int(str(r["observed_on"])[5:7])
                        except (ValueError,TypeError):month=0
                        if 1<=month<=12:strata[y][k][(month,oid)].append(tid)
            else:ref[tid].add(k)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or raw.hash.hexdigest()!=plant["source_summary"]["source_iNat_original_2026_10_27_observations_gzip_SHA256_full"] or total!=288721938 or counts["2025"]!=7827936:
        raise ValueError("Original snapshot full replay mismatch")
    old={grid_key(r):(r["source_year_2025_photo_linked_plants_observation_count"],r["source_2025_plant_species_richness"],r["source_2025_distinct_observers_with_photo_linked_plant_records"]) for r in grid_counts(plant,W)}
    if old!={k:(a[0],len(a[1]),len(a[2])) for k,a in test.items()}:raise ValueError("2025 original exact cell replay failed")
    return ref,tests,strata,{"all_source_rows":total,"eligible_by_year":dict(counts),"source_observation_SHA256":raw.hash.hexdigest(),
                    "exact_2025_all_cell_effort_richness_replay":True}


def strata_arrays(by_year,mode):
    d={}
    for y in YEARS:
        g=defaultdict(list)
        for (month,oid),ids in sorted(by_year[y].items()):
            label=month if mode=="month" else (month,oid)
            g[label].extend((t,oid) for t in ids)
        d[y]=g
    common=sorted(set.intersection(*(set(d[y]) for y in YEARS)))
    n={k:min(len(d[y][k]) for y in YEARS) for k in common}
    return {y:{k:np.asarray(d[y][k],dtype=np.int64) for k in common} for y in YEARS},n
def balanced_draw(arr,n,seed):
    rng=np.random.default_rng(seed)
    if not n:return 0,set(),set()
    s=set();o=set();total=0
    for k,v in n.items():
        a=arr[k]
        x=a if len(a)==v else a[rng.choice(len(a),v,replace=False)]
        s.update(map(int,x[:,0]));o.update(map(int,x[:,1]));total+=v
    return total,s,o
def known_for(ids,caches,old,encoded):
    vv=[]
    for t in sorted(ids):
        a=[]
        for j,m in enumerate(MODES):
            if t not in caches[j]:caches[j][t]=traits(t,old[j],m,encoded)
            a.append(caches[j][t])
        if all(x is not None for x in a):vv.append(a)
    return vv
def average_rows(rows):
    # Averaged over 8 fixed record subsamples; effort/coverage gates
    # applied to EVERY draw before averaging any outcome.
    r=dict(rows[0])
    for f in ("observations","observers","species","known_niche_species","trait_coverage","log_observations","log_observers","log_species",*NICHE_FIELDS):
        r[f]=float(np.mean([x[f] for x in rows]))
    return prepare(r)
def balance_cell(k,old,annual,strata,caches,encoded,mode):
    arrays,n=strata_arrays({y:strata[y].get(k,{}) for y in YEARS},mode)
    diag={"strata":len(n),"matched_records_per_year":sum(n.values())}
    if sum(n.values())<50:return None,False,diag
    out={};hi=True
    for iy,y in enumerate(YEARS):
        draws=[[],[]]
        for d in range(DRAWS):
            a=balanced_draw(arrays[y],n,np.random.SeedSequence([SEED,int(k[0])+100000,int(k[1])+100000,iy,d,int(mode=="observer_month")]))
            vv=known_for(a[1],caches,old,encoded)
            r=[annual_row(annual[y][j],a,vv,m) for j,m in enumerate(MODES)]
            if any(x is None for x in r):return None,False,diag
            hi=hi and high(r[0])
            for j in range(2):draws[j].append(r[j])
        out[y]=[average_rows(x) for x in draws]
    return out,hi,diag
def run(a):
    read=lambda f:json.loads(f.read_text())
    annual,plant,niche=map(read,(a.annual,a.plants,a.niche));q=annual["source_QA"]
    if annual["version"]!="v405-original-annual-niche-coarse" or niche["version"]!="v399-crossfit-niche-coarse-aggregate" or q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original annual source")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong climate/map")
    data={tuple(r["key"]):r["annual"] for r in annual["cells"]}
    keys={k for k,v in data.items() if all(v[y] is not None for y in YEARS)}
    if len(data)!=1323 or len(keys)!=857:raise ValueError("Original v405 geography mismatch")
    land80={tuple(r["key"]) for p in niche["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    client=source_client();canonical,_,_,_,sha=taxonomy_source(client)
    if sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Wrong original taxa")
    refs,tests,strata,newqa=scan_strata(client,canonical,plant,keys)
    if newqa["eligible_by_year"]!=q["eligible_by_year"]:raise ValueError("Original year totals changed")
    climatekeys=set(tests["2025"])
    for ks in refs.values():climatekeys.update(ks)
    land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);geo,_=source_climate_lookup(land,raster,scale,climatekeys,W)
    finally:raster.close();mem.close()
    extra,_,_=load_extra_climate(a.climate,set(geo),W)
    encoded=encode_refs(refs,geo,extra);del refs
    oq=niche["source_QA"]
    if len(encoded[0])!=oq["reference_taxa_with_at_least_20_climate_cells"] or len(geo)!=oq["climate_source_cells"] or len(extra)!=oq["climate_complete_cells"]:raise ValueError("Original reference support changed")
    bg=defaultdict(list)
    for k in keys:
        r=prepare(data[k]["2025"][0]);bg[(r["continent"],r["block"])].append(k)
    out={};cont_cache={}
    for (cont,b),kk in sorted(bg.items()):
        caches=[{},cont_cache.setdefault(cont,{})]
        for k in sorted(kk):
            old=[prepare(r) for r in data[k]["2025"]]
            for y in YEARS:
                a0=tests[y][k];vv=known_for(a0[1],caches,old,encoded)
                for j,m in enumerate(MODES):
                    v=annual_row(old[j],a0,vv,m);prior=data[k][y][j]
                    if v is None or any(abs(v[f]-prior[f])>1e-9 for f in ("observations","observers","species","known_niche_species","trait_coverage",*NICHE_FIELDS)):raise ValueError("Original v405 annual row mismatch")
            rec={"key":k,"balanced":{}}
            for mode in ("month","observer_month"):
                rr,h,diag=balance_cell(k,old,data[k],strata,caches,encoded,mode)
                rec["balanced"][mode]={"annual":rr,"all_draws_high":h,**diag}
            out[k]=rec
        # Memory is bounded to the fixed original evaluation geography.
        for k in kk:
            for y in YEARS:strata[y].pop(k,None)
    a.aggregates.parent.mkdir(parents=True,exist_ok=True)
    a.aggregates.write_text(json.dumps({"version":"v412-original-balanced-annual-coarse","source_QA":newqa,"cells":list(out.values())},allow_nan=False))
    panels=[]
    for gate in (.5,.8):
        for h in (False,True):
            base=sorted(k for k in keys if (gate==.5 or k in land80) and (not h or all(high(data[k][y][0]) for y in YEARS)))
            if len(base)!={(.5,False):857,(.5,True):462,(.8,False):730,(.8,True):387}[gate,h]:raise ValueError("Original v405 panel changed")
            for mode in ("month","observer_month"):
                ks=[k for k in base if out[k]["balanced"][mode]["annual"] is not None and (not h or out[k]["balanced"][mode]["all_draws_high"])]
                counts=[out[k]["balanced"][mode]["matched_records_per_year"] for k in base]
                p={"land":gate,"support":"high_effort" if h else "all","balance":mode,"original_n":len(base),"n_retained":len(ks),"matched_record_quantiles":list(map(float,np.quantile(counts,[0,.05,.5,.95,1]))),"fits":[]}
                for j,m in enumerate(MODES):
                    f={"reference":m,"models":[]}
                    for variant in ("original_same_cells","balanced"):
                        rr={y:[prepare(data[k][y][j] if variant=="original_same_cells" else out[k]["balanced"][mode]["annual"][y][j]) for k in ks] for y in YEARS}
                        for label,fields in (("competing",INTERACTION_FIELDS),("coverage_sensitivity",COVERAGE_FIELDS)):
                            if not ks:
                                z={"status":"no common support","n_cells":0}
                            else:
                                aa,fs=make_change(rr,fields,"2025_minus_earlier_mean");z=fit(aa,fs,SEED+int(gate*10)+h)
                            f["models"].append({"variant":variant,"model":label,"change":z})
                    p["fits"].append(f)
                panels.append(p);print(json.dumps({"v412_panel":p},allow_nan=False),flush=True)
    return {"version":"v412-original-balanced-month-observer-niche-change","source_QA":newqa,"draws":DRAWS,"panels":panels,
      "limits":["record-level balancing is not complete detection correction or matched botanical surveys",
        "same observer-month controls identity/month mix only; taxa and locations within cells can change",
        "unknown observer or invalid month records excluded from balanced samples, not from original replay",
        "fixed reference draws and observation subsamples condition block intervals; not independent replication",
        "strict all-draw gates may strongly select Europe or leave insufficient support; no threshold relaxation",
        "no observed time-varying railway exposure or causal climatic interpretation"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("annual","plants","niche","ne","climate","out","aggregates"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();o=run(a);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
