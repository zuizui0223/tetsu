#!/usr/bin/env python3
"""v401: shared-taxon and fixed reference sample-size sensitivity."""
import argparse,json,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_crossfit_climate_niche_v399 import (
    scan,encode_refs,envelope,W,MIN_REF_CELLS,MIN_REF_BLOCKS,MIN_KNOWN,MIN_COVERAGE,OUTCOMES)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_client,taxonomy_source
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import source_climate_lookup
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import (
    WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file)
from analyze_global_turnover_climate_block_v387 import load_extra_climate
from analyze_global_niche_common_support_v400 import (
    prepare,accepted,high,model,MODES,INTERACTION_FIELDS,COVERAGE_FIELDS)
N_REF=20
DRAWS=16
SEED=40120261125
EXPECTED_COMMON={(0.5,"all"):1326,(0.5,"high_effort"):742,(0.8,"all"):1119,(0.8,"high_effort"):618}
NICHE_FIELDS=OUTCOMES[:-1]
def ref_subset(a,blocks,cont,r,mode):
    if mode==MODES[0]:return a[np.max(np.abs(blocks[a]-np.asarray(r["block"])),axis=1)>1]
    if mode==MODES[1]:return a[cont[a]!=r["continent"]]
    raise ValueError("Unknown reference exclusion")
def summary(x):
    # (n_ref,2) or (draw,n_ref,2). The reference cell dimension is -2.
    q=np.quantile(x,[.05,.95],axis=-2);mu=x.mean(axis=-2)
    return np.stack((mu[...,0],q[0,...,0],q[1,...,0]-q[0,...,0],
                     mu[...,1],q[1,...,1]-q[0,...,1]),axis=-1)
def rarefy(x,seed):
    if len(x)<N_REF:raise ValueError("Reference not supported")
    rng=np.random.default_rng(seed)
    if len(x)==N_REF:return np.repeat(summary(x)[None,:],DRAWS,axis=0)
    ix=np.asarray([rng.choice(len(x),N_REF,replace=False) for _ in range(DRAWS)],int)
    return summary(x[ix])
def traits(tid,r,mode,encoded):
    by,val,blocks,cont=encoded
    a=by.get(tid)
    if a is None:return None
    ix=ref_subset(a,blocks,cont,r,mode)
    if len(ix)<MIN_REF_CELLS or len(np.unique(blocks[ix],axis=0))<MIN_REF_BLOCKS:return None
    x=val[ix];full=summary(x)
    # A reproducible independent species x reference geographic-exclusion seed.
    # No exposure, local outcomes or model estimates enter the random seed.
    if mode==MODES[0]:
        k=(int(r["block"][0])+1000,int(r["block"][1])+1000)
    else:k=(sorted(set(cont)).index(r["continent"])+1,0)
    seed=np.random.SeedSequence([SEED,int(tid),MODES.index(mode),*k])
    return full,rarefy(x,seed),len(ix),len(np.unique(blocks[ix],axis=0))
def metrics(A,r):
    # A is (species, 5) or (draw, species, 5). Last-but-one is species.
    t,p=r["bio1"],r["log1p_bio12"]
    d=np.stack((A[...,0].mean(axis=-1)-t,A[...,1].mean(axis=-1)-t,
        (t<A[...,1]).mean(axis=-1),A[...,2].mean(axis=-1),
        A[...,0].std(axis=-1),p-A[...,3].mean(axis=-1),A[...,4].mean(axis=-1)),axis=-1)
    if d.ndim>1:d=d.mean(axis=0)
    return dict(zip(NICHE_FIELDS,map(float,d)))
def assert_replay(known,old):
    n=len(known)
    if n!=old["known_niche_species"]:raise ValueError("v399 taxon support not reproduced")
    if n:
        calc=metrics(np.asarray([x[0] for x in known]),old)
        if any(abs(calc[k]-old[k])>1e-9 for k in NICHE_FIELDS):
            raise ValueError("Original v399 community metric not reproduced")
def build(rows,other,test,encoded):
    byblock=defaultdict(list)
    for r in rows:byblock[(r["continent"],r["block"])].append(r)
    all_out={};cont_cache={}
    for (c,b),rr in sorted(byblock.items()):
        local_cache={};cont_cache_c=cont_cache.setdefault(c,{})
        for r in rr:
            ids=sorted(test[r["key"]][1]);vals=[{},{}]
            for tid in ids:
                for j,(mode,cache) in enumerate(zip(MODES,(local_cache,cont_cache_c))):
                    if tid not in cache:cache[tid]=traits(tid,r,mode,encoded)
                    if cache[tid] is not None:vals[j][tid]=cache[tid]
            # The original species set and full metrics are exactly re-created,
            # before selecting the intersection of the two reference-known sets.
            old=(r,other[r["key"]])
            for j in range(2):assert_replay(list(vals[j].values()),old[j])
            common=sorted(set(vals[0])&set(vals[1]));n=len(common);cov=n/len(ids)
            info={"key":r["key"],"common_known_species":n,"common_coverage":cov,
                  "original_known_species":[len(v) for v in vals],"total_species":len(ids)}
            variants=[]
            if n>=MIN_KNOWN and cov>=MIN_COVERAGE:
                for j in range(2):
                    vv=[vals[j][tid] for tid in common]
                    base=dict(old[j],known_niche_species=n,trait_coverage=cov)
                    # Recompute coverage×cold rather than retaining original
                    # coverage based on the unequal taxon sets.
                    full=prepare(dict(base,**metrics(np.asarray([v[0] for v in vv]),base)))
                    rare=prepare(dict(base,**metrics(np.stack([v[1] for v in vv],axis=1),base)))
                    diag={"mean_log_original_reference_cells":float(np.mean([math.log1p(v[2]) for v in vv])),
                          "mean_log_original_reference_blocks":float(np.mean([math.log1p(v[3]) for v in vv]))}
                    variants.append({"original":old[j],"shared_full":dict(full,**diag),"shared_rarefied":dict(rare,**diag)})
            info["variants"]=variants;all_out[r["key"]]=info
    return all_out
def compact(m):
    return {"status":m["status"],"n_cells":m["n_cells"],"n_blocks":m["n_blocks"],
            "coef":m.get("coefficients"),"ci":m.get("block_percentile_95_descriptive")}
def run(a):
    read=lambda f:json.loads(f.read_text())
    src=read(a.source);plant=read(a.plants)
    qa=src["source_QA"]
    if src["version"]!="v399-crossfit-niche-coarse-aggregate" or qa["all_source_rows"]!=288721938 or qa["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not qa["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original v399 source")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong original climate/map")
    old={(p["reference_exclusion"],p["min_land_fraction"]):{tuple(r["key"]):prepare(r) for r in p["rows"]} for p in src["panels"]}
    if set(old)!={(m,g) for m in MODES for g in (.5,.8)}:raise ValueError("Missing original panels")
    common={}
    for gate in (.5,.8):
        common[gate]={k for k in set(old[MODES[0],gate])&set(old[MODES[1],gate])
            if accepted(old[MODES[0],gate][k]) and accepted(old[MODES[1],gate][k])}
        if len(common[gate])!=EXPECTED_COMMON[gate,"all"] or sum(high(old[MODES[0],gate][k]) for k in common[gate])!=EXPECTED_COMMON[gate,"high_effort"]:raise ValueError("Exact v400 common geography mismatch")
    if not common[.8]<=common[.5]:raise ValueError("Unexpected land-support nesting")
    for k in common[.8]:
        for mode in MODES:
            if old[mode,.8][k]!=old[mode,.5][k]:raise ValueError("Original nested panel values differ")
    client=source_client();canonical,_,_,_,sha=taxonomy_source(client)
    if sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Wrong original taxa")
    refs,test,newqa=scan(client,canonical,plant)
    if newqa["eligible_by_year"]!=qa["eligible_by_year"]:raise ValueError("Original v399 yearly counts changed")
    keys=set(test)
    for ks in refs.values():keys.update(ks)
    land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);geo,_=source_climate_lookup(land,raster,scale,keys,W)
    finally:raster.close();mem.close()
    extra,_,_=load_extra_climate(a.climate,set(geo),W)
    encoded=encode_refs(refs,geo,extra);del refs
    if len(encoded[0])!=qa["reference_taxa_with_at_least_20_climate_cells"] or len(geo)!=qa["climate_source_cells"] or len(extra)!=qa["climate_complete_cells"]:raise ValueError("Original v399 reference support changed")
    base=[r for k,r in old[MODES[0],.5].items() if k in common[.5]]
    rebuilt=build(base,old[MODES[1],.5],test,encoded)
    # Aggregate artifact contains coarse cell summaries only; taxon IDs,
    # observer IDs, individual positions and reference cell sets not exported.
    a.aggregates.parent.mkdir(parents=True,exist_ok=True)
    a.aggregates.write_text(json.dumps({"version":"v401-original-shared-taxon-rarefaction-coarse","source_QA":newqa,
        "cells":list(rebuilt.values())},allow_nan=False))
    out=[]
    for gate in (.5,.8):
        for effort in ("all","high_effort"):
            keys=[k for k in sorted(common[gate]) if effort=="all" or high(old[MODES[0],gate][k])]
            kept=[k for k in keys if rebuilt[k]["variants"]]
            p={"min_land_fraction":gate,"support":effort,"original_common_cells":len(keys),
               "shared_taxon_supported_cells":len(kept),"shared_mean_coverage":float(np.mean([rebuilt[k]["common_coverage"] for k in keys])),
               "cells_by_continent":dict(Counter(old[MODES[0],gate][k]["continent"] for k in kept)),"fits":[]}
            for j,mode in enumerate(MODES):
                f={"reference_exclusion":mode,"models":{}}
                for variant in ("original","shared_full","shared_rarefied"):
                    rr=[rebuilt[k]["variants"][j][variant] for k in kept]
                    seed=SEED+int(gate*10)+(effort=="high_effort")
                    f["models"][variant]={"competing":model(rr,INTERACTION_FIELDS,seed),
                                          "coverage_sensitivity":model(rr,COVERAGE_FIELDS,seed)}
                p["fits"].append(f)
            out.append(p)
            print(json.dumps({"v401_panel":dict(p,fits=[{"reference_exclusion":f["reference_exclusion"],
                "models":{k:{m:compact(a) for m,a in v.items()} for k,v in f["models"].items()}} for f in p["fits"]])},allow_nan=False),flush=True)
    return {"version":"v401-original-shared-taxa-fixed-20-cell-reference","panels":out,"source_QA":newqa,
            "reference_cells_per_draw":N_REF,"rarefaction_draws":DRAWS,
            "limitations":["common taxa restricts to geographically widespread sufficiently observed plants; not all flora",
               "fixed-size reference samples are technical subsamples of realized occurrence climate; not independent studies",
               "full reference requires at least five blocks; rarefied 20 cells may span fewer blocks; does not equalize geographic extent",
               "bootstrap conditions on averaged fixed rarefaction draws, not all reference uncertainty",
               "coverage adjustment conditions on selection; no imputation of unobserved niches",
               "2025 observed landscape not true rail habitat or physiological tolerance; observational, not causal"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","plants","ne","climate","out","aggregates"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();o=run(a);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
