#!/usr/bin/env python3
"""v403: finite regional-pool null, fixed taxonomy and occupancy composition."""
import argparse,json,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_crossfit_climate_niche_v399 import scan,encode_refs,W,MIN_REF_CELLS,MIN_REF_BLOCKS,regression,TARGET
from analyze_global_niche_shared_taxa_v401 import ref_subset,summary,metrics
from analyze_global_niche_common_support_v400 import MODES,prepare,high,INTERACTION_FIELDS,COVERAGE_FIELDS,COLD
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_client,taxonomy_source
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import source_climate_lookup
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
from analyze_global_turnover_climate_block_v387 import load_extra_climate
EXPECTED={(.5,False):1323,(.5,True):742,(.8,False):1116,(.8,True):618}
LEVELS=("occupancy","family_occupancy","genus_occupancy")
MIN_MOVABLE=.20
RESP=("observed_T_variance","expected_T_variance","excess_T_variance",
      "observed_width_C","expected_width_C","excess_width_C",
      "observed_warm_offset_C","expected_warm_offset_C","excess_warm_offset_C")
def occupancy_bin(n):
    if n<1:raise ValueError("Invalid observed occupancy")
    return 0 if n==1 else 1 if n<=3 else 2 if n<=7 else 3
def strata(taxa,occupancy,labels,level):
    out={}
    for t in taxa:
        # Unknown taxonomy is a unique stratum rather than one invented clade.
        tax=("all",) if level=="occupancy" else ("known",labels[t]) if t in labels else ("unknown",t)
        out[t]=(tax,occupancy_bin(occupancy[t]))
    return out
def expectation(selected,groups,values):
    """Exact expectation under uniform permutation within fixed strata.

    Values are one common reference per species in this block. Population
    variance, not sqrt(expected variance), is the additive estimand.
    """
    pools=defaultdict(list);counts=Counter(groups[t] for t in selected)
    for t,g in groups.items():pools[g].append(t)
    n=len(selected);sum_mu=sum_second=sum_width=var_sum=0.;mov=0;moving_groups=0
    for g,m in counts.items():
        a=np.asarray([values[t] for t in pools[g]],float);N=len(a)
        if m>N:raise ValueError("Impossible stratum count")
        x=a[:,0];mu=float(x.mean());v=float(x.var())
        sum_mu+=m*mu;sum_second+=m*float((x*x).mean());sum_width+=m*float(a[:,2].mean())
        if N>1 and m<N:
            mov+=m;moving_groups+=1
            var_sum+=m*(N-m)/(N-1)*v
    expvar=sum_second/n-(sum_mu**2+var_sum)/(n*n)
    if expvar<-1e-9:raise ValueError("Negative exact expected population variance")
    return max(0.,expvar),sum_width/n,sum_mu/n,mov/n,moving_groups
def reference_values(tid,r,encoded):
    by,val,blocks,cont=encoded;a=by.get(tid)
    if a is None:return None
    out=[]
    for mode in MODES:
        ix=ref_subset(a,blocks,cont,r,mode)
        if len(ix)<MIN_REF_CELLS or len(np.unique(blocks[ix],axis=0))<MIN_REF_BLOCKS:return None
        out.append(summary(val[ix]))
    return out
def verify_row(row,values):
    m=metrics(np.asarray(list(values.values())),row)
    for k,v in m.items():
        if abs(v-row[k])>1e-9:raise ValueError("Original v401 shared full metric mismatch: "+k)
def block_pool(records,test,encoded,family,genus):
    # All original shared-supported land50/all cells in a block define
    # one fixed pool. High effort and land80 use subsets of that pool.
    block_groups=defaultdict(list)
    for rec in records:
        r=rec["variants"][0]["shared_full"]
        block_groups[(r["continent"],tuple(r["block"]))].append(rec)
    out=[]
    for _,rr in sorted(block_groups.items()):
        exemplar=prepare(rr[0]["variants"][0]["shared_full"])
        cache={};cell_sets={};occ=Counter()
        for rec in rr:
            k=tuple(rec["key"]);ss=[]
            for t in sorted(test[k][1]):
                if t not in cache:cache[t]=reference_values(t,exemplar,encoded)
                if cache[t] is not None:ss.append(t)
            if len(ss)!=rec["common_known_species"]:raise ValueError("Original shared taxon count mismatch")
            cell_sets[k]=ss;occ.update(ss)
        tax=sorted(occ)
        labels={"occupancy":{},"family_occupancy":family,"genus_occupancy":genus}
        group_map={lev:strata(tax,occ,labels[lev],lev) for lev in LEVELS}
        for j,mode in enumerate(MODES):
            values={t:cache[t][j] for t in tax}
            for rec in rr:
                k=tuple(rec["key"]);r=prepare(rec["variants"][j]["shared_full"]);ss=cell_sets[k]
                verify_row(r,{t:values[t] for t in ss})
                for lev in LEVELS:
                    ev,ew,em,movable,ng=expectation(ss,group_map[lev],values)
                    ov=r["among_species_temperature_sd_C"]**2;ow=r["temp_width_C"];om=r["warm_offset_C"]
                    expected_offset=em-r["bio1"]
                    z=dict(r,reference_exclusion=mode,null_level=lev,pool_cells=len(rr),pool_taxa=len(tax),
                           movable_fraction=movable,movable_strata=ng,
                           observed_T_variance=ov,expected_T_variance=ev,excess_T_variance=ov-ev,
                           observed_width_C=ow,expected_width_C=ew,excess_width_C=ow-ew,
                           observed_warm_offset_C=om,expected_warm_offset_C=expected_offset,
                           excess_warm_offset_C=om-expected_offset)
                    out.append(z)
    return out
def run(a):
    read=lambda p:json.loads(p.read_text())
    src=read(a.source);plant=read(a.plants);orig=read(a.niche)
    if src["version"]!="v401-original-shared-taxon-rarefaction-coarse" or orig["version"]!="v399-crossfit-niche-coarse-aggregate":raise ValueError("Wrong original artifacts")
    q=src["source_QA"]
    if q["all_source_rows"]!=288721938 or q["eligible_by_year"]["2025"]!=7827936 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original replay")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong climate/map")
    records=[r for r in src["cells"] if r["variants"]]
    if len(records)!=1323 or len({tuple(r["key"]) for r in records})!=len(records):raise ValueError("Original v401 cell count mismatch")
    land80={tuple(r["key"]) for p in orig["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    for (land,effort),n in EXPECTED.items():
        if sum((land==.5 or tuple(r["key"]) in land80) and (not effort or high(r["variants"][0]["shared_full"])) for r in records)!=n:raise ValueError("Original v401 geographic support mismatch")
    client=source_client();canonical,genus,family,_,sha=taxonomy_source(client)
    if sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Original taxa SHA mismatch")
    refs,test,newqa=scan(client,canonical,plant)
    if newqa["eligible_by_year"]!=q["eligible_by_year"]:raise ValueError("Original yearly counts mismatch")
    keys=set(test)
    for ks in refs.values():keys.update(ks)
    land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);geo,_=source_climate_lookup(land,raster,scale,keys,W)
    finally:raster.close();mem.close()
    extra,_,_=load_extra_climate(a.climate,set(geo),W)
    encoded=encode_refs(refs,geo,extra);del refs
    if len(encoded[0])!=orig["source_QA"]["reference_taxa_with_at_least_20_climate_cells"]:raise ValueError("Reference support changed")
    rows=block_pool(records,test,encoded,family,genus)
    # Only cell-level aggregate outcomes, no taxon identities or raw coordinates.
    a.aggregates.parent.mkdir(parents=True,exist_ok=True)
    a.aggregates.write_text(json.dumps({"version":"v403-original-finite-pool-null-coarse","source_QA":newqa,"rows":rows},allow_nan=False))
    out=[]
    for land in (.5,.8):
        for h in (False,True):
            for lev in LEVELS:
                for mode in MODES:
                    base=[r for r in rows if r["reference_exclusion"]==mode and r["null_level"]==lev and
                          (land==.5 or r["key"] in land80) and (not h or high(r))]
                    if len(base)!=EXPECTED[land,h]:raise ValueError("Original source panel support changed")
                    r=[x for x in base if x["pool_cells"]>=3 and x["movable_fraction"]>=MIN_MOVABLE]
                    p={"land":land,"support":"high_effort" if h else "all","null":lev,"reference":mode,
                       "n_before":len(base),"n_after":len(r),
                       "mean_movable_fraction":float(np.mean([x["movable_fraction"] for x in base])),
                       "models":{}}
                    for label,fields in (("competing",INTERACTION_FIELDS),("coverage_sensitivity",COVERAGE_FIELDS)):
                        m=regression(r,"block",fields,(TARGET,COLD),RESP,True)
                        # Same design guarantees additive slopes for observed,
                        # expected and residual; verify to prevent contrast errors.
                        if m.get("coef"):
                            for d in m["coef"].values():
                                for tag in ("T_variance","width_C","warm_offset_C"):
                                    if abs(d["observed_"+tag]-d["expected_"+tag]-d["excess_"+tag])>1e-8:raise ValueError("Null decomposition not additive")
                        p["models"][label]=m
                    out.append(p);print(json.dumps({"v403_panel":p},allow_nan=False),flush=True)
    return {"version":"v403-original-regional-pool-composition-null","source_QA":newqa,"panels":out,
        "limits":["pool is 2025 observed shared-taxon support within a block, includes focal cell; this is a conditional label null not independent availability",
          "taxonomy strata preserve genus/family counts, not phylogenetic distances",
          "occupancy bins preserve approximate local commonness, not exact detection probability",
          "null expectation of variance is exact finite population expectation, not expected SD",
          "no actual assemblage reshuffling or independent replication, no causal railway habitat inference"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","plants","ne","climate","out","aggregates"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();o=run(a);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
