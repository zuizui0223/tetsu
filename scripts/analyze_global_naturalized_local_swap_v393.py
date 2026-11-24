#!/usr/bin/env python3
"""v393: within-region-block naturalization enrichment and local exposure swaps."""
import argparse,json,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_naturalized_count_decomposition_v389 import (
    source_rows,fit,conditional,FIELDS,OUTCOMES,TARGETS,EXPECTED,BOOT)
from analyze_global_naturalized_quality_effort_v390 import PRIOR,spatial_key,high
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows,regression,ADJ
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
from validate_global_rail_local_matched_placebo_v385 import FEATURES,CALIPER,swapped
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
SEED=39320261125
DRAWS=199
def supported_local(rows):
    groups=defaultdict(list)
    for r in rows:groups[spatial_key(r)].append(r)
    good={k for k,v in groups.items() if len(v)>=3 and len({r["log_rail_km"] for r in v})>1}
    return [r for r in rows if spatial_key(r) in good]
def arrays(rows):
    groups=defaultdict(list);blocks=defaultdict(list)
    for i,r in enumerate(rows):
        groups[spatial_key(r)].append(i);blocks[(r["continent"],r["block"])].append(i)
    return np.array([[r[k] for k in FIELDS+OUTCOMES] for r in rows],float),[np.array(v,int) for v in groups.values()],blocks
def match(rows):
    cont=defaultdict(list);group=defaultdict(list)
    for i,r in enumerate(rows):
        cont[r["continent"]].append(i);group[spatial_key(r)].append(i)
    X=np.array([[r[k] for k in FEATURES] for r in rows],float);Z=X.copy()
    for ix in cont.values():
        sd=X[ix].std(axis=0);sd[sd<1e-12]=1.
        Z[ix]=(X[ix]-X[ix].mean(axis=0))/sd
    pairs=[]
    for ix in group.values():
        possible=[]
        for n,i in enumerate(ix):
            for j in ix[n+1:]:
                d=np.abs(Z[i]-Z[j])
                if d.max()<=CALIPER:possible.append((float(d@d),min(i,j),max(i,j)))
        used=set()
        for d,i,j in sorted(possible):
            if i in used or j in used:continue
            used.update((i,j));pairs.append((i,j))
    return pairs
def stats(rows,seed):
    r=supported_local(rows)
    X,g,blocks=arrays(r)
    o={"n_input":len(rows),"n_local_supported":len(r),"n_region_block_strata":len(g),
       "n_blocks":len(blocks),"n_regions":len({x["region"] for x in r}),
       "cells_by_continent":dict(Counter(x["continent"] for x in r))}
    if len(r)<40 or len(g)<10 or len(blocks)<10 or len(r)-len(g)-len(FIELDS)<10:
        o["status"]="insufficient local support";return o
    # Refit original region FE on precisely the local support to separate
    # the effects of removing observations and changing spatial adjustment.
    o["same_support_region_FE"]=fit(r,FIELDS,False)
    weights=block_weights(r)
    co=conditional(X,g,weights,FIELDS)
    if co is None:o["status"]="insufficient within-stratum predictor variation";return o
    o["status"]="local FE estimated";o["region_block_FE"]=co
    bycont=defaultdict(list)
    for b in sorted(blocks):bycont[b[0]].append(b)
    rng=np.random.default_rng(seed)
    samples={k:{t:[] for t in OUTCOMES} for k in TARGETS}
    for _ in range(BOOT):
        w=np.zeros(len(r))
        for bs in bycont.values():
            n=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for b,k in zip(bs,n):
                ix=blocks[b];w[ix]=k/len(ix)
        x=conditional(X,g,w,FIELDS)
        if x is None:continue
        for k in TARGETS:
            for t in OUTCOMES:samples[k][t].append(x[k][t])
    o["block_percentile_95_descriptive"]={k:{t:{"valid":len(v),
        "low":float(np.quantile(v,.025)) if len(v)>=.95*BOOT else None,
        "high":float(np.quantile(v,.975)) if len(v)>=.95*BOOT else None}
        for t,v in ss.items()} for k,ss in samples.items()}
    # This is a technical placebo for the signed local rail coefficient,
    # not a randomization inference p-value. All non-rail values remain.
    pairs=match(r)
    movable=[(i,j) for i,j in pairs if r[i]["rail_km"]!=r[j]["rail_km"]]
    qa={"matched_cells":2*len(pairs),"movable_cells":2*len(movable),
        "movable_fraction":2*len(movable)/len(r),
        "movable_blocks":len({(r[i]["continent"],r[i]["block"]) for i,j in movable})}
    o["matched_swap_support"]=qa
    if qa["movable_fraction"]<.2 or qa["movable_blocks"]<10:
        o["placebo_status"]="insufficient matched exposure mobility";return o
    vals={t:[] for t in OUTCOMES};changed=[]
    j=FIELDS.index("log_rail_km")
    for _ in range(DRAWS):
        fake,n=swapped(r,pairs,rng);changed.append(n)
        Z=X.copy();Z[:,j]=[x["log_rail_km"] for x in fake]
        cf=conditional(Z,g,weights,FIELDS)
        if cf is None:continue
        for t in OUTCOMES:vals[t].append(cf["log_rail_km"][t])
    o["placebo_status"]="conditional technical swaps, exchangeability not established"
    o["placebo_mean_changed_cells"]=float(np.mean(changed))
    o["technical_rail_coefficient_rank"]={t:{
        "actual":co["log_rail_km"][t],"valid":len(v),
        "fake_mean":float(np.mean(v)) if v else None,
        "fake_low_high":[float(min(v)),float(max(v))] if v else None,
        "fake_at_least_actual_plus_one_fraction":(1+sum(x>=co["log_rail_km"][t]-1e-12 for x in v))/(len(v)+1) if len(v)>=.95*DRAWS else None}
        for t,v in vals.items()}
    return o
def run(cells,roads,tree):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:
        raise ValueError("Not original v383 aggregate")
    if roads["source_zip_MD5"]!=MD5 or roads["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":
        raise ValueError("Not original GRIP4")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:
        raise ValueError("Not original EarthEnv")
    panels=[]
    for p in cells["panels"]:
        w,pad=p["width_m"],p["region_guard_degrees"];rr=source_rows(p)
        rd=make_lookup(roads["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        for g in LAND_THRESHOLDS:
            rows,_=enrich_cells(rr,rd,tr,w,g);base=eligible_region_rows(rows)
            n,m,b=EXPECTED[w,pad,g]
            old=fit(base,FIELDS,False)
            if len(rows)!=n or len(base)!=m or abs(old["coefficients"]["log_rail_km"]["log1p_documented"]-PRIOR[w,pad,g][0])>1e-9:
                raise ValueError("Exact v389 original source replay failed")
            out={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":g,"fits":[]}
            for mode in ("all","high_effort"):
                s=base if mode=="all" else eligible_region_rows([r for r in rows if high(r)])
                out["fits"].append({"support":mode,"local_test":stats(s,SEED+w+int(100*pad)+int(g*10)+(1 if mode=="high_effort" else 0))})
            panels.append(out);print(json.dumps({"v393_panel":out},allow_nan=False),flush=True)
    return {"version":"v393-original-local-region-block-naturalized-enrichment-falsification",
            "panels":panels,"bootstrap_draws":BOOT,"swap_draws":DRAWS,"seed":SEED,
            "limits":["same original cross-sectional observations, not independent replication",
             "1000km strata remove broad geography, not all environment or detection differences",
             "population not yet included; historical road/tree proxies",
             "matched swaps preserve local exposures but not genuine railway spatial geometry or conditional exchangeability",
             "technical signed coefficient rank is not a causal p-value",
             "observed GloNAF-documented richness and nonlisted counts, not alien/native abundance or temporal change"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.cells),read(a.road),read(a.tree));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
