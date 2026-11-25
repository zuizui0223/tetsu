#!/usr/bin/env python3
"""v396: outcome-blind local matched mapped-rail-presence contrast."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from analyze_global_naturalized_count_decomposition_v389 import source_rows,OUTCOMES,BOOT,EXPECTED
from analyze_global_naturalized_quality_effort_v390 import spatial_key,high
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_global_naturalized_population_v392 import add_population
from analyze_global_naturalized_local_swap_v393 import supported_local
from analyze_global_naturalized_joint_population_local_v394 import POP_SHA,MATCH,LOCAL_EXPECT,groups,mat,supported
from analyze_global_rail_presence_intensity_v395 import add_exposure,mixed,BASE,coefficients
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
CALIPER=.75
SEED=39620261125
MIN_PAIRS=30
MIN_BLOCKS=10
PRIOR={
 (100000,.3,.5):(211,.07608207867188073,58),
 (100000,.3,.8):(202,.05435549012076447,58),
 (100000,0.,.5):(365,.16760130936403694,123),
 (100000,0.,.8):(350,.1697233713849951,122)}
def standardized(rows):
    X=np.array([[r[k] for k in MATCH] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite matching features")
    Z=X.copy()
    ids=defaultdict(list)
    for i,r in enumerate(rows):ids[r["continent"]].append(i)
    for ix in ids.values():
        sd=X[ix].std(axis=0);sd[sd<1e-12]=1.
        Z[ix]=(X[ix]-X[ix].mean(axis=0))/sd
    return Z
def optimal_pairs(rows,Z):
    """Maximize cardinality, then minimize summed covariate distance."""
    strata=defaultdict(list)
    for i,r in enumerate(rows):strata[spatial_key(r)].append(i)
    pairs=[]
    for ix in strata.values():
        a=sorted((i for i in ix if rows[i]["rail_present"]==0),key=lambda i:rows[i]["key"])
        b=sorted((i for i in ix if rows[i]["rail_present"]==1),key=lambda i:rows[i]["key"])
        if not a or not b:continue
        dif=np.abs(Z[np.asarray(a)][:,None,:]-Z[np.asarray(b)][None,:,:])
        d=(dif*dif).sum(axis=2)
        feasible=dif.max(axis=2)<=CALIPER+1e-12
        # One additional feasible pair is always preferred to any total
        # improvement of squared distance among the remaining feasible pairs.
        penalty=1.+len(a)*len(MATCH)*CALIPER**2
        c=np.full((len(a),len(b)+len(a)),penalty)
        c[:,:len(b)]=np.where(feasible,d,penalty*(len(a)+2))
        ii,jj=linear_sum_assignment(c)
        for i,j in zip(ii,jj):
            if j<len(b):
                if not feasible[i,j]:raise ValueError("Infeasible match")
                pairs.append((a[i],b[j]))
    ids=[i for p in pairs for i in p]
    if len(ids)!=len(set(ids)):raise ValueError("Matched source cell reused")
    return pairs
def pair_weights(blocks):
    c=Counter(blocks)
    w=np.array([1./c[x] for x in blocks],float)
    return w/w.sum()
def summarize(v):
    return {"min":float(np.min(v)),"max":float(np.max(v)),"positive":int((v>0).sum())} if len(v) else None
def compare(rows,seed):
    if not rows:return {"status":"no mixed strata"}
    Z=standardized(rows)
    pairs=optimal_pairs(rows,Z)
    a=[i for i,r in enumerate(rows) if r["rail_present"]==0]
    b=[i for i,r in enumerate(rows) if r["rail_present"]==1]
    qa={"n_mixed_cells":len(rows),"n_mapped_zero":len(a),"n_mapped_positive":len(b),"n_pairs":len(pairs)}
    if not pairs:return dict(qa,status="no pairs within frozen caliper")
    ia,ib=map(np.asarray,zip(*pairs))
    blocks=[(rows[i]["continent"],rows[i]["block"]) for i in ia]
    for i,j in pairs:
        if spatial_key(rows[i])!=spatial_key(rows[j]):raise ValueError("Region/block mismatch")
    d=Z[ib]-Z[ia];w=pair_weights(blocks)
    pre=Z[b].mean(axis=0)-Z[a].mean(axis=0)
    post=np.average(d,axis=0,weights=w)
    qa.update(n_blocks=len(set(blocks)),n_regions=len({rows[i]["region"] for i in ia}),
              pairs_by_continent=dict(Counter(rows[i]["continent"] for i in ia)),
              zero_match_fraction=len(pairs)/len(a),max_single_pair_abs_standardized_diff=float(np.max(np.abs(d))),
              abs_standardized_mean_difference_before={k:float(abs(v)) for k,v in zip(MATCH,pre)},
              abs_blockweighted_paired_mean_difference_after={k:float(abs(v)) for k,v in zip(MATCH,post)},
              median_positive_rail_km=float(np.median([rows[i]["rail_km"] for i in ib])))
    if len(pairs)<MIN_PAIRS or len(set(blocks))<MIN_BLOCKS:
        return dict(qa,status="insufficient prespecified matched support")
    Y=np.array([[r[k] for k in OUTCOMES] for r in rows],float)
    delta=Y[ib]-Y[ia]
    beta=np.average(delta,axis=0,weights=w)
    if abs(beta[0]-beta[1]-beta[2])>1e-9:raise ValueError("Nonadditive count contrasts")
    qa["matched_blockweighted_difference"]=dict(zip(OUTCOMES,map(float,beta)))
    # Standardized mean residual imbalance is displayed; no fit is chosen
    # by sign or by the bootstrap interval.
    g=defaultdict(list)
    for i,k in enumerate(blocks):g[k].append(i)
    bycont=defaultdict(list)
    for k in sorted(g):bycont[k[0]].append(k)
    rng=np.random.default_rng(seed)
    draws=[]
    for _ in range(BOOT):
        ww=np.zeros(len(pairs))
        for bs in bycont.values():
            mult=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for k,m in zip(bs,mult):ww[g[k]]=m/len(g[k])
        draws.append(np.average(delta,axis=0,weights=ww))
    draws=np.asarray(draws)
    qa["block_percentile_95_descriptive"]={t:{"low":float(np.quantile(draws[:,j],.025)),
        "high":float(np.quantile(draws[:,j],.975)),"valid":len(draws)} for j,t in enumerate(OUTCOMES)}
    loo=[]
    for k in sorted(g):
        keep=[i for i,b in enumerate(blocks) if b!=k]
        if len(keep)>=MIN_PAIRS and len({blocks[i] for i in keep})>=MIN_BLOCKS:
            loo.append(np.average(delta[keep],axis=0,weights=pair_weights([blocks[i] for i in keep])))
    qa["leave_one_block_out"]={"n_valid":len(loo),"ranges":{t:summarize(np.array(loo)[:,j]) if loo else None for j,t in enumerate(OUTCOMES)}}
    qa["status"]="estimated matched association, not causal"
    return qa
def run(cells,roads,tree,pop):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:raise ValueError("Wrong original cells")
    if roads["source_zip_MD5"]!=MD5 or roads["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong roads")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree")
    if pop["version"]!="v391-original-JRC-GHS-POP-2020-1km-global-counts" or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong population")
    out=[]
    for p in cells["panels"]:
        w,pad=p["width_m"],p["region_guard_degrees"];rr=source_rows(p)
        rd=make_lookup(roads["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
        for land in LAND_THRESHOLDS:
            k=w,pad,land
            m,_=enrich_cells(rr,rd,tr,w,land);z,_=add_population(m,pp,rd)
            if len(m)!=EXPECTED[k][0] or len(z)!=len(m):raise ValueError("Original source support mismatch")
            panel={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":land,"fits":[]}
            for i,mode in enumerate(("all","high_effort")):
                r=eligible_region_rows(z if i==0 else [x for x in z if high(x)])
                local=[add_exposure(x) for x in supported_local(r)]
                if len(local)!=LOCAL_EXPECT[k][i][0]:raise ValueError("v395 original local support mismatch")
                s=mixed(local)
                if w==100000:
                    n,d,nh=PRIOR[k]
                    if len(s)!=(n if i==0 else nh):raise ValueError("v395 mixed support mismatch")
                    if i==0:
                        g,_=groups(s);fields=BASE+("rail_present",)
                        old=coefficients(mat(s,fields),g,block_weights(s),fields,("rail_present",))
                        if old is None or abs(old["rail_present"]["log1p_documented"]-d)>1e-9:raise ValueError("v395 original binary coefficient mismatch")
                seed=SEED+w+int(100*pad)+int(10*land)+i
                panel["fits"].append({"support":mode,"matched_presence":compare(s,seed)})
            out.append(panel);print(json.dumps({"v396_panel":panel},allow_nan=False),flush=True)
    return {"version":"v396-original-mapped-presence-optimal-local-pairing","panels":out,
            "matching_features":MATCH,"caliper":CALIPER,"minimum_pairs":MIN_PAIRS,
            "minimum_blocks":MIN_BLOCKS,"bootstrap_draws":BOOT,
            "limitations":["mapped absence is not verified absence",
            "same source cross-sectional observational data; matched cells differ in unmeasured habitat and detection",
            "matched pairs selected by measured covariates, not randomized",
            "small support is a result and is not repaired by changing the caliper",
            "unlisted species are not native; no causal or temporal inference"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","population","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    r=run(read(a.cells),read(a.road),read(a.tree),read(a.population));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,allow_nan=False))
if __name__=="__main__":main()
