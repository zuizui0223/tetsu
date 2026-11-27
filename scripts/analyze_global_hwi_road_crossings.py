#!/usr/bin/env python3
"""Actual road crossings as matched HWI comparator; no new biological scan."""
import argparse,json,math
from collections import Counter,defaultdict
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr
from analyze_world_railstatus_direct_boundary_beta import CROSS,INACTIVE,block_weights
from scan_global_grip4_road_crossings import candidates,combine
from pathlib import Path
from analyze_world_railstatus_bird_conditional_hwi import (
    cells_from_pairs,YEARS,W,read,adapt_pair,split_pair,build_pairs,
    fit,FIELDS,END,ACTIVE,CATS,verify_original_sources,verify_replay,
    EXPECTED_SHA,MD5,POP_SHA,ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,
    read_source_builtup,source_continents,worldclim_tif,climate_scale,
    source_climate_lookup,make_lookup,original_cell_rows,inject_built,
    inject_status,enrich_cells,join_population,AXES,HWI,pair_key,MIN_GUILD_HALF,SOURCE_SHA,
    standardized,guild_pair)
from scan_world_railstatus_bird_historical_hwi import PERIOD,PERIODS
from analyze_world_railstatus_bird_temporal_influence import prepare,deletion,summarize

def weights(p,axes):
    raw={}
    for ax in axes:
        gs=AXES[ax]
        if not all(g in p["guilds"] for g in gs):continue
        n=min(c[2] for g in gs for side in ("a","b") for c in p["guilds"][g]["split_counts_"+side])
        if n<MIN_GUILD_HALF[2]:raise ValueError("Original guild gate violated")
        raw[ax]=n
    total=sum(raw.values())
    return {k:v/total for k,v in raw.items()} if total else {}

def influence(rows,fields,target,endpoint):
    z=fit(rows,target,endpoint,fields)
    if z["status"]!="estimated":return z
    m=prepare(rows,fields,target,endpoint)
    if m is None:raise ValueError("Influence replay inconsistent")
    vals={k:[] for k in m["terms"]}
    for key in m["keys"]:
        q=deletion(m,{key})
        if q is not None:
            for k in vals:vals[k].append(q[k])
    z["deletion_diagnostic"]={k:summarize(v) for k,v in vals.items()}
    z["leave_continent_out"]={c:deletion(m,{k for k in m["keys"] if k[0]==c}) for c in sorted({k[0] for k in m["keys"]})}
    z["information"]={k:{q:v for q,v in e.items() if q not in ("scores","ci","estimate")} for k,e in m["effects"].items()}
    return z

def models(rows,conditional=False):
    targets=("g0_jaccard","g1_jaccard","difference_jaccard","difference_simpson")
    if conditional:targets+=("raw_difference_jaccard","standardized_minus_raw_jaccard")
    z={}
    for mode,fields,endpoint in (("block_FE",FIELDS,False),("endpoint_FE",END,True)):
        z[mode]={t:(influence(rows,fields,t,endpoint) if t=="difference_jaccard" else fit(rows,t,endpoint,fields)) for t in targets}
        if conditional:
            a,b,c=(z[mode][t] for t in ("difference_jaccard","raw_difference_jaccard","standardized_minus_raw_jaccard"))
            if all(m["status"]=="estimated" for m in (a,b,c)) and abs(a["effects"][ACTIVE]["estimate"]-b["effects"][ACTIVE]["estimate"]-c["effects"][ACTIVE]["estimate"])>1e-8:raise ValueError("Conditional decomposition failed")
    return z

def panel(rows,gg):
    raw=[z for r in rows if (z:=guild_pair(r,gg[r["a"],r["b"]],HWI)) is not None]
    keys={(r["a"],r["b"]) for r in raw}
    out={"n_raw":len(raw),"unadjusted":models(raw),"conditioned":{}}
    for label,axes in (("habitat_x_migration",tuple(AXES)),
                       ("resident_only",tuple(a for a in AXES if a.endswith("_m1"))),
                       ("dense_resident",("h1_m1",))):
        w={k:v for k in sorted(keys) if (v:=weights(gg[k],axes))}
        subset=[standardized(r,gg[r["a"],r["b"]],w[r["a"],r["b"]]) for r in rows if (r["a"],r["b"]) in w]
        out["conditioned"][label]={"n":len(subset),"strata_per_pair":dict(Counter(map(len,w.values()))),
            "pairs_by_stratum":dict(Counter(a for v in w.values() for a in v)),"models":models(subset,True)}
    return out
def fit_road(rows,target,endpoint,fields,road_field):
    blocks=defaultdict(list)
    for i,r in enumerate(rows):blocks[(r["continent"],r["block"])].append(i)
    qa={"n_pairs":len(rows),"n_blocks":len(blocks),"n_endpoints":len({r[k] for r in rows for k in ("a","b")}),
        "n_cross_active":sum(r[ACTIVE]>0 for r in rows),"n_cross_inactive":sum(r[INACTIVE]>0 for r in rows)}
    if len(blocks)<20 or len(rows)<100:return dict(qa,status="insufficient support")
    w=block_weights(rows);X=np.array([[r[k] for k in fields]+[r[target]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite pair")
    if endpoint:
        nodes={k:i for i,k in enumerate(sorted({r[q] for r in rows for q in ("a","b")}))}
        ii=np.repeat(np.arange(len(rows)),2);jj=np.array([nodes[r[k]] for r in rows for k in ("a","b")])
        D=csr_matrix((np.repeat(np.sqrt(w),2),(ii,jj)),shape=(len(rows),len(nodes)))
        X*=np.sqrt(w)[:,None]
        for j in range(X.shape[1]):
            sol=lsmr(D,X[:,j],atol=1e-12,btol=1e-12,maxiter=4000)
            if sol[1] not in (0,1,2,4,5):raise ValueError("Endpoint FE projection failed")
            X[:,j]-=D@sol[0]
        if np.max(np.abs(D.T@X))>1e-7:raise ValueError("Endpoint FE residual check")
    else:
        for ix in blocks.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
        X*=np.sqrt(w)[:,None]
    A=X[:,:-1];y=X[:,-1];norm=np.sqrt((A*A).sum(axis=0))
    # Relative to pre-projection raw feature norm; avoid rescaling numerical FE noise.
    orig=np.sqrt(np.sum(w[:,None]*np.array([[r[k] for k in fields] for r in rows],float)**2,axis=0))
    keep=norm>np.maximum(1e-9,orig*1e-8);names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE,road_field}<=set(names):return dict(qa,status="no residual variation in one or more crossing types")
    Z=A[:,keep]/norm[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return dict(qa,status="rank deficient")
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@y);e=y-Z@b;coef=b/norm[keep]
    scores=np.array([Z[ix].T@e[ix] for ix in blocks.values()])@inv.T/norm[keep]
    idx=names.index(ACTIVE)
    residual_blocks=int(sum(np.linalg.norm(Z[ix,idx])>1e-8 for ix in blocks.values()))
    qa["n_blocks_residual_active"]=residual_blocks
    qa["n_blocks_residual_road"]=int(sum(np.linalg.norm(Z[ix,names.index(road_field)])>1e-8 for ix in blocks.values()))
    if residual_blocks<20:return dict(qa,status="insufficient residual crossing support")
    def effect(terms):
        if not set(terms)<=set(names):return None
        if any(sum(np.linalg.norm(Z[ix,names.index(k)])>1e-8 for ix in blocks.values())<20 for k in terms):return None
        c=np.array([terms.get(k,0.) for k in names]);v=float(c@coef)
        se=math.sqrt(len(blocks)/(len(blocks)-1)*float(np.sum((scores@c)**2)))
        return {"estimate":v,"interval95_approx":[v-1.96*se,v+1.96*se]}
    return dict(qa,status="estimated",condition_number=float(s[0]/s[-1]),
                effects={c:effect({c:1.}) for c in (*CROSS,road_field)},active_minus_inactive=effect({ACTIVE:1.,INACTIVE:-1.}),active_minus_road=effect({ACTIVE:1.,road_field:-1.}),inactive_minus_road=effect({INACTIVE:1.,road_field:-1.}))


def compare_road(rows,by):
    rows=[dict(r,low_minus_high_jaccard=-r["difference_jaccard"],low_minus_high_simpson=-r["difference_simpson"],
               cross_highway=float(by[1].get((r["a"],r["b"]),0)>0),
               cross_primary_road=float(by[2].get((r["a"],r["b"]),0)>0)) for r in rows]
    out={"n_pairs":len(rows),"crossing_joint_counts":dict(Counter(str((int(r[ACTIVE]),int(r[INACTIVE]),int(r["cross_highway"]),int(r["cross_primary_road"]))) for r in rows)),"models":{}}
    for mode,fields,endpoint in (("block_FE",FIELDS,False),("endpoint_FE",END,True)):
        out["models"][mode]={}
        for target in ("low_minus_high_jaccard","low_minus_high_simpson"):
            out["models"][mode][target]={"prior_same_support":fit(rows,target,endpoint,fields)}
            for rf in ("cross_highway","cross_primary_road"):
                out["models"][mode][target][rf]=fit_road(rows,target,endpoint,fields+(rf,),rf)
    return out
def road_panel(rows,gg,by):
    raw=[z for r in rows if (z:=guild_pair(r,gg[r["a"],r["b"]],HWI)) is not None]
    keys={(r["a"],r["b"]) for r in raw}
    ws={k:v for k in sorted(keys) if (v:=weights(gg[k],tuple(AXES)))}
    cond=[standardized(r,gg[r["a"],r["b"]],ws[r["a"],r["b"]]) for r in rows if (r["a"],r["b"]) in ws]
    return {"raw":compare_road(raw,by),"habitat_migration_standardized":compare_road(cond,by)}

def main():
    p=argparse.ArgumentParser()
    for n in ("historical","conditional","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","road-crossings","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();hist,conditional,annual,src,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.historical,a.conditional,a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if hist["version"]!="railstatus-avonet-heldback-2018-22-hwi-pairs-v1" or not hist["exact_2025_full_and_joint_guild_replay"] or hist["source"]["source_SHA256"]!=SOURCE_SHA or hist["source"]["heldback_years"]!=list(range(2018,2023)) or hist["width"]!=W or hist["guild_half_gate"]!=list(MIN_GUILD_HALF) or conditional["version"]!="railstatus-avonet-conditional-hwi-pairs-v1" or conditional["source"]["source_SHA256"]!=SOURCE_SHA or any(len({hist[k],conditional[k],src[k],annual[k]})!=1 for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong frozen historical source")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Invalid geometry")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")

    source={y:{pair_key(p):p for p in annual["pairs"][y]["birds"]} for y in YEARS[:2]}
    source["2025"]={pair_key(p):p for p in src["pairs"][str(W)]["birds"] if "split_similarity" in p}
    common=set.intersection(*(set(v) for v in source.values()))
    if len(common)!=1645 or hist["n_frozen_candidates"]!=1645:raise ValueError("Prior geography changed")
    gg={y:{pair_key(p):p for p in hist["pairs"][y]} for y in PERIODS}
    if set(gg["2025"])!=common or not set(gg[PERIOD])<=common:raise ValueError("Heldback geography changed")
    prior_cond={y:{pair_key(p):p for p in conditional["pairs"][y]} for y in YEARS}
    if any(set(prior_cond[y])!=common for y in YEARS):raise ValueError("Prior conditional geography changed")
    for k,p in gg["2025"].items():
        if p["full"]!=source["2025"][k] or p["guilds"]!=prior_cond["2025"][k]["guilds"]:raise ValueError("2025 exact pair replay failed")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);ras,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(ras);clim,_=source_climate_lookup(land,ras,scale,{k for pair in common for k in pair},W)
    finally:ras.close();mem.close()
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)]);rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)]);pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(W)][c]} for c in CATS}
    by,receipts=combine(a.road_crossings,candidates(a.historical))
    panels=[]
    for landgate in (.5,.8):
        rr={}
        for y in PERIODS:
            raw=[p["full"] for k,p in sorted(gg[y].items())]
            rows,_,_=original_cell_rows(cells_from_pairs(raw),rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
            rows,_=inject_built(rows,built[str(W)]);rows,_=inject_status(rows,status,W);rows,_=enrich_cells(rows,rd,tr,W,landgate);rows,_=join_population(rows,pp,rd)
            pairs=build_pairs([adapt_pair(p) for p in raw],rows,el)
            rr[y]=[split_pair(r,gg[y][r["a"],r["b"]]["full"]) for r in pairs]
        if len(rr["2025"])!={.5:819,.8:649}[landgate]:raise ValueError("Prior 2025 geographic support changed")
        key_raw={ (r["a"],r["b"]) for r in rr["2025"] if all(all(g in prior_cond[y][r["a"],r["b"]]["guilds"] for g in HWI) for y in YEARS)}
        if len(key_raw)!={.5:537,.8:424}[landgate]:raise ValueError("Original all-year HWI support changed")
        old=[guild_pair(r,gg["2025"][r["a"],r["b"]],HWI) for r in rr["2025"] if (r["a"],r["b"]) in key_raw]
        replay=fit(old,"difference_jaccard",False,FIELDS)
        if replay["status"]!="estimated" or abs(replay["effects"][ACTIVE]["estimate"]-{.5:-.02039421041806581,.8:-.0010014545921946557}[landgate])>1e-9:raise ValueError("2025 HWI slope not reproduced")
        native=rr[PERIOD]
        # The primary sample is the early-period eligible subset of the
        # fixed 1645 candidate geography, never selected on HWI effects.
        result={"land":landgate,"n_full_historical":len(native),"n_2025_replay":len(rr["2025"]),"n_prior_hwi_geography":len(key_raw),
                "historical_native":road_panel(native,gg[PERIOD],by),
                "benchmark_2025_original_hwi_geography":road_panel([r for r in rr["2025"] if (r["a"],r["b"]) in key_raw],gg["2025"],by)}
        if landgate==.5:
            p=result["historical_native"]
            for key,n,expected in (("raw",693,.03170),("habitat_migration_standardized",396,.00431)):
                q=p[key]["models"]["block_FE"]["low_minus_high_jaccard"]["prior_same_support"]
                if p[key]["n_pairs"]!=n or q["status"]!="estimated" or abs(q["effects"][ACTIVE]["estimate"]-expected)>0.000006:raise ValueError("Historical HWI baseline not replayed")
        panels.append(result)
        print(json.dumps({"hwi_road_comparison":result},allow_nan=False),flush=True)
    obj={"version":"hwi-grip4-crossing-comparator-v1","source_replay":checks,"road_receipts":receipts,"panels":panels,
         "primary":"2018-22 pooled, 200km land50 native, habitat x migration standardized, low-minus-high HWI Jaccard, block FE, active minus highway direct coefficient contrast",
         "limits":["GRIP4 2018 compiled sources vs Overture 2026 geometry; mapped crossing, not historical service","binary within-feature crossing is conservative, segmentation and map completeness differ","rail-positive endpoints only, not worldwide road-only landscapes","no road or rail movement measured","historical HWI result already known, new comparator not independent biology","inadequate residual variation means unidentified"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
