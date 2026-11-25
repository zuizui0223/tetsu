#!/usr/bin/env python3
"""Posthoc audit of the weaker-coverage family residual candidate.
Direct low-minus-high contrast on common geographic blocks, with paired
block scores and exact replay of previous high/broad coverage panels.
"""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_world_railstatus_growthform import (
    FIELDS,TARGETS,FULL,POP,ACTIVE,INACTIVE,POP_SHA,join_population,inject_status,
    local_support,CATS,LOG,read,make_lookup,enrich_cells,read_source_builtup,
    ORIGINAL_NE_SHA,WORLDCLIM_SHA,EXPECTED_SHA,MD5,sha_file,source_continents,
    worldclim_tif,climate_scale,source_climate_lookup,original_cell_rows,
    verify_original_sources,inject_built,block_weights,GRID,flexible,EXTRA,
    SHA,W,attach,fit)
OUT=("observed","expected","excess")
EXPECTED={
    (.5,.5):(1498,-.0013339714089861657),
    (.5,.25):(3002,-.004339441611927267),
    (.8,.5):(1269,-.0016225748335649897),
    (.8,.25):(2629,-.003119844914440245)}
SOURCE_HASH="822b46b1cc8bf0d902d02258516075f71c49b868886817df4f96337838df2cf4"

def groupkey(r):return r["continent"],r["block"]
def split_common(rows):
    groups=defaultdict(lambda:[[],[]])
    for r in rows:
        if r["coverage"]<.25:raise ValueError("Unexpected source coverage")
        groups[groupkey(r)][int(r["coverage"]>=.5)].append(r)
    keys=sorted(k for k,v in groups.items() if min(map(len,v))>=2)
    lo=[r for k in keys for r in groups[k][0]]
    hi=[r for k in keys for r in groups[k][1]]
    return keys,lo,hi

def arm(rows,fields,keys):
    if len(rows)<100 or len(keys)<20 or len(rows)-len(keys)-len(fields)<30:return None
    g=defaultdict(list)
    for i,r in enumerate(rows):g[groupkey(r)].append(i)
    if set(g)!=set(keys):raise ValueError("Nonidentical common block support")
    w=block_weights(rows)
    X=np.array([[r[k] for k in fields+OUT] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite design")
    for ix in g.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None];A=X[:,:len(fields)];Y=X[:,len(fields):]
    scale=np.sqrt(np.sum(A*A,axis=0));keep=scale>1e-10
    names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE}<=set(names):return None
    Z=A[:,keep]/scale[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return None
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@Y);e=Y-Z@b
    coef=b/scale[keep,None]
    score=np.array([(inv@(Z[g[k]].T@e[g[k]]))/scale[keep,None] for k in keys])
    out={}
    for label,terms in TARGETS.items():
        c=np.array([terms.get(k,0.) for k in names]);v=c@coef
        if abs(v[0]-v[1]-v[2])>1e-8:raise ValueError("Null decomposition failed")
        out[label]=(v,np.einsum("i,bio->bo",c,score))
    return out

def describe(rows):
    if not rows:return {"n":0}
    return {"n":len(rows),"n_blocks":len({groupkey(r) for r in rows}),
            "continents":dict(Counter(r["continent"] for r in rows)),
            "coverage_quantiles":list(map(float,np.quantile([r["coverage"] for r in rows],[0,.25,.5,.75,1]))),
            "mean_observations":float(np.mean([r["observations"] for r in rows])),
            "mean_species":float(np.mean([r["species"] for r in rows])),
            "mean_active_km":float(np.mean([r["active_normal_km"] for r in rows])),
            "mean_inactive_km":float(np.mean([r["inactive_normal_km"] for r in rows]))}
def compare(rows,fields):
    keys,lo,hi=split_common(rows)
    q={"low":describe(lo),"high":describe(hi),"n_common_blocks":len(keys)}
    a,b=arm(lo,fields,keys),arm(hi,fields,keys)
    if a is None or b is None:return dict(q,status="insufficient or rank-deficient common support")
    # Both arms give equal weight to the same blocks. Match the block scores
    # for a direct low-minus-high contrast, not separate significance tests.
    out={}
    for k in TARGETS:
        al,sl=a[k];bh,sh=b[k]
        z={}
        for label,v,u in (("low",al,sl),("high",bh,sh),("low_minus_high",al-bh,sl-sh)):
            se=np.sqrt(len(keys)/(len(keys)-1)*np.sum(u*u,axis=0))
            z[label]={t:{"estimate":float(x),"interval95_approx":[float(x-1.96*y),float(x+1.96*y)]}
                       for t,x,y in zip(OUT,v,se)}
        out[k]=z
    # Reproduce previous arm estimator on selected support.
    for label,rr in (("low",lo),("high",hi)):
        replay=fit(rr,fields,OUT)
        if replay["status"]!="estimated":raise ValueError("Independent arm replay failed")
        for k in TARGETS:
            for t in OUT:
                if abs(replay["effects"][k][t]["estimate"]-out[k][label][t]["estimate"])>1e-9:
                    raise ValueError("Arm coefficient mismatch")
    return dict(q,status="estimated",effects=out)

def audit(m,traits,land):
    # Exact previously observed support/coefficients; not a new primary.
    prev={}
    for c in (.5,.25):
        r=local_support(attach(m,traits,c,True));z=fit(r,FIELDS,OUT);n,b=EXPECTED[land,c]
        if len(r)!=n or z["status"]!="estimated" or abs(z["effects"]["active_minus_inactive"]["excess"]["estimate"]-b)>1e-9:
            raise ValueError("Previous woodiness panel not replayed")
        prev[str(c)]=z
    r=attach(m,traits,.25,True)
    lo=[x for x in r if x["coverage"]<.5];hi=[x for x in r if x["coverage"]>=.5]
    out={"previous":prev,"native_low":fit(local_support(lo),FIELDS,OUT),"native_high":fit(local_support(hi),FIELDS,OUT),
         "native_low_support":describe(lo),"native_high_support":describe(hi),
         "common_blocks":compare(r,FIELDS),
         "common_blocks_flexible_effort":compare([flexible(x) for x in r],FIELDS+EXTRA+("coverage",)),
         "common_blocks_high_effort":compare([x for x in r if x["observations"]>=200 and x["observers"]>=20],FIELDS),
         "common_blocks_both_positive":compare([x for x in r if x["active_normal_km"]>0 and x["inactive_normal_km"]>0],FIELDS)}
    return out

def main():
    p=argparse.ArgumentParser()
    for n in ("traits","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();src,plant,rail,status,tree,road,pop=map(read,(a.traits,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail)
    if src["version"]!="railstatus-growthform-2025-coarse-v1" or src["source"]["commit"]!=SHA or src["source_replay"]["rows"]!=288721938 or src["source_replay"]["plants"]!=7827936:raise ValueError("Wrong growthform source replay")
    if src["source"]["taxonomy_SHA256"]!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Taxonomy mismatch")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    if not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong railway status")
    if src["source"]["source_SHA256"]!=SOURCE_HASH:raise ValueError("Pinned growthform SHA mismatch")
    traits={tuple(r["key"]):r for r in src["cells"]}
    if len(traits)!=src["n_cells"]:raise ValueError("Duplicate trait cells")
    original=plant["source_grids"][str(W)][GRID]
    keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for r in original}
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);clim,_=source_climate_lookup(land,raster,scale,keys,W)
    finally:raster.close();mem.close()
    rows,_,_=original_cell_rows(original,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
    rows,_=inject_built(rows,built[str(W)]);rows,worst=inject_status(rows,status,W)
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    panels=[]
    for land in (.5,.8):
        m,_=enrich_cells(rows,rd,tr,W,land);m,_=join_population(m,pp,rd)
        p={"land":land,"base_n":len(m),**audit(m,traits,land)}
        panels.append(p);print(json.dumps({"woodiness_coverage_audit":p},allow_nan=False),flush=True)
    result={"version":"railstatus-woodiness-common-block-coverage-audit-v1","source_replay":checks,
         "growthform_source":src["source"],"rail_km_replay_error":worst,"panels":panels,
         "primary":"posthoc diagnostic of family residual: 100km land50, common blocks with >=2 cells in each coverage arm, >=20 family-known species; direct (low coverage 25-<50% minus high coverage >=50%) difference of active-minus-inactive slopes on excess woodiness",
         "limits":["coverage groups are not randomized or identical cells; common blocks do not remove within-block selection","family expectation is global external annotated pool, not regional available flora","conditional on known taxa; no imputation or full-flora identification","coarse observed landscape, not trackbed plants","approximate intervals and correlated sensitivity panels; not new confirmation"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
