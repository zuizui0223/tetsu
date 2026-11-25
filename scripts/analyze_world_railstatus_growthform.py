#!/usr/bin/env python3
"""Coverage-gated woodiness association on actual global growthform source.
Primary: active versus inactive railway log length contrast on recorded
known-woodiness fraction, within spatial blocks with environment/effort.
"""
import argparse,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from analyze_world_railstatus_population_local import (
    FULL,POP,ACTIVE,INACTIVE,POP_SHA,join_population,inject_status,local_support,
    CATS,LOG,read,make_lookup,enrich_cells,read_source_builtup,ORIGINAL_NE_SHA,
    WORLDCLIM_SHA,EXPECTED_SHA,MD5,sha_file,source_continents,worldclim_tif,
    climate_scale,source_climate_lookup,original_cell_rows,verify_original_sources,
    inject_built,block_weights,GRID)
from analyze_world_railstatus_cross_taxa_population_local import flexible,EXTRA
from scan_world_railstatus_growthform import SHA,W
FIELDS=FULL+(POP,)
TARGETS={"active":{ACTIVE:1.},"inactive":{INACTIVE:1.},"active_minus_inactive":{ACTIVE:1.,INACTIVE:-1.}}

def attach(rows,traits,coverage,family=False):
    out=[]
    for r in rows:
        t=traits.get(r["key"])
        if t is None:raise ValueError("Trait source grid missing original cell")
        if any(t[k]!=r[k] for k in ("species","observations","observers")):raise ValueError("Source count mismatch")
        n=t["family_reference_known"] if family else t["known_woodiness"]
        if n<20 or n/r["species"]<coverage:continue
        obs=t["family_observed_woody"] if family else t["woody"]
        expected=t["family_expected_woody"] if family else 0.
        if not 0<=obs<=n or not 0<=expected<=n:raise ValueError("Impossible woodiness fraction")
        out.append(dict(r,coverage=n/r["species"],observed=obs/n,expected=expected/n,excess=(obs-expected)/n))
    return out

def fit(rows,fields=FIELDS,targets=("observed",)):
    g=defaultdict(list)
    for i,r in enumerate(rows):g[(r["continent"],r["block"])].append(i)
    qa={"n_cells":len(rows),"n_blocks":len(g),"continents":dict(Counter(r["continent"] for r in rows))}
    if len(rows)<100 or len(g)<20 or len(rows)-len(g)-len(fields)<30:return dict(qa,status="insufficient support")
    w=block_weights(rows)
    X=np.array([[r[k] for k in fields+targets] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite design")
    for ix in g.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None];A=X[:,:len(fields)];Y=X[:,len(fields):]
    scale=np.sqrt((A*A).sum(axis=0));keep=scale>1e-10
    names=[k for k,v in zip(fields,keep) if v]
    if not {ACTIVE,INACTIVE}<=set(names):return dict(qa,status="insufficient status variation")
    Z=A[:,keep]/scale[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return dict(qa,status="rank deficient")
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@Y);e=Y-Z@b
    coef=b/scale[keep,None]
    score=np.array([(inv@(Z[ix].T@e[ix]))/scale[keep,None] for ix in g.values()])
    effects={}
    for label,terms in TARGETS.items():
        c=np.array([terms.get(k,0.) for k in names]);v=c@coef
        if len(targets)==3 and abs(v[0]-v[1]-v[2])>1e-8:raise ValueError("Family decomposition failed")
        u=np.einsum("i,bio->bo",c,score);se=np.sqrt(len(g)/(len(g)-1)*(u*u).sum(axis=0))
        effects[label]={k:{"estimate":float(a),"interval95_approx":[float(a-1.96*b),float(a+1.96*b)]} for k,a,b in zip(targets,v,se)}
    return dict(qa,status="estimated",effects=effects)

def evaluate(rows,family):
    rr=local_support(rows);t=("observed","expected","excess") if family else ("observed",)
    # Same-support effect estimate and a flexible effort + coverage sensitivity.
    return {"n_before_local":len(rows),"primary":fit(rr,FIELDS,t),
            "coverage_effort":fit([flexible(r) for r in rr],FIELDS+EXTRA+("coverage",),t),
            "high_effort":fit(local_support([r for r in rr if r["observations"]>=200 and r["observers"]>=20]),FIELDS,t),
            "both_positive":fit(local_support([r for r in rr if r["active_normal_km"]>0 and r["inactive_normal_km"]>0]),FIELDS,t),
            "leave_continent_out":{c:fit(local_support([r for r in rr if r["continent"]!=c]),FIELDS,t) for c in sorted({r["continent"] for r in rr})}}

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
        for coverage in (.5,.25):
            z={"land":land,"coverage":coverage,"base_n":len(m),
               "known_woodiness":evaluate(attach(m,traits,coverage),False),
               "family_reference":evaluate(attach(m,traits,coverage,True),True)}
            panels.append(z);print(json.dumps({"growthform_panel":z},allow_nan=False),flush=True)
    result={"version":"railstatus-growthform-association-v1","source_replay":checks,
            "growthform_source":src["source"],"overall":src["overall"],"coverage":src["coverage"],
            "rail_km_replay_error":worst,"panels":panels,
            "primary":"100km, land>=50%, >=20 known and >=50% woodiness coverage, local block FE, active minus inactive slope on proportion woody among source-known species; negative predicted by disturbance filtering",
            "limits":["no actual railway plant observations or management intensity","family reference external global pool, not regional available species","conditional on known taxa and sampling effort; no missing-value imputation","mapped active candidate, not actual train frequency","posthoc mechanism candidate, approximate block intervals"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
