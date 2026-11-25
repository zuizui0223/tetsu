#!/usr/bin/env python3
"""v418 pair the original absolute and within-stratum HWI model estimates.
Same geographic pairs, same year/block clusters, same selected species pool.
Different trait labels define different outcomes and effort covariates.
"""
import argparse,json,math
import numpy as np
from pathlib import Path
from analyze_world_railstatus_bird_hwi_conditional_v416 import (
    cells_from_pairs,subset,YEARS,W,read,adapt_pair,split_pair,build_pairs,
    CATS,verify_original_sources,verify_replay,EXPECTED_SHA,MD5,POP_SHA,
    ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,read_source_builtup,source_continents,
    worldclim_tif,climate_scale,source_climate_lookup,make_lookup,original_cell_rows,
    inject_built,inject_status,enrich_cells,join_population,pair_key,MIN_GUILD_HALF,
    AXES,guild_pair,SOURCE_SHA,FIELDS,END)
from analyze_world_railstatus_bird_temporal_influence import prepare,average,TERMS,deletion,summarize
PERIODS={"validation_2023_2024":YEARS[:2],"descriptive_all_years":YEARS}
DIFFS={"absolute_minus_global":(-1.,1.,0.),"relative_minus_absolute":(0.,-1.,1.),
       "relative_minus_global":(-1.,0.,1.)}
def paired_average(models,years,term,weights):
    """Cluster sandwich of weighted model differences with equal-year averaging."""
    keys=models[0][years[0]]["keys"]
    if len(keys)<2 or any(m[y]["keys"]!=keys for m in models for y in years):
        raise ValueError("Different spatial clusters")
    b=[];scores=[]
    for y in years:
        b.append(sum(w*m[y]["effects"][term]["estimate"] for w,m in zip(weights,models)))
        scores.append(sum(w*m[y]["effects"][term]["scores"] for w,m in zip(weights,models)))
    S=np.asarray(scores,float);cov=len(keys)/(len(keys)-1)*S@S.T
    v=float(np.mean(b));se=math.sqrt(max(0.,float(np.sum(cov))))/len(years)
    return {"estimate":v,"interval95_approx":[v-1.96*se,v+1.96*se],"annual_covariance":cov.tolist(),"years":list(years)}
def audit(by,axes,fields,target,endpoint,prior):
    models=[]
    for ax in axes:
        m={y:prepare(by[ax][y],fields,target,endpoint) for y in YEARS}
        old=prior[ax]["endpoint_FE" if endpoint else "block_FE"][target]
        ok=all(v is not None for v in m.values())
        if ok!=(old["status"]=="estimated"):raise ValueError("Prior support changed")
        if ok:
            for label,years in PERIODS.items():
                for term in TERMS:
                    x=average(m,years,term);z=old[label][term]
                    if abs(x["estimate"]-z["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(x["interval95_approx"],z["interval95_approx"]))>1e-9:
                        raise ValueError("Exact v416 mean/interval replay failed")
        models.append(m)
    if any(v is None for m in models for v in m.values()):
        return {"status":"insufficient residual or rank support","estimated_axes":[ax for ax,m in zip(axes,models) if all(v is not None for v in m.values())]}
    keys=models[0]["2025"]["keys"]
    pairs=lambda rows:{(r["a"],r["b"]) for r in rows}
    if any(pairs(by[ax][y])!=pairs(by[axes[0]]["2025"]) for ax in axes for y in YEARS):
        raise ValueError("Unmatched original geometry")
    out={"status":"estimated","n_pairs":models[0]["2025"]["n_pairs"],"n_blocks":len(keys),"periods":{}}
    for label,years in PERIODS.items():
        z={}
        for name,w in DIFFS.items():
            z[name]={term:paired_average(models,years,term,w) for term in TERMS}
        for term in TERMS:
            v=z["relative_minus_global"][term]["estimate"]-z["absolute_minus_global"][term]["estimate"]-z["relative_minus_absolute"][term]["estimate"]
            if abs(v)>1e-10:raise ValueError("Point contrast identity failed")
        out["periods"][label]=z
    # Diagnostic leave-block / leave-continent of the primary paired difference.
    def omit(drop):
        d=[[deletion(m[y],drop) for y in YEARS] for m in models]
        if any(v is None for m in d for v in m):return None
        return {t:float(np.mean([b[t]-a[t] for a,b in zip(d[1],d[2])])) for t in TERMS}
    out["relative_minus_absolute_leave_one_block_out"]={}
    vals=[z for k in keys if (z:=omit({k})) is not None]
    for t in TERMS:out["relative_minus_absolute_leave_one_block_out"][t]=summarize([z[t] for z in vals])
    out["relative_minus_absolute_leave_continent_out"]={c:omit({k for k in keys if k[0]==c}) for c in sorted({k[0] for k in keys})}
    return out
def main():
    p=argparse.ArgumentParser()
    for n in ("guild","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","verified","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();guild,annual,src,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.guild,a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    verified=read(a.verified)
    if verified["source_run"]!=45584519015 or len(verified["panels"])!=2:raise ValueError("Wrong verified v416 source")
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if guild['source']['source']['source_SHA256']!=SOURCE_SHA:raise ValueError('Verified AVONET source SHA changed')
    if guild["version"]!="v416-conditional-HWI-pair-guilds" or not guild["exact_original_six_guild_replay"] or not guild["exact_full_pair_replay"] or guild["width"]!=W or guild["guild_half_gate"]!=list(MIN_GUILD_HALF) or any(guild[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong AVONET pair source")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Invalid geometry")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    source={y:{pair_key(p):p for p in annual["pairs"][y]["birds"]} for y in YEARS[:2]}
    source["2025"]={pair_key(p):p for p in src["pairs"][str(W)]["birds"] if "split_similarity" in p}
    common=set.intersection(*(set(v) for v in source.values()))
    if len(common)!=guild["n_all_year_common_candidates"]:raise ValueError("Candidate geography mismatch")
    gg={y:{pair_key(p):p for p in guild["pairs"][y]} for y in YEARS}
    if any(set(gg[y])!=common for y in YEARS):raise ValueError("Guild geography changed")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);ras,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(ras);clim,_=source_climate_lookup(land,ras,scale,{k for pair in common for k in pair},W)
    finally:ras.close();mem.close()
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)]);rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)]);pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(W)][c]} for c in CATS}
    panels=[]
    for landgate in (.5,.8):
        rr={}
        for y in YEARS:
            raw=[source[y][k] for k in sorted(common)]
            rows,_,_=original_cell_rows(cells_from_pairs(raw),rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
            rows,_=inject_built(rows,built[str(W)]);rows,_=inject_status(rows,status,W);rows,_=enrich_cells(rows,rd,tr,W,landgate);rows,_=join_population(rows,pp,rd)
            full=build_pairs([adapt_pair(p) for p in raw],rows,el)
            rr[y]=[split_pair(r,source[y][r["a"],r["b"]]) for r in full]
        n={len(v) for v in rr.values()}
        if n!={819 if landgate==.5 else 649}:raise ValueError("Prior annual bird common-pair support not reproduced")

        prior=next(p for p in verified["panels"] if p["land"]==landgate)
        if prior["n_full_common"]!=len(rr["2025"]):raise ValueError("Prior full support mismatch")
        axis_rows={};axis_keys={}
        for ax,gs in AXES.items():
            by={y:[z for r in rr[y] if (z:=guild_pair(r,gg[y][r["a"],r["b"]],gs)) is not None] for y in YEARS}
            keys=set.intersection(*({(r["a"],r["b"]) for r in v} for v in by.values()))
            axis_keys[ax]=keys;axis_rows[ax]=by
        panel={"land":landgate,"comparisons":{}}
        for key in ("family","family_ecology"):
            axes=("HWI",key+"_absolute",key+"_relative")
            keys=set.intersection(*(axis_keys[a] for a in axes))
            p=prior["common_support"][key]
            if len(keys)!=p["n_pairs"]:raise ValueError("Prior common support mismatch")
            by={ax:{y:subset(axis_rows[ax][y],keys) for y in YEARS} for ax in axes}
            panel["comparisons"][key]={"n_pairs":len(keys),"axes":axes,
                "models":{mode:{t:audit(by,axes,fields,t,endpoint,p["axes"]) for t in ("difference_jaccard","difference_simpson")}
                           for mode,fields,endpoint in (("block_FE",FIELDS,False),("endpoint_FE",END,True))}}
        panels.append(panel);print(json.dumps({"v418_panel":panel},allow_nan=False),flush=True)
    out={"version":"v418-paired-absolute-relative-HWI-comparison","source_run":45584519015,"source_replay":checks,"panels":panels,
        "primary":"family_ecology relative-minus-absolute active-minus-inactive high-minus-low Jaccard coefficient, 200km land50, 2023–25 paired covariance on exact v416 common geography",
        "limits":["posthoc comparison of differing trait labels and guild-specific effort designs, not decomposition of causal phylogeny",
                  "absolute and relative share selected source taxa but not observed community composition or equal absolute trait gap",
                  "both fits must meet original residual block and rank gates; no threshold relaxation",
                  "block sandwich normal approximate intervals, no study-wide multiple testing control, same discovery data"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
