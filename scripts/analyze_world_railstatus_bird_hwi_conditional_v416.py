#!/usr/bin/env python3
"""v416 conditional HWI group temporal audit with common geography."""
import argparse,json
from pathlib import Path
from analyze_world_railstatus_bird_avonet_connectivity import (
    cells_from_pairs,subset,YEARS,W,read,adapt_pair,split_pair,build_pairs,
    fit,PRED,ENDPOINT,ACTIVE,CATS,verify_original_sources,verify_replay,
    EXPECTED_SHA,MD5,POP_SHA,ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,
    read_source_builtup,source_continents,worldclim_tif,climate_scale,
    source_climate_lookup,make_lookup,original_cell_rows,inject_built,
    inject_status,enrich_cells,join_population,AXES,pair_key,MIN_GUILD_HALF,
    EFF_PAIR,guild_pair)
from analyze_world_railstatus_bird_temporal_influence import annual_audit
from scan_world_railstatus_bird_hwi_conditional_v416 import AXES,EXTRA_AXES
EXPECTED={
 (.5,"habitat"):(534,[.005745947215563904,-.01752348160790691,.010232779017653923]),
 (.5,"migration"):(620,[.010489419301141487,.005026548792989028,-.002715273438424339]),
 (.5,"HWI"):(537,[-.0428307676734211,-.02673629693574759,-.02039421041806581]),
 (.8,"habitat"):(410,None),
 (.8,"migration"):(482,[.01912770360659595,.006182228592814591,-.005074486608240293]),
 (.8,"HWI"):(424,[-.02690674145753878,-.010750591092105724,-.0010014545921946557]),
}
EXPECTED[.8,"HWI"]=(424,[-.026906541757003827,-.010750591092105724,-.0010014545921946557])
SOURCE_SHA="a8f0e0c889a05485f306029433cd45068836324e2bc6cf5b1ed8e34417f3d5fe"
FIELDS=PRED+tuple(f"g{i}_"+k for i in (0,1) for k in EFF_PAIR)
END=ENDPOINT+tuple(f"g{i}_"+k for i in (0,1) for k in EFF_PAIR if k.startswith("diff_"))

def audits(by):
    return {mode:{target:annual_audit(by,fields,target,endpoint) for target in
                  ("g0_jaccard","g1_jaccard","difference_jaccard","difference_simpson")}
            for mode,fields,endpoint in (("block_FE",FIELDS,False),("endpoint_FE",END,True))}
def main():
    p=argparse.ArgumentParser()
    for n in ("guild","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();guild,annual,src,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.guild,a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
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

        p={"land":landgate,"n_full_common":len(rr["2025"]),"axes":{}}
        axis_rows={};axis_keys={}
        for axis,gs in AXES.items():
            by={y:[z for r in rr[y] if (z:=guild_pair(r,gg[y][r["a"],r["b"]],gs)) is not None] for y in YEARS}
            common_g=set.intersection(*({(r["a"],r["b"]) for r in v} for v in by.values()))
            axis_keys[axis]=common_g
            by={y:subset(by[y],common_g) for y in YEARS};axis_rows[axis]=by
            n=len(common_g)
            m=audits(by)
            if (landgate,axis) in EXPECTED:
                old_n,expected=EXPECTED[landgate,axis]
                if n!=old_n:raise ValueError("Original guild support changed")
                replay=m["block_FE"]["difference_jaccard"]
                if expected is None:
                    if replay["status"]=="estimated":raise ValueError("Original unsupported panel unexpectedly estimated")
                else:
                    if replay["status"]!="estimated":raise ValueError("Original estimated panel now unsupported")
                    for y,v in zip(YEARS,expected):
                        if abs(replay["years"][y]["active"]["estimate"]-v)>1e-9:raise ValueError("Prior guild contrast not reproduced")
            p["axes"][axis]={"groups":gs,"n_all_year_common":n,"audit":m}
        # Compare original global split, absolute HWI in identical selected
        # species pool, and relative HWI in same family/ecology strata.
        p["common_support_comparisons"]={}
        for key in ("family","family_ecology"):
            axes=("HWI",key+"_absolute",key+"_relative")
            keys=set.intersection(*(axis_keys[a] for a in axes))
            p["common_support_comparisons"][key]={"n_pairs":len(keys),"axes":{
                a:audits({y:subset(axis_rows[a][y],keys) for y in YEARS}) for a in axes}}
        panels.append(p)
        def compact(z):
            if z["status"]!="estimated":return z
            return {k:v for k,v in z.items() if k!="years"}
        def brief(m):
            return {mode:{k:compact(v) for k,v in fits.items()} for mode,fits in m.items()}
        print(json.dumps({"v416_panel":{"land":landgate,"n_full_common":p["n_full_common"],
            "axes":{a:{"groups":v["groups"],"n":v["n_all_year_common"],"audit":brief(v["audit"])} for a,v in p["axes"].items()},
            "common_support":{key:{"n_pairs":v["n_pairs"],"axes":{a:brief(m) for a,m in v["axes"].items()}} for key,v in p["common_support_comparisons"].items()}}},allow_nan=False),flush=True)
    out={"version":"v416-conditional-HWI-connectivity-audit","source_replay":checks,"source":guild["source"],"conditional_audit":guild["conditional_audit"],"panels":panels,
         "primary":"relative HWI conditioned within family x habitat density x migration, active-minus-inactive high-minus-low Jaccard on common geography, 200km land50 block FE, three-year paired covariance; endpoint FE and Simpson, land80, family-only and preceding 2023-24 sensitivity",
         "limits":["posthoc hypothesis on same discovery data; not independent confirmation",
                   "balance is in external source matched species pool, not observed local community or full phylogeny",
                   "relative and absolute labels represent different trait contrasts; same geography comparison does not ensure same local detection",
                   "existing residual 20-block gate unchanged; do not claim if insufficient",
                   "no direct movement or actual train schedule data, no causal inference"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
