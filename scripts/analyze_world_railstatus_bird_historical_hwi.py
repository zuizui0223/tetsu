#!/usr/bin/env python3
"""Frozen historical pooled-period test. Do not reopen 2023-25 discovery.
Reuse previous spatial estimator and conditional standardization.
"""
import argparse,json,math
from collections import Counter
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
def main():
    p=argparse.ArgumentParser()
    for n in ("historical","conditional","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
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
                "historical_native":panel(native,gg[PERIOD]),
                "historical_on_original_hwi_geography":panel([r for r in native if (r["a"],r["b"]) in key_raw],gg[PERIOD])}
        panels.append(result)
        print(json.dumps({"historical_hwi_panel":result},allow_nan=False),flush=True)
    obj={"version":"railstatus-avonet-heldback-2018-22-hwi-connectivity-v1","source_replay":checks,"source":hist["source"],"panels":panels,
         "primary":"2018-22 pooled observations on frozen 1645 candidate geography, 200km land50, high-minus-low HWI within shared habitat x migration strata; block FE active crossing Jaccard, direct status contrast companion",
         "limits":["pooled multi-year incidence, not same temporal grain as annual discovery or independent ecosystem replication","retrospective heldback observation years on same platform and geography; earlier records may have informed unrelated work","fixed 2026 railway map, not historical operation","same observers may recur over years but groups globally disjoint","trait mix controlled within supported strata only; taxonomy, range size and detectability remain","no threshold tuning or new temporal window if unsupported"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
