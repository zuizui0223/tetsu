#!/usr/bin/env python3
"""AVONET guild direct contrasts on identical year-common pairs.
Test habitat-open minus dense first; migration and HWI secondary.
"""
import argparse,json
from pathlib import Path
from analyze_world_railstatus_observer_annual_connectivity import (
    cells_from_pairs,subset,YEARS,W,PRIOR_BIRD,read,adapt_pair,split_pair,
    build_pairs,fit,PRED,ENDPOINT,ACTIVE,CATS,verify_original_sources,verify_replay,
    EXPECTED_SHA,MD5,POP_SHA,ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,
    read_source_builtup,source_continents,worldclim_tif,climate_scale,source_climate_lookup,
    make_lookup,original_cell_rows,inject_built,inject_status,enrich_cells,join_population)
from analyze_world_railstatus_observer_split_connectivity import symmetric_effort,EFF_PAIR
from scan_world_railstatus_bird_avonet_guilds import AXES,pair_key,MIN_GUILD_HALF
def guild_pair(base,p,gs):
    if not all(g in p["guilds"] for g in gs):return None
    out=dict(base)
    for i,g in enumerate(gs):
        z=p["guilds"][g]
        for k,v in symmetric_effort(z["split_counts_a"],z["split_counts_b"]).items():out[f"g{i}_"+k]=v
        for metric,key in (("jaccard","jaccard"),("simpson","simpson_similarity")):
            values=z["split_similarity"]
            out[f"g{i}_{metric}"]=(values["01"][key]+values["10"][key])/2
    for metric in ("jaccard","simpson"):out["difference_"+metric]=out["g1_"+metric]-out["g0_"+metric]
    return out
def models(rows):
    fields=PRED+tuple(f"g{i}_"+k for i in (0,1) for k in EFF_PAIR)
    endpoint=ENDPOINT+tuple(f"g{i}_"+k for i in (0,1) for k in EFF_PAIR if k.startswith("diff_"))
    result={}
    for mode,pred in (("block_FE",fields),("endpoint_FE",endpoint)):
        result[mode]={t:fit(rows,t,mode=="endpoint_FE",pred) for t in ("g0_jaccard","g1_jaccard","difference_jaccard","g0_simpson","g1_simpson","difference_simpson")}
        for metric in ("jaccard","simpson"):
            a,b,c=(result[mode][t+"_"+metric] for t in ("g0","g1","difference"))
            if all(v["status"]=="estimated" for v in (a,b,c)):
                if abs(b["effects"][ACTIVE]["estimate"]-a["effects"][ACTIVE]["estimate"]-c["effects"][ACTIVE]["estimate"])>1e-8:raise ValueError("Paired coefficient identity failed")
    return result
def main():
    p=argparse.ArgumentParser()
    for n in ("guild","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();guild,annual,src,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.guild,a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if guild["version"]!="railstatus-avonet-bird-guild-observer-beta-v1" or not guild["exact_full_pair_replay"] or guild["width"]!=W or guild["guild_half_gate"]!=list(MIN_GUILD_HALF) or any(guild[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong AVONET pair source")
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
        for axis,gs in AXES.items():
            by={y:[z for r in rr[y] if (z:=guild_pair(r,gg[y][r["a"],r["b"]],gs)) is not None] for y in YEARS}
            common_g=set.intersection(*({(r["a"],r["b"]) for r in v} for v in by.values()))
            p["axes"][axis]={"groups":gs,"n_native":{y:len(v) for y,v in by.items()},"n_all_year_common":len(common_g),
                             "annual_common":{y:models(subset(by[y],common_g)) for y in YEARS},
                             "annual_native":{y:models(v) for y,v in by.items()}}
        panels.append(p);print(json.dumps({"avonet_guild_panel":p},allow_nan=False),flush=True)
    obj={"version":"railstatus-avonet-bird-guild-connectivity-v1","source_replay":checks,"source":guild["source"],"panels":panels,
         "primary":"200km land50, same all-year common pair support, direct open minus dense disjoint Jaccard active crossing effect, endpoint FE; repeat all three years; migration/HWI and block FE, Simpson and land80 secondary",
         "limits":["no direct movement or railway habitat data","open vs dense is source habitat density, not forest dependence measured locally","high HWI proxy not realized dispersal","shared species and geography over years","guild-specific detection and sample size, subset selection; no causal inference"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
