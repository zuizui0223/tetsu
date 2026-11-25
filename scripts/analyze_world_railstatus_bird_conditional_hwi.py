#!/usr/bin/env python3
"""Fixed common-stratum standardized HWI comparisons. No threshold tuning.
Primary: identical habitat-density x migration mix within both HWI groups
and across all years. Conditioning is on recorded support, not causal.
"""
import argparse,json,math
from collections import Counter
from pathlib import Path
from analyze_world_railstatus_bird_avonet_connectivity import (
    cells_from_pairs,subset,YEARS,W,read,adapt_pair,split_pair,build_pairs,
    fit,PRED,ENDPOINT,ACTIVE,CATS,verify_original_sources,verify_replay,
    EXPECTED_SHA,MD5,POP_SHA,ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,
    read_source_builtup,source_continents,worldclim_tif,climate_scale,
    source_climate_lookup,make_lookup,original_cell_rows,inject_built,
    inject_status,enrich_cells,join_population,EFF_PAIR,guild_pair)
from analyze_world_railstatus_bird_temporal_influence import annual_audit
from scan_world_railstatus_bird_conditional_hwi import AXES,HWI,pair_key,MIN_GUILD_HALF,SOURCE_SHA
FIELDS=PRED+tuple(f"g{i}_"+k for i in (0,1) for k in EFF_PAIR)
END=ENDPOINT+tuple(f"g{i}_"+k for i in (0,1) for k in EFF_PAIR if k.startswith("diff_"))
# Frozen min over endpoint/observer half/group/year known species counts:
# independent of shared/union overlap and crossing status.
def common_stratum_weights(byyear,axes):
    raw={}
    for ax in axes:
        gs=AXES[ax]
        if not all(all(g in byyear[y]["guilds"] for g in gs) for y in YEARS):continue
        n=min(c[2] for y in YEARS for g in gs for side in ("a","b") for c in byyear[y]["guilds"][g]["split_counts_"+side])
        if n<MIN_GUILD_HALF[2]:raise ValueError("Invalid frozen guild gate")
        raw[ax]=n
    s=sum(raw.values())
    return {k:v/s for k,v in raw.items()} if s else {}

def standardized(row,p,weights):
    if not weights:return None
    parts={k:guild_pair(row,p,AXES[k]) for k in weights}
    if any(v is None for v in parts.values()) or abs(sum(weights.values())-1)>1e-10:
        raise ValueError("Fixed stratum support changed")
    z=dict(row)
    for i in (0,1):
        for k in (*EFF_PAIR,"jaccard","simpson"):
            f=f"g{i}_"+k;z[f]=math.fsum(weights[a]*v[f] for a,v in parts.items())
    for metric in ("jaccard","simpson"):
        z["difference_"+metric]=z["g1_"+metric]-z["g0_"+metric]
    raw=guild_pair(row,p,HWI)
    if raw is None:raise ValueError("Original pooled HWI pair not supported")
    for metric in ("jaccard","simpson"):
        z["raw_difference_"+metric]=raw["difference_"+metric]
        z["standardized_minus_raw_"+metric]=z["difference_"+metric]-z["raw_difference_"+metric]
    return z

def audit(rows,axes,gg):
    keys=set.intersection(*({(r["a"],r["b"]) for r in rows[y]} for y in YEARS))
    weights={k:w for k in sorted(keys) if (w:=common_stratum_weights({y:gg[y][k] for y in YEARS},axes))}
    by={y:[standardized(r,gg[y][r["a"],r["b"]],weights[r["a"],r["b"]]) for r in rows[y] if (r["a"],r["b"]) in weights] for y in YEARS}
    summary={"n_all_year_common":len(weights),"number_of_strata_per_pair":dict(Counter(map(len,weights.values()))),
             "n_pairs_contributing_by_stratum":dict(Counter(a for w in weights.values() for a in w)),
             "mean_stratum_weights":{a:sum(w.get(a,0.) for w in weights.values())/len(weights) for a in axes} if weights else {}}
    # Standardized and raw differences use the same design, to diagnose
    # attenuation from mixing strata; not the earlier unrestricted design.
    for mode,fields,endpoint in (("block_FE",FIELDS,False),("endpoint_FE",END,True)):
        summary[mode]={t:annual_audit(by,fields,t,endpoint) for t in
         ("g0_jaccard","g1_jaccard","difference_jaccard","difference_simpson","raw_difference_jaccard","standardized_minus_raw_jaccard")}
        z=summary[mode]
        if all(z[t]["status"]=="estimated" for t in ("difference_jaccard","raw_difference_jaccard","standardized_minus_raw_jaccard")):
            for y in YEARS:
                v=[z[t]["years"][y]["active"]["estimate"] for t in ("difference_jaccard","raw_difference_jaccard","standardized_minus_raw_jaccard")]
                if abs(v[0]-v[1]-v[2])>1e-8:raise ValueError("Conditional contrast decomposition failed")
    return summary

def main():
    p=argparse.ArgumentParser()
    for n in ("guild","original_guild","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();guild,original,annual,src,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.guild,a.original_guild,a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if guild["version"]!="railstatus-avonet-conditional-hwi-pairs-v1" or not guild["exact_hwi_pair_replay"] or guild["source"]["source_SHA256"]!=SOURCE_SHA or original["version"]!="railstatus-avonet-bird-guild-observer-beta-v1" or original["source"]["source"]["source_SHA256"]!=SOURCE_SHA or not guild["exact_full_pair_replay"] or guild["width"]!=W or guild["guild_half_gate"]!=list(MIN_GUILD_HALF) or any(guild[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong AVONET pair source")
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

        p={"land":landgate,"n_full_common":len(rr["2025"]),
           "all_hwi_raw":{},"conditioned":{}}
        origgg={y:{pair_key(x):x for x in original["pairs"][y]} for y in YEARS}
        if any(set(origgg[y])!=common for y in YEARS):raise ValueError("Original guild support changed")
        for y in YEARS:
            for k in common:
                for g in HWI:
                    if gg[y][k]["guilds"].get(g)!=origgg[y][k]["guilds"].get(g):
                        raise ValueError("Original HWI aggregate replay failed")
        raw={y:[z for r in rr[y] if (z:=guild_pair(r,gg[y][r["a"],r["b"]],HWI)) is not None] for y in YEARS}
        keys=set.intersection(*({(r["a"],r["b"]) for r in v} for v in raw.values()))
        raw={y:subset(raw[y],keys) for y in YEARS}
        q=annual_audit(raw,FIELDS,"difference_jaccard",False)
        expected_n,expected={.5:(537,[-.0428307676734211,-.02673629693574759,-.02039421041806581]),
                              .8:(424,[-.026906541757003827,-.010750591092105724,-.0010014545921946557])}[landgate]
        if len(keys)!=expected_n or q["status"]!="estimated" or any(abs(q["years"][y]["active"]["estimate"]-b)>1e-9 for y,b in zip(YEARS,expected)):
            raise ValueError("Original HWI result not replayed")
        p["all_hwi_raw"]=q
        # Only use pairs in the original HWI comparison. Main all-strata
        # standardization plus fixed resident and dense-resident sensitivities.
        rr0={y:subset(rr[y],keys) for y in YEARS}
        for label,axes in (("habitat_x_migration",tuple(AXES)),
                           ("resident_only",tuple(a for a in AXES if a.endswith("_m1"))),
                           ("dense_resident",("h1_m1",))):
            p["conditioned"][label]=audit(rr0,axes,gg)
        panels.append(p)
        def compact(v):
            return {k:x for k,x in v.items() if k not in ("years","validation_leave_one_block_out","validation_leave_continent_out")} if v["status"]=="estimated" else v
        print(json.dumps({"conditional_hwi_panel":{"land":landgate,"n_raw":len(keys),
            "conditioned":{l:{k:( {t:compact(v) for t,v in z.items()} if k in ("block_FE","endpoint_FE") else z) for k,z in v.items()} for l,v in p["conditioned"].items()}}},allow_nan=False),flush=True)
    out={"version":"railstatus-avonet-conditional-hwi-connectivity-v1","source_replay":checks,"source":guild["source"],"panels":panels,
         "primary":"200km land50, high minus low HWI on disjoint Jaccard after standardizing both guilds to identical habitat-density x migration mixture per pair, across all three years fixed strata and weights; block FE, direct active minus inactive as companion; no change to residual block gate",
         "limits":["posthoc after all-year HWI discovery; no independent confirmation",
                   "weight by minimum known species count across years, guilds and observer halves at endpoints; selected observed support, not population standardized",
                   "within-stratum group differences may reflect taxonomy, range size, detectability or other traits",
                   "equal mix within each pair, not common mix across all geography",
                   "no direct movement, rail habitat, management or operation data",
                   "endpoint FE or stricter supports may remain unidentified; do not relax"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
