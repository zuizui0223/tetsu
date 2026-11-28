"""Secondary physical network reach test on frozen HWI observations."""
import argparse,json,math
from pathlib import Path
from analyze_world_railstatus_bird_historical_hwi import (
    cells_from_pairs,YEARS,W,read,adapt_pair,split_pair,build_pairs,
    fit,FIELDS,END,ACTIVE,CATS,verify_original_sources,verify_replay,
    EXPECTED_SHA,MD5,POP_SHA,ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,
    read_source_builtup,source_continents,worldclim_tif,climate_scale,
    source_climate_lookup,make_lookup,original_cell_rows,inject_built,
    inject_status,enrich_cells,join_population,AXES,HWI,pair_key,MIN_GUILD_HALF,SOURCE_SHA,
    standardized,guild_pair,weights,PERIOD,PERIODS)
from analyze_world_railstatus_plant_insect_bird_beta import score_fit
from analyze_world_railstatus_plant_bird_disjoint import combine,interval
A="log_component_active_reach";I="log_component_inactive_reach"
INACTIVE="cross_inactive_normal"
def estimate(rows,endpoint=False,metric="jaccard"):
    fields=(END if endpoint else FIELDS)+(A,I)
    target="difference_"+metric
    f=score_fit(rows,fields,target,(A,I,ACTIVE,INACTIVE),endpoint)
    qa={"n_pairs":len(rows),"n_blocks":len({(r["continent"],r["block"]) for r in rows})}
    if f is None:return dict(qa,status="insufficient support or rank")
    # Existing guild ordering high-minus-low, negate to display low-minus-high.
    effects={k:interval(combine(f,{k:-1.})) for k in (A,I,ACTIVE,INACTIVE)}
    effects["low_minus_high_active_minus_inactive_reach"]=interval(combine(f,{A:-1.,I:1.}))
    return dict(qa,status="estimated",effects=effects)

def panel_reach(rr,lookup):
    out=[]
    for r in rr:
        k=(r["a"],r["b"])
        if k not in lookup:continue
        z=lookup[k]
        if any(bool(z[c]["crossings"])!=bool(r["cross_"+c]) for c in ("active_normal","inactive_normal")):
            raise ValueError("network crossing flag differs")
        out.append(dict(r,**{A:math.log1p(z["active_normal"]["reach_km"]),
                           I:math.log1p(z["inactive_normal"]["reach_km"])}))
    return {"n_original":len(rr),"n_topology_supported":len(out),
            "same_support_original":fit(out,"difference_jaccard",False,FIELDS),
            "primary_reach":estimate(out),"endpoint_FE":estimate(out,True),
            "simpson":estimate(out,metric="simpson")}

def main():
    p=argparse.ArgumentParser()
    for n in ("reach","historical","conditional","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
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
    net=read(a.reach)
    if net["version"]!="railstatus-component-reach-v1" or not net["exact_crossing_and_grid_km_replay"] or net["n_candidate_edges"]!=1645:
        raise ValueError("wrong network reach source")
    lookup={pair_key(p):p for p in net["edges"]}
    if len(lookup)!=net["n_supported_edges"] or not set(lookup)<=common:
        raise ValueError("wrong topology support")
    panels=[]
    for landgate in (.5,.8):
        raw=[p["full"] for k,p in sorted(gg[PERIOD].items())]
        rows,_,_=original_cell_rows(cells_from_pairs(raw),rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
        rows,_=inject_built(rows,built[str(W)]);rows,_=inject_status(rows,status,W)
        rows,_=enrich_cells(rows,rd,tr,W,landgate);rows,_=join_population(rows,pp,rd)
        pairs=build_pairs([adapt_pair(p) for p in raw],rows,el)
        rr=[split_pair(r,gg[PERIOD][r["a"],r["b"]]["full"]) for r in pairs]
        full=[z for r in rr if (z:=guild_pair(r,gg[PERIOD][r["a"],r["b"]],HWI)) is not None]
        if len(full)!={.5:693,.8:551}[landgate]:
            raise ValueError("historical raw HWI support changed")
        old=fit(full,"difference_jaccard",False,FIELDS)
        if landgate==.5 and (old["status"]!="estimated" or abs(old["effects"][ACTIVE]["estimate"]+.03170)>.0000051):
            raise ValueError("prior historical raw result not replayed")
        conditioned=[]
        for r in rr:
            k=(r["a"],r["b"]);p0=gg[PERIOD][k]
            if not all(g in p0["guilds"] for g in HWI):continue
            w=weights(p0,tuple(AXES))
            if w:conditioned.append(standardized(r,p0,w))
        if len(conditioned)!={.5:396,.8:317}[landgate]:
            raise ValueError("historical conditioned support changed")
        prior_model=fit(conditioned,"difference_jaccard",False,FIELDS)
        if landgate==.5 and (prior_model["status"]!="estimated" or abs(prior_model["effects"][ACTIVE]["estimate"]+.00431)>.0000051):
            raise ValueError("prior historical standardized result not replayed")
        z={"land":landgate,"native_raw_n":len(full),"native_conditioned_n":len(conditioned),
            "raw":panel_reach(full,lookup),"conditioned":panel_reach(conditioned,lookup)}
        panels.append(z);print(json.dumps({"historical_hwi_reach_panel":z},allow_nan=False),flush=True)
    obj={"version":"railstatus-hwi-component-reach-association-v1","source_replay":checks,
         "geometry_support":{k:v for k,v in net.items() if k!="edges"},"panels":panels,
         "primary":"2018-22 pooled, land50, original habitat x migration standardized low-minus-high HWI similarity. log1p same-component reach active-minus-inactive conditional on original crossing flags and environment",
         "limits":["secondary distinct exposure of prior observed communities, not independent biological replication",
                   "association with mapped components, not actual movement or path",
                   "map status 2026 vs past bird records, no traffic frequency",
                   "missing connector support selective; same-support original association reported",
                   "no new HWI threshold, date window, species scan or eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
