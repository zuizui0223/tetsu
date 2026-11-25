#!/usr/bin/env python3
"""Retrospective annual replication of plant-bird community similarity coupling.
Reuses original annual pair aggregates, never rescans occurrences.
"""
import argparse,json,hashlib
from pathlib import Path
from collections import Counter
import analyze_world_railstatus_observer_split_connectivity as s
from analyze_world_railstatus_plant_insect_bird_beta import (
    BIO,ORIENT,EXPECTED,attach_agb,join,fields,setup,score_fit,model,cv,combine,interval)
from analyze_world_railstatus_observer_annual_connectivity import cells_from_pairs
from scan_world_railstatus_observer_annual_beta import candidate_support,pair_key,W, GROUPS
YEARS=("2023","2024","2025")

def annual_fit(rows,epoch=2024,metric="jaccard",endpoint=False,flex=True):
    rr=setup(rows,GROUPS,epoch,flex);base=fields(GROUPS,endpoint,flex)
    qa={"n_pairs":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
        "n_endpoints":len({r[k] for r in rr for k in ("a","b")})}
    out={}
    for i,j in ORIENT:
        plant=f"plants_{i}_{metric}";target=f"birds_{j}_{metric}"
        q=score_fit(rr,base+(plant,),target,(plant,),endpoint)
        if q is None:return dict(qa,status="insufficient support or rank"),None
        out[i,j]=q[plant]
    pooled=combine(out,{ij:.5 for ij in ORIENT})
    return dict(qa,status="estimated",effect=interval(pooled)),pooled

def matched(by):
    common=set.intersection(*({(r["a"],r["b"]) for r in by[y]} for y in YEARS))
    out={}
    for y in YEARS:
        out[y]=[r for r in by[y] if (r["a"],r["b"]) in common]
        if len(out[y])!=len(common):raise ValueError("Duplicate or missing annual pair")
    # Same geography must keep same block and continent in every year.
    d=[{(r["a"],r["b"]):(r["continent"],r["block"]) for r in out[y]} for y in YEARS]
    if any(x!=d[0] for x in d[1:]):raise ValueError("Annual geography changed")
    return out

def annual_panel(by,epoch=2024,metric="jaccard",endpoint=False,flex=True):
    rr=matched(by);result={"n_common_pairs":len(rr["2025"]),"years":{}}
    f={}
    for y in YEARS:
        result["years"][y],f[y]=annual_fit(rr[y],epoch,metric,endpoint,flex)
    if all(f[y] is not None for y in YEARS[:2]):
        earlier=combine(f,{y:.5 for y in YEARS[:2]})
        result["earlier_mean"]=interval(earlier)
        if f["2025"] is not None:
            result["all_year_mean"]=interval(combine(f,{y:1/3 for y in YEARS}))
            result["2025_minus_earlier"]=interval(combine({"now":f["2025"],"earlier":earlier},{"now":1.,"earlier":-1.}))
    return result

def replay(rows,old,gate):
    if len(rows)!=old["n_pairs"]:raise ValueError("Prior pair support changed")
    for lab,kw in (("primary",{}),("flexible_effort",{"flex":True}),("simpson",{"metric":"simpson"}),("AGB2020",{"epoch":2020}),("endpoint_FE",{"endpoint":True})):
        z=model(rows,GROUPS,**kw);q=old[lab]
        if z["n_pairs"]!=q["n_pairs"] or z["status"]!=q["status"]:raise ValueError("Prior model support changed")
        if z["status"]=="estimated":
            for k,v in z["effects"].items():
                u=q["effects"][k]
                if abs(v["estimate"]-u["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],u["interval95_approx"]))>1e-9:
                    raise ValueError("Prior model effect replay failed")

def main():
    p=argparse.ArgumentParser()
    for k in ("annual","beta_model","split","prior","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    src,prior,edge,agb,am,plant,rail,status,tree,road,pop,annual,beta=map(s.read,(a.split,a.prior,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population,a.annual,a.beta_model))
    checks=s.verify_original_sources(plant,rail);s.verify_replay(prior,plant)
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_gate_observations_observers_species"]!=[25,5,10] or src["split_seed"]!="tetsu-railstatus-observer-split-v1":raise ValueError("Wrong frozen split")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")) or any(src["qa"][g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Frozen source mismatch")
    if annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or annual["width_m"]!=W or annual["qa"]["rows"]!=288721938 or any(annual[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong annual source")
    if beta["version"]!="railstatus-disjoint-community-composition-coupling-v1" or beta["biomass_source_SHA256"]!=agb["source_sha256"] or any(beta["source_SHA256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","edges","agb","agb_model")):raise ValueError("Wrong verified composition model")
    old={(r["width_m"],r["land"]):r for r in beta["panels"]}
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or am["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or am["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in s.CATS):raise ValueError("Wrong geometry")
    if tree["source_file_SHA256"]!=s.EXPECTED_SHA or road["source_zip_MD5"]!=s.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong road/tree")
    if s.sha_file(a.ne)!=s.ORIGINAL_NE_SHA or s.sha_file(a.climate)!=s.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=s.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong climate/map/pop")
    cand,cells=candidate_support(src)
    if annual["candidate_pairs"]!={g:len(cand[g]) for g in GROUPS}:raise ValueError("Frozen candidates changed")
    source={y:{g:annual["pairs"][y][g] for g in GROUPS} for y in YEARS[:2]}
    source["2025"]={g:src["pairs"][str(W)][g] for g in GROUPS}
    for y in YEARS:
        for g in GROUPS:
            keys={pair_key(r) for r in source[y][g] if "split_similarity" in r}
            if not keys<=cand[g] or (y!="2025" and len(keys)!=len(source[y][g])):raise ValueError("Annual candidates changed")
    grid=prior["grids"][str(W)]
    keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for g in GROUPS for r in grid[g]}|set.union(*cells.values())
    built,_=s.read_source_builtup(a.built);land,_=s.source_continents(a.ne);ras,mem,_=s.worldclim_tif(a.climate)
    try:
        scale,_=s.climate_scale(ras);clim,_=s.source_climate_lookup(land,ras,scale,keys,W)
    finally:ras.close();mem.close()
    tr=s.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    rd=s.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    pp=s.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    ag20=s.make_lookup(agb["grids"]["2020"][str(W)]);ag24=s.make_lookup(agb["grids"]["2024"][str(W)])
    el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(W)][c]} for c in s.CATS}
    geo={};raw={}
    for y in YEARS:
        geo[y]={};raw[y]={}
        for g in GROUPS:
            raw[y][g]={pair_key(r):r for r in source[y][g]}
            if len(raw[y][g])!=len(source[y][g]):raise ValueError("Duplicate source pair")
            records=[s.convert(z) for z in grid[g]] if y=="2025" else cells_from_pairs(source[y][g])
            rr,_,_=s.original_cell_rows(records,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
            rr,_=s.inject_built(rr,built[str(W)]);rr,_=s.inject_status(rr,status,W);geo[y][g]=rr
    panels=[]
    for gate in (.5,.8):
        by={};counts={}
        for y in YEARS:
            groups={};ag={}
            for g in GROUPS:
                rr,_=s.enrich_cells(geo[y][g],rd,tr,W,gate);rr,_=s.join_population(rr,pp,rd)
                full=s.build_pairs([s.adapt_pair(p) for p in source[y][g]],rr,el)
                groups[g]=[z for r in full if (z:=s.split_pair(r,raw[y][g][r["a"],r["b"]])) is not None]
                for k,v in attach_agb(rr,ag20,ag24,rd).items():
                    if k in ag and ag[k]!=v:raise ValueError("AGB grid mismatch")
                    ag[k]=v
            counts[y]={g:len(groups[g]) for g in GROUPS}
            by[y]=join(groups,raw[y],ag,GROUPS)
        if tuple(counts["2025"][g] for g in GROUPS)!=EXPECTED[W,gate][:2]:raise ValueError("Prior source pair counts changed")
        replay(by["2025"],old[W,gate]["plants_birds"],gate)
        z={"width_m":W,"land":gate,"native_pair_counts":counts,"native_joint_counts":{y:len(by[y]) for y in YEARS},
           "primary_flexible":annual_panel(by),"linear_effort":annual_panel(by,flex=False),
           "endpoint_FE":annual_panel(by,endpoint=True,flex=False),
           "simpson":annual_panel(by,metric="simpson",flex=False),"AGB2020":annual_panel(by,epoch=2020),
           "leave_continent_out":{c:annual_panel({y:[r for r in by[y] if r["continent"]!=c] for y in YEARS}) for c in sorted({r["continent"] for r in by["2025"]})}}
        panels.append(z);print(json.dumps({"annual_beta_coupling_panel":z},allow_nan=False),flush=True)
    out={"version":"railstatus-annual-disjoint-plant-bird-beta-v1","source_replay":checks,
         "input_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("annual","beta_model","split","edges","agb","agb_model")},
         "panels":panels,"primary":"200km land80, earlier 2023/24 mean plant-bird Jaccard association on all-year common pairs, flexible effort, AGB and geography; paired year/block scores",
         "limits":["2025-selected railway-positive geography; retrospective, not independent ecosystem validation","same platform, static railway geography and biomass","same-side endpoint observer overlap adjusted; bird observer half disjoint from plant","annual insects unavailable; no mediation or causal rail claim","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
