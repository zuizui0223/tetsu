#!/usr/bin/env python3
"""Does railway status modify observer-disjoint plant–bird richness coupling?
Use exactly the verified prior sample; no raw occurrence scan.
All four half orientations and their paired block-score covariance retained.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,ORIENTATIONS,TERMS,EXPECTED,extract,attach,linear,interval,combine,one,temporal)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import RICH_BUILT
MODS=(b.ACTIVE,b.INACTIVE,b.BIO,"log_subway","log_tram_lightrail","log_road_density")
COMP=(RICH_BUILT,"tree_consensus_fullcell",b.POP,"bio1","bio1_sq","abs_lat","abs_lat_sq")
PRIMARY_TERMS=("plant_x_active","plant_x_inactive","plant_x_active_minus_inactive","plant_x_AGB","plant_x_road")
BASE_MODS={f"int_{k}":k for k in MODS}
NAMES={"plant_x_active":"int_"+b.ACTIVE,"plant_x_inactive":"int_"+b.INACTIVE,
       "plant_x_AGB":"int_"+b.BIO,"plant_x_road":"int_log_road_density"}

def design(rows,i,epoch,flex,competing):
    rr=b.with_epoch(b.local_support(rows),epoch);mods=MODS+(COMP if competing else ())
    for r in rr:
        r["plant_centered"]=r[f"plants_{i}_species"]-math.log1p(50.)
        # Main effects at plant richness 50; slope modifiers unaffected by centering.
        for k in mods:r["int_"+k]=r["plant_centered"]*r[k]
        if competing:r["plant_centered_sq"]=r["plant_centered"]**2
    fields=GEO+(b.BIO,b.IA,b.II)+EFF+(FLEX if flex else ())+("plant_centered",)+tuple("int_"+k for k in mods)
    if competing:fields+=("plant_centered_sq",)
    if len(set(fields))!=len(fields):raise ValueError("Duplicate design columns")
    return rr,fields

def one_interaction(rows,epoch=2024,flex=False,competing=False):
    rr=b.local_support(rows)
    qa={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
        "continents":dict(Counter(r["continent"] for r in rr))}
    # Preserve at least 20 local blocks with each status varying.
    groups=defaultdict(list)
    for r in rr:groups[r["continent"],r["block"]].append(r)
    support={k:sum(len({r[k] for r in g})>1 for g in groups.values()) for k in (b.ACTIVE,b.INACTIVE)}
    qa["status_varying_blocks"]=support
    if min(support.values(),default=0)<20:return dict(qa,status="insufficient status variation"),None
    fits={}
    for i,j in ORIENTATIONS:
        designrows,fields=design(rr,i,epoch,flex,competing)
        f={k:linear(designrows,fields,f"birds_{j}_species",v) for k,v in NAMES.items()}
        if any(v is None for v in f.values()):return dict(qa,status="insufficient support or rank"),None
        f["plant_x_active_minus_inactive"]=combine(f,{"plant_x_active":1.,"plant_x_inactive":-1.})
        fits[i,j]=f
    paired={t:{k:combine({ij:v[t] for ij,v in fits.items()},w) for k,w in TERMS.items()} for t in PRIMARY_TERMS}
    out=dict(qa,status="estimated",contrasts={t:{k:interval(v) for k,v in q.items()} for t,q in paired.items()})
    return out,paired

def evaluate_modifiers(rows):
    rr=b.local_support(rows)
    return {"n_cells":len(rr),"primary":one_interaction(rr)[0],
            "flexible_effort":one_interaction(rr,flex=True)[0],
            "competing_environment_and_nonlinear_plant":one_interaction(rr,flex=True,competing=True)[0],
            "AGB2020_same_cells":one_interaction(rr,2020)[0],
            "high_effort":one_interaction([r for r in rr if min(r["observations"],r["bird_observations"])>=200 and min(r["observers"],r["bird_observers"])>=20])[0],
            "both_status_positive":one_interaction([r for r in rr if min(r["active_normal_km"],r["inactive_normal_km"])>0])[0],
            "leave_continent_out":{c:one_interaction([r for r in rr if r["continent"]!=c])[0] for c in sorted({r["continent"] for r in rr})}}

def temporal_modifiers(by):
    keys=set.intersection(*({r["key"] for r in by[y]} for y in YEARS))
    by={y:b.local_support([r for r in by[y] if r["key"] in keys]) for y in YEARS}
    if any({r["key"] for r in by[y]}!={r["key"] for r in by["2025"]} for y in YEARS):raise ValueError("Annual geography changed")
    out={"n_common_cells":len(by["2025"])}
    for label,flex in (("linear",False),("flexible",True)):
        results={};scores={}
        for y in YEARS:results[y],scores[y]=one_interaction(by[y],flex=flex)
        z={"years":results}
        if all(scores.values()):
            for k,ys in (("earlier_mean",YEARS[:2]),("all_year_mean",YEARS)):
                z[k]={t:{c:interval(combine({y:scores[y][t][c] for y in ys},{y:1/len(ys) for y in ys})) for c in TERMS} for t in PRIMARY_TERMS}
        out[label]=z
    return out

def replay_check(new,old):
    if new["status"]!=old["status"] or new["n_cells"]!=old["n_cells"]:raise ValueError("Previous split sample changed")
    if new["status"]=="estimated":
        for k,e in new["contrasts"].items():
            p=old["contrasts"][k]
            if abs(e["estimate"]-p["estimate"])>1e-9 or max(abs(x-y) for x,y in zip(e["interval95_approx"],p["interval95_approx"]))>1e-9:raise ValueError("Previous split slope/interval changed")

def main():
    p=argparse.ArgumentParser()
    for n in ("agb","agb_model","disjoint","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();agb,old,prior,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.agb,a.agb_model,a.disjoint,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if prior["version"]!="railstatus-observer-disjoint-plant-bird-biomass-v1" or prior["biomass_source_sha256"]!=agb["source_sha256"] or any(prior["source_artifact_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("agb","agb_model","split","annual")):raise ValueError("Prior disjoint artifact mismatch")
    prior_panels={(p["width_m"],p["land"]):p for p in prior["panels"]}
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or old["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or old["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass replay mismatch")
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_seed"]!=SALT.decode() or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or src["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Frozen observer split invalid")
    if annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or annual["width_m"]!=200000 or annual["qa"]["rows"]!=288721938 or any(len({src[k],annual[k],taxa[k]})!=1 for k in ("obs_sha256","taxa_sha256")) or annual["split_seed"]!=src["split_seed"]:raise ValueError("Annual source mismatch")
    cand,_=candidate_support(src)
    if annual["candidate_pairs"]!={g:len(cand[g]) for g in GROUPS}:raise ValueError("Frozen annual candidates changed")
    for y in YEARS[:2]:
        for g in GROUPS:
            if not {pair_key(p) for p in annual["pairs"][y][g]}<=cand[g]:raise ValueError("Annual candidates changed")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen map/climate/pop mismatch")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);raster,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(raster)
        for w in b.GRID_WIDTHS:
            grid=taxa["grids"][str(w)];yr=YEARS if w==200000 else ("2025",)
            ps={y:{g:src["pairs"][str(w)][g] if y=="2025" else annual["pairs"][y][g] for g in GROUPS} for y in yr}
            maps={y:{g:extract(ps[y][g]) for g in GROUPS} for y in yr}
            prior_counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in GROUPS}
            for g in GROUPS:
                if any(prior_counts[g].get(k)!=v[0] for k,v in maps["2025"][g].items()):raise ValueError("2025 exact source replay failed")
            keys={b.key(z) for g in GROUPS for z in grid[g]}|{k for y in yr for g in GROUPS for k in maps[y][g]}
            clim,_=b.source_climate_lookup(land,raster,scale,keys,w)
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            geo={}
            for y in yr:
                geo[y]={}
                for g in GROUPS:
                    raw=[b.convert(z) for z in grid[g]] if y=="2025" else cells_from_pairs([p for p in ps[y][g] if "split_similarity" in p])
                    rows,_,_=b.original_cell_rows(raw,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                    rows,_=b.inject_built(rows,built[str(w)]);rows,_=b.inject_status(rows,status,w);geo[y][g]=rows
            for gate in b.LAND_THRESHOLDS:
                by={}
                for y in yr:
                    rows={}
                    for g in GROUPS:
                        m,_=b.enrich_cells(geo[y][g],rd,tr,w,gate);m,_=b.join_population(m,pp,rd)
                        rows[g]=[b.flexible(r) for r in m]
                    joined,_=b.join(rows,ag20,ag24,w,rd)
                    bc={r["key"]:r["species"] for r in rows["birds"]}
                    joined=[dict(r,bird_species=bc[r["key"]]) for r in joined]
                    if y=="2025":
                        n,expected=EXPECTED[w][gate];m=b.models(joined,2024)
                        q=m["bird_plant_richness_association"]["effects"]["plant_richness"]["bird_log_richness"]["estimate"]
                        if m["n"]!=n or abs(q-expected)>5.1e-6:raise ValueError("Prior full-count plant/bird coefficient not replayed")
                    by[y]=attach(joined,maps[y])
                oldp=prior_panels[w,gate]
                for name,flex in (("primary",False),("flexible_effort",True)):
                    replay_check(one(by["2025"],flex=flex)[0],oldp["2025"][name])
                z={"width_m":w,"land":gate,"2025":evaluate_modifiers(by["2025"])}
                if w==200000:
                    replay=temporal(by)
                    for label in ("linear","flexible"):
                        for y in YEARS:replay_check(replay[label]["years"][y],oldp["common_annual"][label]["years"][y])
                    z["common_annual"]=temporal_modifiers(by)
                    if z["common_annual"]["n_common_cells"]!=oldp["common_annual"]["n_common_cells"]:raise ValueError("Annual support changed")
                panels.append(z);print(json.dumps({"plant_bird_coupling_modifier_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-observer-disjoint-plant-bird-coupling-modifiers-v1","source_replay":checks,"prior_disjoint_SHA256":hashlib.sha256(a.disjoint.read_bytes()).hexdigest(),"biomass_source_sha256":agb["source_sha256"],"panels":panels,
      "primary":"100km land80 2025 same 1640 cells, mean opposite-observer-side plant richness x (log active km minus log inactive km coefficient contrast) in joint covariate model. Plant x woody AGB and road modifiers also reported; 200km, earlier common cells, nonlinear effort/competing environmental moderators secondary.",
      "limits":["cross-sectional conditional richness association, not biotic interaction or causation","same minus cross diagnostic, not bias amount","same platform, geographic and access/detection confounding remain","static 2026 railway map","woody biomass not productivity or bird biomass","posthoc proposed extension; earlier years already used for main coupling but not this modifier","multiple correlated sensitivity panels"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
