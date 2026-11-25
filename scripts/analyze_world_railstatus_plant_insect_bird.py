#!/usr/bin/env python3
"""Does insect richness add information about bird richness on matched cells?
Plant and insect predictors share one observer half; bird response uses the
opposite half. No claim of three pairwise-disjoint observer groups or mediation.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,EXPECTED,extract,attach,linear,interval,combine,one)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import fit_predict,metrics,CONTINENTS
from analyze_world_railstatus_plant_bird_coupling_modifiers import replay_check
CROSS=((0,1),(1,0))
ADDED=("insects","fungi")
def extra_attach(rows,maps,full,groups):
    out=[]
    for r in rows:
        k=r["key"]
        if any(k not in maps[g] for g in groups):continue
        z=dict(r)
        for g in groups:
            f,hs=maps[g][k]
            if full[g].get(k)!=f:raise ValueError("Exact added-taxon full count mismatch")
            z[g+"_full_observations"],z[g+"_full_observers"],_=f
            for i,h in enumerate(hs):
                for name,n in zip(("observations","observers","species"),h):
                    x=f"{g}_{i}_{name}";z[x]=math.log1p(n)
                    if name!="species":z[x+"_sq"]=z[x]**2;z[x+"_cu"]=z[x]**3
                z[f"{g}_{i}_effort_interaction"]=z[f"{g}_{i}_observations"]*z[f"{g}_{i}_observers"]
        out.append(z)
    return out

def fields_for(groups,flex):
    eff=tuple(f"{g}_{i}_{k}" for g in groups for i in (0,1) for k in ("observations","observers"))
    extra=tuple(f"{k}_{p}" for k in eff for p in ("sq","cu"))+tuple(f"{g}_{i}_effort_interaction" for g in groups for i in (0,1))
    fields=GEO+(b.BIO,)+EFF+eff+(FLEX+extra if flex else ())
    if len(set(fields))!=len(fields):raise ValueError("Duplicate design")
    return fields

def cv(rows,fields,groups):
    # Across continents, no local FE intercept fitted on test outcomes.
    # Baseline and extended models have exactly the same effort covariates,
    # same training/test rows and frozen ridge. Predictor side differs from bird.
    folds=[]
    for c in CONTINENTS:
        te=[r for r in rows if r["continent"]==c];tr=[r for r in rows if r["continent"]!=c]
        if len(te)<50 or len(tr)<max(300,4*(len(fields)+len(groups)+1)):continue
        z=[]
        for i,j in CROSS:
            plant=f"plants_{i}_species";target=f"birds_{j}_species"
            full=fields+(plant,)+tuple(f"{g}_{i}_species" for g in groups)
            scores={"full":metrics(te,fit_predict(tr,te,full,target,False),target)[0]}
            for g in ("plants",*groups):
                pred=tuple(k for k in full if k!=f"{g}_{i}_species")
                scores["without_"+g]=metrics(te,fit_predict(tr,te,pred,target,False),target)[0]
            z.append(scores)
        m={k:sum(v[k] for v in z)/2 for k in z[0]}
        folds.append({"continent":c,"n_test":len(te),"n_train":len(tr),"full_MSE":m["full"],
                      "incremental_dMSE":{g:m["without_"+g]-m["full"] for g in ("plants",*groups)}})
    return {"n_folds":len(folds),"adequate_folds":len(folds)>=4,"folds":folds,
            "mean_delta_MSE":{g:sum(v["incremental_dMSE"][g] for v in folds)/len(folds) for g in ("plants",*groups)} if folds else {},
            "n_positive_folds":{g:sum(v["incremental_dMSE"][g]>0 for v in folds) for g in ("plants",*groups)}}

def model(rows,groups,epoch=2024,flex=False,do_cv=False):
    rr=b.with_epoch(b.local_support(rows),epoch);fields=fields_for(groups,flex)
    qa={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
        "continents":dict(Counter(r["continent"] for r in rr))}
    orientations={}
    for i,j in CROSS:
        plant=f"plants_{i}_species";target=f"birds_{j}_species"
        base=fields+(plant,);full=base+tuple(f"{g}_{i}_species" for g in groups)
        z={"plant_without_added_richness":linear(rr,base,target,plant),
           "plant_with_added_richness":linear(rr,full,target,plant)}
        z.update({g:linear(rr,full,target,f"{g}_{i}_species") for g in groups})
        if any(v is None for v in z.values()):return dict(qa,status="insufficient support or rank")
        z["plant_slope_attenuation"]=combine(z,{"plant_without_added_richness":1.,"plant_with_added_richness":-1.})
        if len(groups)==2:z["insects_minus_fungi"]=combine(z,{"insects":1.,"fungi":-1.})
        orientations[i,j]=z
    effects={k:interval(combine({ij:v[k] for ij,v in orientations.items()},{ij:.5 for ij in CROSS})) for k in orientations[CROSS[0]]}
    if abs(effects["plant_without_added_richness"]["estimate"]-effects["plant_with_added_richness"]["estimate"]-effects["plant_slope_attenuation"]["estimate"])>1e-9:
        raise ValueError("Paired attenuation identity failed")
    out=dict(qa,status="estimated",effects=effects)
    if do_cv:out["continent_holdout"]=cv(rr,fields,groups)
    return out

def panel(rows,groups):
    rr=b.local_support(rows)
    return {"groups":list(groups),"n_supported":len(rr),"primary":model(rr,groups,do_cv=True),
      "flexible_effort":model(rr,groups,flex=True,do_cv=True),"AGB2020_same_cells":model(rr,groups,2020),
      "high_effort":model([r for r in rr if min(r["observations"],r["bird_observations"],*(r[g+"_full_observations"] for g in groups))>=200 and min(r["observers"],r["bird_observers"],*(r[g+"_full_observers"] for g in groups))>=20],groups),
      "leave_continent_out":{c:model([r for r in rr if r["continent"]!=c],groups) for c in sorted({r["continent"] for r in rr})}}

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
            grid=taxa["grids"][str(w)];yr=("2025",)
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
                oldp=prior_panels[w,gate]["2025"]
                replay_check(one(by["2025"])[0],oldp["primary"])
                replay_check(one(by["2025"],flex=True)[0],oldp["flexible_effort"])
                maps_extra={g:extract(src["pairs"][str(w)][g]) for g in ADDED}
                full_extra={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in ADDED}
                if any(full_extra[g].get(k)!=v[0] for g in ADDED for k,v in maps_extra[g].items()):raise ValueError("Added taxa full replay changed")
                r=by["2025"]
                z={"width_m":w,"land":gate,"insect":panel(extra_attach(r,maps_extra,full_extra,("insects",)),("insects",)),
                   "fungi_comparator":panel(extra_attach(r,maps_extra,full_extra,("fungi",)),("fungi",)),
                   "both_same_cells":panel(extra_attach(r,maps_extra,full_extra,ADDED),ADDED)}
                panels.append(z);print(json.dumps({"plant_insect_bird_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-observer-disjoint-plant-insect-bird-biomass-v1","source_replay":checks,
       "prior_disjoint_SHA256":hashlib.sha256(a.disjoint.read_bytes()).hexdigest(),"biomass_source_sha256":agb["source_sha256"],"panels":panels,
       "primary":"100km land80 2025 plant and insect richness in same observer half, bird in opposite half. Matched block-FE conditional insect slope; paired change in plant slope after adding insect richness; holdout predictive increments as distinct evidence. Joint efforts all six halves; primary linear, flexible and 200km sensitivity.",
       "limits":["plant and insect observers may overlap; only bird versus predictors disjoint",
                 "all insects not prey assemblage; fungi comparator not negative control; no causal mediation",
                 "same platform, residual environment and accessibility confounding",
                 "prediction conditional on observed other-taxon richness and effort, not independent range validation",
                 "whole landscape observed richness, not abundance or bird biomass","no eBird data used"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
