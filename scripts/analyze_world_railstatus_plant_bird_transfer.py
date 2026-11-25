#!/usr/bin/env python3
"""Frozen leave-continent-out predictive increment on full observer-disjoint support.
Reuse the source reconstruction of the verified disjoint association; no raw scan.
"""
import argparse, json, hashlib
from pathlib import Path
from analyze_world_railstatus_plant_bird_disjoint import (
    b, YEARS, GROUPS, GEO, EFF, FLEX, EXPECTED, SALT, MIN_FULL, MIN_SPLIT,
    cells_from_pairs, candidate_support, pair_key, extract, attach, one)
from analyze_world_railstatus_plant_bird_coupling_modifiers import replay_check
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    fit_predict, metrics, CONTINENTS, RIDGE)
from collections import Counter

def transfer(rows, donor="plants", recipient="birds", flex=True, epoch=2024):
    if {donor, recipient} != {"plants", "birds"}:
        raise ValueError("Only the two frozen directions are allowed")
    rows = b.with_epoch(b.local_support(rows),epoch)
    base = GEO + EFF + (FLEX if flex else ())
    if b.BIO in base or len(set(base)) != len(base):
        raise ValueError("Duplicated predictors")
    folds=[]
    for c in CONTINENTS:
        train=[r for r in rows if r["continent"]!=c]
        test=[r for r in rows if r["continent"]==c]
        if len(test)<50 or len(train)<max(300,4*(len(base)+2)):
            continue
        scores=[]
        for i,j in ((0,1),(1,0)):
            x=f"{donor}_{i}_species"
            y=f"{recipient}_{j}_species"
            fields={"environment":base, "biomass":base+(b.BIO,),
                    "richness":base+(x,), "both":base+(b.BIO,x)}
            if any(y in f or f.count(x)>1 for f in fields.values()):
                raise ValueError("Target leakage")
            scores.append({k: metrics(test,fit_predict(train,test,f,y,False),y)[0]
                           for k,f in fields.items()})
        m={k:sum(z[k] for z in scores)/2 for k in scores[0]}
        folds.append({"continent":c,"n_test":len(test),"n_train":len(train),
                      "n_test_blocks":len({r["block"] for r in test}),
                      "mse":m,"orientation_mse":scores,
                      "richness_given_biomass":m["biomass"]-m["both"],
                      "biomass_given_richness":m["richness"]-m["both"],
                      "richness_without_biomass":m["environment"]-m["richness"],
                      "biomass_without_richness":m["environment"]-m["biomass"]})
    terms=("richness_given_biomass","biomass_given_richness","richness_without_biomass","biomass_without_richness")
    return {"n_cells":len(rows),"continent_counts":dict(Counter(r["continent"] for r in rows)),
            "donor":donor,"recipient":recipient,"flexible_effort":flex,"AGB_epoch":epoch,
            "ridge_lambda":RIDGE,"n_folds":len(folds),"adequate_folds":len(folds)>=4,
            "folds":folds,"mean_delta_MSE":{t:sum(f[t] for f in folds)/len(folds) for t in terms} if folds else {},
            "n_positive_folds":{t:sum(f[t]>0 for f in folds) for t in terms},
            "relative_richness_gain_given_biomass":1-sum(f["mse"]["both"] for f in folds)/sum(f["mse"]["biomass"] for f in folds) if folds and sum(f["mse"]["biomass"] for f in folds)>0 else None}

def panel(rows):
    rr=b.local_support(rows)
    high=[r for r in rr if min(r["observations"],r["bird_observations"])>=200 and min(r["observers"],r["bird_observers"])>=20]
    return {"primary_flexible":transfer(rr),"linear_effort":transfer(rr,flex=False),
            "AGB2020":transfer(rr,epoch=2020),"high_effort":transfer(high),
            "reverse_flexible":transfer(rr,donor="birds",recipient="plants"),
            "reverse_linear":transfer(rr,donor="birds",recipient="plants",flex=False)}

def main():
    p=argparse.ArgumentParser()
    for n in ("disjoint","agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    prior=b.read(a.disjoint)
    if prior["version"]!="railstatus-observer-disjoint-plant-bird-biomass-v1" or prior["biomass_source_sha256"]!=agb["source_sha256"] or any(prior["source_artifact_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("agb","agb_model","split","annual")):
        raise ValueError("Verified prior artifact mismatch")
    old_panels={(p["width_m"],p["land"]):p for p in prior["panels"]}
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
                previous=old_panels[w,gate]["2025"]
                replay_check(one(by["2025"])[0],previous["primary"])
                replay_check(one(by["2025"],flex=True)[0],previous["flexible_effort"])
                z={"width_m":w,"land":gate,"2025":panel(by["2025"])}
                if w==200000:
                    common=set.intersection(*({r["key"] for r in by[y]} for y in YEARS))
                    common_rows={y:b.local_support([r for r in by[y] if r["key"] in common]) for y in YEARS}
                    if any({r["key"] for r in common_rows[y]}!={r["key"] for r in common_rows["2025"]} for y in YEARS):
                        raise ValueError("Annual support changed")
                    z["common_annual"]={y:panel(common_rows[y]) for y in YEARS}
                panels.append(z);print(json.dumps({"disjoint_plant_bird_biomass_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-observer-disjoint-plant-bird-transfer-v1","prior_disjoint_sha256":hashlib.sha256(a.disjoint.read_bytes()).hexdigest(),"source_replay":checks,"biomass_source_sha256":agb["source_sha256"],"source_artifact_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("agb","agb_model","split","annual")},"panels":panels,
      "primary":"100km land80 2025 leave-continent-out plant-richness incremental MSE gain beyond woody AGB and environment, flexible joint four-half effort, normalized block-weighted ridge, mean two disjoint orientations",
      "limits":["same platform and shared access/detection/environment despite disjoint identities","split-supported adjacent-pair endpoint selection; no representative global claim","heldout folds share training data; no inferential p-values, no causal decomposition","2020/2024 AGB is woody, including zero ocean pixels; land80 primary","earlier annual panels common geography, static map/AGB, separate within-year continent holdout; not a forecast","observed richness, not trophic interactions or bird biomass"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
