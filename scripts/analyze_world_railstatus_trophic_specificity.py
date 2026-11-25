#!/usr/bin/env python3
"""Frozen matched resource-specificity contrast using existing diet and insect aggregates."""
import argparse,json,math,hashlib
from collections import Counter
from pathlib import Path
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,EXPECTED,extract,attach,linear,interval,combine,one,temporal)
from analyze_world_railstatus_plant_bird_diet import (
    diet_attach,DEFF,DFLEX,AXES,replay,one_diet)
from analyze_world_railstatus_plant_insect_bird import extra_attach,fields_for,model
from scan_world_railstatus_bird_diet_cells import SOURCE_SHA
CROSS=((0,1),(1,0))
TERMS=("insect_invertivore_minus_fruit_seed","plant_fruit_seed_minus_invertivore","insect_invertivore","insect_fruit_seed","plant_invertivore","plant_fruit_seed")

def matched(rows,epoch=2024,flex=True):
    rr=b.with_epoch(b.local_support(rows),epoch)
    fields=fields_for(("insects",),flex)+DEFF+(DFLEX if flex else ())
    if len(set(fields))!=len(fields):raise ValueError("Duplicate predictors")
    qa={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
        "continents":dict(Counter(r["continent"] for r in rr))}
    out={}
    for i,j in CROSS:
        plant=f"plants_{i}_species";insect=f"insects_{i}_species"
        pred=fields+(plant,insect)
        p={g:linear(rr,pred,f"d{g}_{j}_species",plant) for g in (0,1)}
        q={g:linear(rr,pred,f"d{g}_{j}_species",insect) for g in (0,1)}
        if any(z is None for z in (*p.values(),*q.values())):
            return dict(qa,status="insufficient support or rank")
        out[i,j]={"insect_invertivore_minus_fruit_seed":combine(q,{1:1.,0:-1.}),
                  "plant_fruit_seed_minus_invertivore":combine(p,{0:1.,1:-1.}),
                  "insect_invertivore":q[1],"insect_fruit_seed":q[0],
                  "plant_invertivore":p[1],"plant_fruit_seed":p[0]}
    effects={t:interval(combine({ij:z[t] for ij,z in out.items()},{ij:.5 for ij in CROSS})) for t in TERMS}
    for contrast,positive,negative in ((TERMS[0],TERMS[2],TERMS[3]),(TERMS[1],TERMS[5],TERMS[4])):
        if abs(effects[contrast]["estimate"]-effects[positive]["estimate"]+effects[negative]["estimate"])>1e-10:
            raise ValueError("Paired contrast identity failed")
    return dict(qa,status="estimated",effects=effects)

def panel(rows,records,insect_map,insect_full):
    rr=b.local_support(extra_attach(diet_attach(rows,records,AXES[0]),{"insects":insect_map},{"insects":insect_full},("insects",)))
    high=[r for r in rr if min(r["observations"],r["bird_observations"],r["insects_full_observations"])>=200 and min(r["observers"],r["bird_observers"],r["insects_full_observers"])>=20]
    return {"n_supported":len(rr),"primary_flexible":matched(rr),
            "linear_effort":matched(rr,flex=False),"AGB2020":matched(rr,epoch=2020),
            "high_effort":matched(high),
            "leave_continent_out":{c:matched([r for r in rr if r["continent"]!=c]) for c in sorted({r["continent"] for r in rr})}}

def compare_replay(new,old):
    if new["status"]!=old["status"] or new["n_cells"]!=old["n_cells"]:
        raise ValueError("Previous model support changed")
    if new["status"]!="estimated":return
    for k,v in new["effects"].items():
        z=old["effects"][k]
        if abs(v["estimate"]-z["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],z["interval95_approx"]))>1e-9:
            raise ValueError("Previous coefficient or interval failed replay")

def main():
    p=argparse.ArgumentParser()
    for n in ("diet","diet_model","insect_model","disjoint","agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();diet,prior,agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.diet,a.disjoint,a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    diet_model=b.read(a.diet_model);insect_model=b.read(a.insect_model)
    if diet_model["version"]!="railstatus-avonet-diet-plant-biomass-coupling-v1" or diet_model["diet_source"]!=diet["source"] or diet_model["biomass_source_sha256"]!=agb["source_sha256"]:
        raise ValueError("Frozen diet model source mismatch")
    if insect_model["version"]!="railstatus-observer-disjoint-plant-insect-bird-biomass-v1" or insect_model["prior_disjoint_SHA256"]!=hashlib.sha256(a.disjoint.read_bytes()).hexdigest() or insect_model["biomass_source_sha256"]!=agb["source_sha256"]:
        raise ValueError("Frozen insect model source mismatch")
    diet_panels={(p["width_m"],p["land"]):p for p in diet_model["panels"]}
    insect_panels={(p["width_m"],p["land"]):p for p in insect_model["panels"]}
    if diet['version']!='railstatus-avonet-trophic-disjoint-cell-counts-v1' or not diet['exact_full_and_split_cell_replay'] or diet['source']['source_SHA256']!=SOURCE_SHA or any(diet[k]!=src[k] for k in ('obs_sha256','taxa_sha256','split_seed')):raise ValueError('Wrong trophic source')
    if prior['version']!='railstatus-observer-disjoint-plant-bird-biomass-v1' or prior['biomass_source_sha256']!=agb['source_sha256'] or any(prior['source_artifact_sha256'][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ('agb','agb_model','split','annual')):raise ValueError('Wrong previous disjoint model')
    prior_panels={(p['width_m'],p['land']):p for p in prior['panels']}
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
            insect_map=extract(src["pairs"][str(w)]["insects"])
            insect_full={b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid["insects"]}
            if any(insect_full.get(k)!=v[0] for k,v in insect_map.items()):
                raise ValueError("Exact insect full count replay failed")
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
                replay(one(by["2025"])[0],oldp["2025"]["primary"])
                replay(one(by["2025"],flex=True)[0],oldp["2025"]["flexible_effort"])
                records={y:{tuple(r["key"]):r for r in diet["grids"][str(w)][y]} for y in yr}
                r=by["2025"];record=records["2025"]
                diet_rows=b.local_support(diet_attach(r,record,AXES[0]))
                insect_rows=b.local_support(extra_attach(r,{"insects":insect_map},{"insects":insect_full},("insects",)))
                old_d=diet_panels[w,gate]["2025"]["fruit_seed_vs_invertivore"]
                old_i=insect_panels[w,gate]["insect"]
                for flex,label in ((False,"primary"),(True,"flexible_effort")):
                    compare_replay(one_diet(diet_rows,flex=flex)[0],old_d[label])
                    compare_replay(model(insect_rows,("insects",),flex=flex),old_i[label])
                z={"width_m":w,"land":gate,"2025":panel(r,record,insect_map,insect_full)}
                panels.append(z);print(json.dumps({"trophic_specificity_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-observer-disjoint-trophic-specificity-v1","source_replay":checks,
        "input_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("diet","diet_model","insect_model","disjoint","split")},
        "panels":panels,
        "primary":"2025 200km land80 flexible effort: direct invertivore minus fruit_seed coefficient on log insect richness, controlling opposite-side plant richness, woody AGB, mapped rail and other geography, all full and guild efforts. Mean disjoint orientations, paired block covariance.",
        "limits":["same iNaturalist platform; plant and insect predictors share observer half, opposite to birds",
                  "all insects are not measured prey; AVONET dietary category is not realized diet or trophic flux",
                  "guild habitat/taxonomy/detectability confounding remains; not a causal effect",
                  "2025 only, no independent replication; no eBird or railway-specific modification test"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
