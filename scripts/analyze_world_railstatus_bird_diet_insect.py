#!/usr/bin/env python3
"""Direct food-guild contrast in conditional insect richness association.
No new observation scan. Same joint cells/design and opposite observer side.
"""
import argparse,json,math,hashlib
from collections import Counter
from pathlib import Path
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_diet import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,ORIENTATIONS,TERMS,EXPECTED,extract,attach,linear,interval,combine,one,temporal,
    SOURCE_SHA,GUILDS,GATE,AXES,diet_attach,DEFF,DFLEX,one_diet,replay)
from analyze_world_railstatus_plant_insect_bird import extra_attach,fields_for
CROSS=((0,1),(1,0))
TARGETS=("plant_without_insects","plant_with_insects","insect_richness","woody_AGB","plant_attenuation")
GUILD_NAMES=("fruit_seed_or_granivore","invertivore","invertivore_minus_comparator")

def verify_diet(rows,old):
    for name,flex in (("primary",False),("flexible_effort",True)):
        now,_=one_diet(rows,flex=flex);prev=old[name]
        if now["status"]!=prev["status"] or now["n_cells"]!=prev["n_cells"]:raise ValueError("Diet support changed")
        if now["status"]=="estimated":
            for m,g in now["effects"].items():
                for k,orient in g.items():
                    for lab,z in orient.items():
                        q=prev["effects"][m][k][lab]
                        if abs(z["estimate"]-q["estimate"])>1e-9 or max(abs(x-y) for x,y in zip(z["interval95_approx"],q["interval95_approx"]))>1e-9:
                            raise ValueError("Prior diet coefficient/interval mismatch")

def estimate(rows,epoch=2024,flex=False):
    rr=b.with_epoch(b.local_support(rows),epoch)
    fields=fields_for(("insects",),flex)+DEFF+(DFLEX if flex else ())
    if len(set(fields))!=len(fields):raise ValueError("Duplicate effort columns")
    qa={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
        "continents":dict(Counter(r["continent"] for r in rr))}
    fs={}
    for i,j in ORIENTATIONS:
        plant=f"plants_{i}_species";insect=f"insects_{i}_species"
        base=fields+(plant,);full=base+(insect,)
        gfits={}
        for g in (0,1):
            target=f"d{g}_{j}_species"
            z={"plant_without_insects":linear(rr,base,target,plant),
               "plant_with_insects":linear(rr,full,target,plant),
               "insect_richness":linear(rr,full,target,insect),
               "woody_AGB":linear(rr,full,target,b.BIO)}
            if any(v is None for v in z.values()):return dict(qa,status="insufficient support or rank")
            z["plant_attenuation"]=combine(z,{"plant_without_insects":1.,"plant_with_insects":-1.})
            gfits[g]=z
        fs[i,j]={t:{GUILD_NAMES[0]:gfits[0][t],GUILD_NAMES[1]:gfits[1][t],
                   GUILD_NAMES[2]:combine({g:gfits[g][t] for g in (0,1)},{0:-1.,1:1.})} for t in TARGETS}
    result={t:{g:{k:interval(combine({ij:v[t][g] for ij,v in fs.items()},w)) for k,w in TERMS.items()} for g in GUILD_NAMES} for t in TARGETS}
    for g in GUILD_NAMES:
        for k in TERMS:
            a,c,d=(result[t][g][k]["estimate"] for t in ("plant_without_insects","plant_with_insects","plant_attenuation"))
            if abs(a-c-d)>1e-9:raise ValueError("Paired attenuation identity failed")
    return dict(qa,status="estimated",effects=result)

def panel(rows,gs):
    rr=b.local_support(rows)
    return {"groups":list(gs),"n_supported":len(rr),"primary":estimate(rr),"flexible_effort":estimate(rr,flex=True),
         "AGB2020_same_cells":estimate(rr,2020),
         "high_effort":estimate([r for r in rr if min(r["observations"],r["bird_observations"],r["insects_full_observations"])>=200 and min(r["observers"],r["bird_observers"],r["insects_full_observers"])>=20]),
         "leave_continent_out":{c:estimate([r for r in rr if r["continent"]!=c]) for c in sorted({r["continent"] for r in rr})}}

def main():
    p=argparse.ArgumentParser()
    for n in ("diet","diet_model","disjoint","agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();diet,dm,prior,agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.diet,a.diet_model,a.disjoint,a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if diet['version']!='railstatus-avonet-trophic-disjoint-cell-counts-v1' or not diet['exact_full_and_split_cell_replay'] or diet['source']['source_SHA256']!=SOURCE_SHA or any(diet[k]!=src[k] for k in ('obs_sha256','taxa_sha256','split_seed')):raise ValueError('Wrong trophic source')
    if prior['version']!='railstatus-observer-disjoint-plant-bird-biomass-v1' or prior['biomass_source_sha256']!=agb['source_sha256'] or any(prior['source_artifact_sha256'][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ('agb','agb_model','split','annual')):raise ValueError('Wrong previous disjoint model')
    prior_panels={(p['width_m'],p['land']):p for p in prior['panels']}
    if dm['version']!='railstatus-avonet-diet-plant-biomass-coupling-v1' or dm['diet_source']!=diet['source'] or dm['biomass_source_sha256']!=agb['source_sha256']:raise ValueError('Previous diet model mismatch')
    diet_panels={(p['width_m'],p['land']):p for p in dm['panels']}
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
                oldp=prior_panels[w,gate]
                replay(one(by["2025"])[0],oldp["2025"]["primary"])
                replay(one(by["2025"],flex=True)[0],oldp["2025"]["flexible_effort"])
                records={y:{tuple(r["key"]):r for r in diet["grids"][str(w)][y]} for y in yr}
                insect_map={"insects":extract(src["pairs"][str(w)]["insects"])}
                insect_full={"insects":{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid["insects"]}}
                if any(insect_full["insects"].get(k)!=v[0] for k,v in insect_map["insects"].items()):raise ValueError("Full insect census replay failed")
                z={"width_m":w,"land":gate,"2025":{}}
                # Primary pooled fruit/seed vs invertivore; separately report
                # granivore and frugivore without changing their support gates.
                for gs in AXES:
                    k=gs[0]+"_vs_"+gs[1]
                    rr=diet_attach(by["2025"],records["2025"],gs)
                    verify_diet(rr,diet_panels[w,gate]["2025"][k])
                    z["2025"][k]=panel(extra_attach(rr,insect_map,insect_full,("insects",)),gs)
                panels.append(z)
                def compact(q):
                    if q["status"]!="estimated":return q
                    return {"status":q["status"],"n_cells":q["n_cells"],"n_blocks":q["n_blocks"],"cross":{m:{g:v["cross"] for g,v in z.items()} for m,z in q["effects"].items()}}
                print(json.dumps({"diet_insect_panel":{"width_m":w,"land":gate,"2025":{k:{j:(compact(v) if isinstance(v,dict) and "status" in v else {c:compact(q) for c,q in v.items()} if j=="leave_continent_out" else v) for j,v in r.items()} for k,r in z["2025"].items()}}},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-avonet-diet-insect-conditional-coupling-v1","source_replay":checks,"diet_source":diet["source"],"biomass_source_sha256":agb["source_sha256"],
         "source_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("diet","diet_model","split","disjoint","agb","agb_model")},"panels":panels,
         "primary":"200km land80 2025, direct invertivore minus pooled fruit/seed bird insect-richness slope conditional on plant richness, woody AGB, same geography and joint full-taxon/guild/half efforts. Both opposite observer orientations averaged with paired block scores; flexible effort mandatory robustness, granivore/frugivore separately secondary.",
         "limits":["all recorded insect richness, not prey diversity/abundance, AVONET guild not realized diet",
                   "plant and insect predictor observers may overlap; bird opposite side cannot share identities",
                   "guild differences may reflect habitat, range size, pool size, detectability or saturation",
                   "cross-sectional associations and slope attenuation, not causal mediation",
                   "annual insect aggregates absent; no earlier-period insect contrast",
                   "woody AGB not productivity or bird biomass, railway specificity not demonstrated","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
