#!/usr/bin/env python3
"""Frozen earlier-year test of seed vs invertivore insect-richness contrast."""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_bird_diet_insect import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,ORIENTATIONS,TERMS,EXPECTED,extract,attach,linear,interval,combine,one,temporal,
    SOURCE_SHA,GUILDS,GATE,AXES,diet_attach,DEFF,DFLEX,one_diet,replay,verify_diet,
    extra_attach,fields_for,TARGETS,GUILD_NAMES)
from scan_world_railstatus_annual_insect_cells import candidate_cells
COMPARISONS=(("granivore","invertivore"),("fruit_seed","invertivore"))

def estimate_scores(rows,epoch=2024,flex=False):
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
            if any(v is None for v in z.values()):return dict(qa,status="insufficient support or rank"),None
            z["plant_attenuation"]=combine(z,{"plant_without_insects":1.,"plant_with_insects":-1.})
            gfits[g]=z
        fs[i,j]={t:{GUILD_NAMES[0]:gfits[0][t],GUILD_NAMES[1]:gfits[1][t],
                   GUILD_NAMES[2]:combine({g:gfits[g][t] for g in (0,1)},{0:-1.,1:1.})} for t in TARGETS}
    scores={t:{g:{k:combine({ij:v[t][g] for ij,v in fs.items()},w) for k,w in TERMS.items()} for g in GUILD_NAMES} for t in TARGETS}
    result={t:{g:{k:interval(v) for k,v in q.items()} for g,q in r.items()} for t,r in scores.items()}
    for g in GUILD_NAMES:
        for k in TERMS:
            a,c,d=(result[t][g][k]["estimate"] for t in ("plant_without_insects","plant_with_insects","plant_attenuation"))
            if abs(a-c-d)>1e-9:raise ValueError("Paired attenuation identity failed")
    return dict(qa,status="estimated",effects=result),scores

def verify_previous(new,old):
    if new["status"]!=old["status"] or new["n_cells"]!=old["n_cells"]:raise ValueError("Previous insect-guild model support changed")
    if new["status"]!="estimated":return
    for t,r in new["effects"].items():
        for g,q in r.items():
            for k,v in q.items():
                z=old["effects"][t][g][k]
                if abs(z["estimate"]-v["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(z["interval95_approx"],v["interval95_approx"]))>1e-9:raise ValueError("Previous insect-guild effect replay failed")

def annual_model(native):
    maps={y:{r["key"]:r for r in v} for y,v in native.items()}
    if any(len(maps[y])!=len(native[y]) for y in YEARS):raise ValueError("Repeated cell")
    keys=set.intersection(*(set(v) for v in maps.values()))
    rr={y:b.local_support([maps[y][k] for k in sorted(keys)]) for y in YEARS}
    if any({r["key"] for r in rr[y]}!={r["key"] for r in rr["2025"]} for y in YEARS):raise ValueError("Common cell geography changed")
    for k in keys:
        if len({(maps[y][k]["continent"],maps[y][k]["block"]) for y in YEARS})!=1:raise ValueError("Common block changed")
    out={"n_common_cells":len(rr["2025"]),"n_native_cells":{y:len(b.local_support(v)) for y,v in native.items()},"models":{}}
    for label,epoch,flex in (("linear",2024,False),("flexible",2024,True),("AGB2020",2020,False)):
        m={};f={}
        for y in YEARS:m[y],f[y]=estimate_scores(rr[y],epoch,flex)
        q={"years":m}
        for name,ys in (("earlier_mean",YEARS[:2]),("all_year_mean",YEARS)):
            if not all(f[y] is not None for y in ys):q[name]=None;continue
            q[name]={t:{g:{k:interval(combine({y:f[y][t][g][k] for y in ys},{y:1/len(ys) for y in ys})) for k in TERMS} for g in GUILD_NAMES} for t in TARGETS}
        out["models"][label]=q
    out["earlier_leave_continent_out"]={}
    for c in sorted({r["continent"] for r in rr["2025"]}):
        ff={y:estimate_scores([r for r in rr[y] if r["continent"]!=c],flex=True)[1] for y in YEARS[:2]}
        out["earlier_leave_continent_out"][c]=interval(combine({y:ff[y]["insect_richness"][GUILD_NAMES[2]]["cross"] for y in ff},{y:.5 for y in ff})) if all(v is not None for v in ff.values()) else None
    return out

def compact(z):
    def effects(q):
        if q is None:return None
        return {t:{g:v["cross"] for g,v in r.items()} for t,r in q.items()}
    out={"width_m":z["width_m"],"land":z["land"],"comparisons":{}}
    for name,r in z["comparisons"].items():
        out["comparisons"][name]={"n_common_cells":r["n_common_cells"],"n_native_cells":r["n_native_cells"],"earlier_leave_continent_out":r["earlier_leave_continent_out"],"models":{}}
        for lab,m in r["models"].items():
            out["comparisons"][name]["models"][lab]={"years":{y:{"status":v["status"],"n_cells":v["n_cells"],"cross":effects(v.get("effects"))} for y,v in m["years"].items()},"earlier_mean":effects(m["earlier_mean"]),"all_year_mean":effects(m["all_year_mean"])}
    return out

def main():
    p=argparse.ArgumentParser()
    for n in ("insect_annual","prior_insect","diet","diet_model","disjoint","agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();ia,pi,diet,dm,prior,agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.insect_annual,a.prior_insect,a.diet,a.diet_model,a.disjoint,a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if ia["version"]!="railstatus-annual-insect-disjoint-cells-v1" or not ia["exact_2025_full_and_half_replay"] or ia["width_m"]!=200000 or tuple(ia["years"])!=YEARS or any(ia[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or any(ia["source_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","diet")):raise ValueError("Annual insect source mismatch")
    if pi["version"]!="railstatus-avonet-diet-insect-conditional-coupling-v1" or pi["diet_source"]!=diet["source"] or pi["biomass_source_sha256"]!=agb["source_sha256"] or any(pi["source_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("diet","diet_model","split","disjoint","agb","agb_model")):raise ValueError("Previous insect diet model mismatch")
    insect_panels={p["land"]:p for p in pi["panels"] if p["width_m"]==200000}
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
        for w in (200000,):
            grid=taxa["grids"][str(w)];yr=YEARS
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
                # Verify the 2025 annual insect source is precisely the original
                # eligible insect endpoints on the diet candidate geography.
                cand_insect,expected_insect=candidate_cells(src,diet)
                ia_maps={y:{tuple(r["key"]):(tuple(r["full"]),tuple(tuple(x) for x in r["halves"])) for r in ia["cells"][y]} for y in YEARS}
                if ia_maps["2025"]!=expected_insect or any(not set(ia_maps[y])<=cand_insect[y] for y in YEARS):raise ValueError("Annual insect geography/2025 replay mismatch")
                z={"width_m":w,"land":gate,"comparisons":{}}
                for gs in COMPARISONS:
                    label=gs[0]+"_vs_"+gs[1]
                    native={}
                    for y in YEARS:
                        rr=diet_attach(by[y],records[y],gs)
                        if y=="2025":verify_diet(rr,diet_panels[w,gate]["2025"][label])
                        im={"insects":ia_maps[y]};ic={"insects":{k:v[0] for k,v in ia_maps[y].items()}}
                        native[y]=extra_attach(rr,im,ic,("insects",))
                    # Same native 2025 model, all effect/contrast estimates and intervals.
                    for name,flex in (("primary",False),("flexible_effort",True)):
                        new,_=estimate_scores(native["2025"],flex=flex)
                        verify_previous(new,insect_panels[gate]["2025"][label][name])
                    z["comparisons"][label]=annual_model(native)
                panels.append(z);print(json.dumps({"diet_insect_temporal_panel":compact(z)},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-avonet-diet-insect-temporal-v1","source_replay":checks,"diet_source":diet["source"],"biomass_source_sha256":agb["source_sha256"],
            "source_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("insect_annual","prior_insect","diet","diet_model","split","disjoint","agb","agb_model")},"panels":panels,
            "primary":"200km land80 fixed 2023–25 common cells, earlier 2023/24 equal-year mean invertivore minus granivore insect richness slope, flexible joint effort, opposite observer sides; paired spatial block covariance. Primary before reading earlier insect data. Pooled fruit/seed comparison retained as secondary.",
            "limits":["retrospective years held back for insect-guild contrast, not independent ecosystems",
                      "earlier plant/bird diet data used previously; only insect-guild contrast new",
                      "2025 selected geography, same observation platform, shared access and detection confounding",
                      "annual insect counts not actual prey abundance or diet interactions",
                      "static biomass and rail map; not causal or rail-specific","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
