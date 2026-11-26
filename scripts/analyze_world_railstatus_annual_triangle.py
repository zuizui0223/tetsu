#!/usr/bin/env python3
"""Common-cell plant–insect–bird spatial triangle, with earlier-year replay.

Each edge uses opposite observer halves; there are not three mutually
disjoint observer groups. All three half-effort sets and climate controls
are identical across the three edges. Bird conditional slopes are a distinct
diagnostic; do not interpret them as causal mediation.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_bird_diet_insect_temporal import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,ORIENTATIONS,TERMS,EXPECTED,extract,attach,linear,interval,combine,
    one,temporal,SOURCE_SHA,replay,extra_attach,fields_for,candidate_cells)
from analyze_world_railstatus_coupling_climate import climate_extra,attach_climate,fit_model,check,EXTRA
EDGES=(("plants","insects"),("plants","birds"),("insects","birds"))
CROSS=((0,1),(1,0))
KEYS=tuple(g+"__"+h for g,h in EDGES)
COVS=("insects","plant_with_added_richness","plant_slope_attenuation")

def fit_triangle(rows,flex=True,extra=EXTRA,epoch=2024):
    rr=b.with_epoch(b.local_support(rows),epoch)
    fields=fields_for(("insects",),flex)+extra
    qa={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
        "continents":dict(Counter(r["continent"] for r in rr))}
    scores={}
    for g,h in EDGES:
        fs={}
        for i,j in ORIENTATIONS:
            x=f"{g}_{i}_species";y=f"{h}_{j}_species"
            fs[i,j]=linear(rr,fields+(x,),y,x)
        if any(v is None for v in fs.values()):return dict(qa,status="insufficient support or rank"),None
        for lab,weights in TERMS.items():
            scores[g+"__"+h+"__"+lab]=combine(fs,weights)
    joint,joint_scores=fit_model(rr,("insects",),flex,extra)
    # fit_model currently uses 2024 woody AGB; only use joint results for 2024.
    if epoch==2024:
        if joint_scores is None:return dict(qa,status="insufficient joint support or rank"),None
        for k in COVS:scores["bird_joint__"+k]=joint_scores[k]
    return dict(qa,status="estimated",effects={k:interval(v) for k,v in scores.items()}),scores

def verify_native(rows,old):
    for label,flex in (("linear_effort",False),("flexible_effort",True)):
        q,_=fit_model(rows,("insects",),flex,EXTRA)
        check(q,old["models"][label]["extra_climate"])

def annual_panel(by):
    maps={y:{r["key"]:r for r in v} for y,v in by.items()}
    if any(len(maps[y])!=len(by[y]) for y in YEARS):raise ValueError("Repeated cell")
    keys=set.intersection(*(set(v) for v in maps.values()))
    rr={y:b.local_support([maps[y][k] for k in sorted(keys)]) for y in YEARS}
    if any({r["key"] for r in rr[y]}!={r["key"] for r in rr["2025"]} for y in YEARS):raise ValueError("Common support drift")
    for k in keys:
        if len({(maps[y][k]["continent"],maps[y][k]["block"]) for y in YEARS})!=1:raise ValueError("Common block drift")
    out={"n_common_cells":len(rr["2025"]),"n_native_cells":{y:len(b.local_support(v)) for y,v in by.items()},"models":{}}
    for label,flex,extra,epoch in (("flexible_climate",True,EXTRA,2024),("linear_climate",False,EXTRA,2024),
                                  ("flexible_original",True,(),2024),("flexible_climate_AGB2020",True,EXTRA,2020)):
        yrs={};fs={}
        for y in YEARS:yrs[y],fs[y]=fit_triangle(rr[y],flex,extra,epoch)
        z={"years":yrs}
        for name,ys in (("earlier_mean",YEARS[:2]),("all_year_mean",YEARS)):
            if all(fs[y] is not None for y in ys):
                z[name]={k:interval(combine({y:fs[y][k] for y in ys},{y:1/len(ys) for y in ys})) for k in fs[ys[0]]}
            else:z[name]=None
        out["models"][label]=z
    out["earlier_leave_continent_out"]={}
    for c in sorted({r["continent"] for r in rr["2025"]}):
        fs={y:fit_triangle([r for r in rr[y] if r["continent"]!=c])[1] for y in YEARS[:2]}
        out["earlier_leave_continent_out"][c]={k:interval(combine({y:fs[y][k] for y in fs},{y:.5 for y in fs})) for k in fs[YEARS[0]] if k.endswith("__cross") or k.startswith("bird_joint__")} if all(v is not None for v in fs.values()) else None
    return out

def compact(z):
    out={"land":z["land"],"width_m":z["width_m"],"n_common_cells":z["n_common_cells"],"n_native_cells":z["n_native_cells"],"models":{}}
    def f(q):return {k:v for k,v in q.items() if k.endswith("__cross") or k.startswith("bird_joint__")} if q is not None else None
    for lab,r in z["models"].items():
        out["models"][lab]={"years":{y:{"status":v["status"],"n_cells":v["n_cells"],"effects":f(v.get("effects"))} for y,v in r["years"].items()},"earlier_mean":f(r["earlier_mean"]),"all_year_mean":f(r["all_year_mean"])}
    out["earlier_leave_continent_out"]=z["earlier_leave_continent_out"]
    return out

def main():
    p=argparse.ArgumentParser()
    for n in ("insect_annual","triad_climate","diet","disjoint","agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();ia,tc,diet,prior,agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.insect_annual,a.triad_climate,a.diet,a.disjoint,a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if ia["version"]!="railstatus-annual-insect-disjoint-cells-v1" or not ia["exact_2025_full_and_half_replay"] or ia["width_m"]!=200000 or tuple(ia["years"])!=YEARS or any(ia[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or any(ia["source_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","diet")):raise ValueError("Annual insect source mismatch")
    if tc["version"]!="railstatus-observer-disjoint-richness-extra-climate-v1" or tc["worldclim_sha256"]!=b.sha_file(a.climate):raise ValueError("Prior climate result mismatch")
    climate_panels={p["land"]:p for p in tc["panels"] if p["width_m"]==200000}
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
        for w in (200000,):
            grid=taxa["grids"][str(w)];yr=YEARS
            ps={y:{g:src["pairs"][str(w)][g] if y=="2025" else annual["pairs"][y][g] for g in GROUPS} for y in yr}
            maps={y:{g:extract(ps[y][g]) for g in GROUPS} for y in yr}
            prior_counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in GROUPS}
            for g in GROUPS:
                if any(prior_counts[g].get(k)!=v[0] for k,v in maps["2025"][g].items()):raise ValueError("2025 exact source replay failed")
            keys={b.key(z) for g in GROUPS for z in grid[g]}|{k for y in yr for g in GROUPS for k in maps[y][g]}
            clim,_=b.source_climate_lookup(land,raster,scale,keys,w)
            ext,qa=climate_extra(a.climate,raster,clim,w)
            print(json.dumps({"triangle_climate_source":qa}),flush=True)
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
                insect_map={"insects":extract(src["pairs"][str(w)]["insects"])}
                insect_full={"insects":{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid["insects"]}}
                if any(insect_full["insects"].get(k)!=v[0] for k,v in insect_map["insects"].items()):raise ValueError("Full insect census replay failed")
                # Verify the 2025 annual insect source is precisely the original
                # eligible insect endpoints on the diet candidate geography.
                cand_insect,expected_insect=candidate_cells(src,diet)
                ia_maps={y:{tuple(r["key"]):(tuple(r["full"]),tuple(tuple(x) for x in r["halves"])) for r in ia["cells"][y]} for y in YEARS}
                if ia_maps["2025"]!=expected_insect or any(not set(ia_maps[y])<=cand_insect[y] for y in YEARS):raise ValueError("Annual insect geography/2025 replay mismatch")
                native={}
                for y in YEARS:
                    im={"insects":ia_maps[y]};ic={"insects":{k:v[0] for k,v in ia_maps[y].items()}}
                    native[y]=attach_climate(extra_attach(by[y],im,ic,("insects",)),ext)
                verify_native(native["2025"],climate_panels[gate]["insect"])
                z={"width_m":w,"land":gate,**annual_panel(native)}
                panels.append(z);print(json.dumps({"annual_triangle_panel":compact(z)},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    out={"version":"railstatus-observer-disjoint-annual-climate-triangle-v1","source_replay":checks,
         "source_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("insect_annual","triad_climate","diet","split","annual","disjoint","agb","agb_model","climate")},"panels":panels,
         "primary":"200km land80, fixed 2023-25 common cells, earlier 2023/24 mean plant-richness association with insect-richness response, opposite observer sides, flexible joint effort across all six halves, BIO1 latitude and BIO4/BIO12/BIO15 plus AGB and geography/rail controls, paired spatial block covariance.",
         "limits":["pairwise disjoint observers, not three mutually disjoint groups","2025 candidate geography, same original iNat platform and access patterns",
                   "cross-sectional slopes, not within-cell population changes or causal trophic mediation",
                   "grid-center climate normals and static 2026 railway map, woody AGB not prey or bird biomass",
                   "no rail-specific effect established, no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
