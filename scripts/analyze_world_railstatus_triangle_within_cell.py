#!/usr/bin/env python3
"""Fixed annual three-taxon panel: spatial versus within-cell richness association.
Reuse exact previous source and common-cell geography, no new occurrence scan.
"""
import argparse,json,hashlib
from pathlib import Path
from collections import Counter
from analyze_world_railstatus_annual_triangle import (
    b,cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,ORIENTATIONS,TERMS,EXPECTED,extract,attach,linear,interval,combine,
    one,temporal,SOURCE_SHA,replay,extra_attach,fields_for,candidate_cells,
    climate_extra,attach_climate,fit_triangle,verify_native,EXTRA,EDGES,CROSS)
from analyze_world_railstatus_plant_bird_within_cell import fit_panel,common_rows

def panel_fields(flex):
    # Static geographic fields are exactly absorbed by cell fixed effects.
    ff=fields_for(("insects",),flex)
    result=tuple(k for k in ff if k not in GEO+(b.BIO,))
    if len(set(result))!=len(result) or any(not k.startswith(("plants_","birds_","insects_")) for k in result):raise ValueError("Invalid dynamic effort fields")
    return result

def estimate(by,trend="year",flex=True):
    rows=next(iter(by.values()))
    qa={"n_cells":len(rows),"n_years":len(by),"n_blocks":len({(r["continent"],r["block"]) for r in rows})}
    fields=panel_fields(flex)
    effects={}
    for g,h in EDGES:
        fs={}
        for i,j in ORIENTATIONS:
            x=f"{g}_{i}_species";y=f"{h}_{j}_species"
            fs[i,j]=fit_panel(by,fields+(x,),y,x,trend)
        if any(v is None for v in fs.values()):
            effects[g+"__"+h]=None;continue
        effects[g+"__"+h]={k:interval(combine(fs,w)) for k,w in TERMS.items()}
    # Same plant/insect observer side, opposite bird side. Estimate conditional
    # association and attenuation only; no path product or causal mediation.
    fs={}
    for i,j in CROSS:
        p=f"plants_{i}_species";q=f"insects_{i}_species";y=f"birds_{j}_species"
        a=fit_panel(by,fields+(p,),y,p,trend)
        c=fit_panel(by,fields+(p,q),y,p,trend)
        d=fit_panel(by,fields+(p,q),y,q,trend)
        if any(v is None for v in (a,c,d)):fs={};break
        fs[i,j]={"plant_without_insects":a,"plant_with_insects":c,"insects":d,"plant_attenuation":combine({"a":a,"c":c},{"a":1.,"c":-1.})}
    joint={k:interval(combine({ij:v[k] for ij,v in fs.items()},{ij:.5 for ij in CROSS})) for k in fs[CROSS[0]]} if fs else None
    if joint and abs(joint["plant_without_insects"]["estimate"]-joint["plant_with_insects"]["estimate"]-joint["plant_attenuation"]["estimate"])>1e-9:raise ValueError("Paired attenuation identity failed")
    return dict(qa,status="estimated" if all(v is not None for v in effects.values()) and joint is not None else "partial or insufficient support",
                pairwise=effects,joint=joint)

def replay_spatial(rr,old):
    if len(rr["2025"])!=old["n_common_cells"]:raise ValueError("Prior common support changed")
    for label,flex,extra,epoch in (("flexible_climate",True,EXTRA,2024),("linear_climate",False,EXTRA,2024),
                                  ("flexible_original",True,(),2024),("flexible_climate_AGB2020",True,EXTRA,2020)):
        for y in YEARS:
            new,_=fit_triangle(rr[y],flex,extra,epoch);prev=old["models"][label]["years"][y]
            if new["status"]!=prev["status"] or new["n_cells"]!=prev["n_cells"]:raise ValueError("Prior spatial model support changed")
            if new["status"]!="estimated":continue
            if set(new["effects"])!=set(prev["effects"]):raise ValueError("Spatial effect columns changed")
            for k,v in new["effects"].items():
                q=prev["effects"][k]
                if abs(v["estimate"]-q["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],q["interval95_approx"]))>1e-9:raise ValueError("Previous spatial coefficients not replayed")

def analyze(native,old):
    rr=common_rows(native);replay_spatial(rr,old)
    specs={"cell_year_flex":("year",True),"cell_year":("year",False),
           "cell_block_year_flex":("block_year",True),"cell_block_year":("block_year",False)}
    return {"n_common_cells":len(rr["2025"]),"continents":dict(Counter(r["continent"] for r in rr["2025"])),
            "models":{k:estimate(rr,*v) for k,v in specs.items()},
            "periods":{ys[0]+"_"+ys[1]:{k:estimate({y:rr[y] for y in ys},*v) for k,v in specs.items()} for ys in (YEARS[:2],YEARS[1:])},
            "leave_continent_out":{c:estimate({y:[r for r in rr[y] if r["continent"]!=c] for y in YEARS},"year",True) for c in sorted({r["continent"] for r in rr["2025"]})},
            "prior_spatial_earlier_mean":old["models"]["flexible_climate"]["earlier_mean"]}

def main():
    p=argparse.ArgumentParser()
    for n in ("triangle","insect_annual","triad_climate","diet","disjoint","agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();previous,ia,tc,diet,prior,agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.triangle,a.insect_annual,a.triad_climate,a.diet,a.disjoint,a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if previous["version"]!="railstatus-observer-disjoint-annual-climate-triangle-v1" or any(previous["source_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("insect_annual","triad_climate","diet","split","annual","disjoint","agb","agb_model","climate")):raise ValueError("Prior triangle artifact mismatch")
    previous_panels={p["land"]:p for p in previous["panels"] if p["width_m"]==200000}
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
                z={"width_m":w,"land":gate,**analyze(native,previous_panels[gate])}
                panels.append(z);print(json.dumps({"triangle_within_cell_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    out={"version":"railstatus-observer-disjoint-triangle-within-cell-v1","source_replay":checks,
         "triangle_sha256":hashlib.sha256(a.triangle.read_bytes()).hexdigest(),"panels":panels,
         "primary":"200km land80 common 532 cells 2023–25, cell and year FE with flexible effort all six halves. New plant-insect and insect-bird opposite-side slopes; plant-bird prior diagnostic on exact common triad sample. Report direct same-minus-cross contrasts, not comparison of separate p-values. Block-year and two-year subperiods fixed sensitivities.",
         "limits":["within-cell recorded richness changes, not measured population dynamics",
                   "static climate, biomass and rail exposure absorbed, time-varying environment and access may remain",
                   "pairwise disjoint observers only; predictor plant/insect observer half may overlap",
                   "three years, common support selected from 2025 geography, spatial-block intervals approximate",
                   "spatial versus temporal slopes are distinct estimands, no causal mediation or rail effect","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
