#!/usr/bin/env python3
"""One fixed earlier-period height comparison; retain direct and fixed-total estimands."""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter
from analyze_world_railstatus_verified_plant_size import (
    b,GEO,CORE,CLIMATE,EXTRA,SHARE,TOTAL,linear,combine,interval,extract,climate_extra,
    GROUPS,WIDTHS,TRAITS,LIFE_SHA,LEP_SHA,SIZE_SHA,SOURCE_SHA,attach,estimate,replay)
from scan_world_railstatus_plant_size_temporal import YEARS,frozen_candidates,compare
HEIGHT="height_nonwoody"

def check_effects(new,old):
    if new["status"]!=old["status"] or new["n_cells"]!=old["n_cells"]:raise ValueError("Previous height support changed")
    if new["status"]=="estimated":
        for k,v in new["effects"].items():
            q=old["effects"][k]
            if abs(v["estimate"]-q["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],q["interval95_approx"]))>1e-9:raise ValueError("Previous height effect changed")

def common(by):
    maps={y:{r["key"]:r for r in by[y]} for y in YEARS}
    if any(len(maps[y])!=len(by[y]) for y in YEARS):raise ValueError("Repeated cells")
    keys=set.intersection(*(set(v) for v in maps.values()))
    for k in keys:
        if len({(maps[y][k]["continent"],maps[y][k]["block"]) for y in YEARS})!=1:raise ValueError("Block drift")
    return {y:b.local_support([maps[y][k] for k in sorted(keys)]) for y in YEARS}

def temporal(by,flex=True,share=False,epoch=2024,both=False):
    rr=common({y:[r for r in by[y] if not both or min(r["active_normal_km"],r["inactive_normal_km"])>0] for y in YEARS})
    fs={};years={}
    for y in YEARS:years[y],fs[y]=estimate(rr[y],flex,share,epoch)
    out={"n_common_cells":len(rr["2025"]),"years":years,"earlier_mean":None}
    if all(fs[y] is not None for y in YEARS[:2]):
        out["earlier_mean"]={k:interval(combine({y:fs[y][k] for y in YEARS[:2]},{y:.5 for y in YEARS[:2]})) for k in fs[YEARS[0]]}
    return out

def ready_year(ready,source):
    out=[]
    for r in ready:
        k=r["key"]
        if k not in source:continue
        o,b0,s=source[k]["full"]
        out.append(dict(r,observations=o,observers=b0,species=s,
                        log_observations=math.log1p(o),log_observers=math.log1p(b0),log_species=math.log1p(s)))
    return out

def evaluate_height(ready,maps,gate,old):
    for label,flex,share in (("primary",True,False),("linear_effort",False,False),("fixed_total_active_share",True,True)):
        q,_=estimate(attach(ready,maps["2025"],HEIGHT,gate),flex,share);check_effects(q,old[label])
    env={y:ready_year(ready,maps[y]) for y in YEARS}
    by={y:attach(env[y],maps[y],HEIGHT,gate) for y in YEARS}
    q={"native":{y:len(b.local_support(v)) for y,v in by.items()},
       "primary":temporal(by),"linear":temporal(by,False),"fixed_total":temporal(by,share=True),
       "both_positive":temporal(by,both=True),"both_positive_fixed_total":temporal(by,share=True,both=True),
       "AGB2020":temporal(by,epoch=2020),
       "family_centered":temporal({y:attach(env[y],maps[y],HEIGHT,gate,family=True) for y in YEARS}),
       "leave_continent_out":{c:temporal({y:[r for r in by[y] if r["continent"]!=c] for y in YEARS}) for c in sorted({r["continent"] for r in common(by)["2025"]})}}
    # Both observer halves evaluated on identical all-year common cells; no
    # averaging of independent standard errors, shared spatial scores retained.
    halves={(y,i):attach(env[y],maps[y],HEIGHT,gate,side=i) for y in YEARS for i in (0,1)}
    keys=set.intersection(*({r["key"] for r in v} for v in halves.values()))
    ff={};hh={}
    for k,v in halves.items():
        hh[k],ff[k]=estimate([r for r in v if r["key"] in keys])
    out={"n_common_cells":len(b.local_support([r for r in halves["2025",0] if r["key"] in keys])),"earlier_mean":None,
         "years":{y:{str(i):hh[y,i] for i in (0,1)} for y in YEARS}}
    early=[(y,i) for y in YEARS[:2] for i in (0,1)]
    if all(ff[k] is not None for k in early):
        out["earlier_mean"]={t:interval(combine({k:ff[k][t] for k in early},{k:.25 for k in early})) for t in ff[early[0]]}
    q["observer_halves"]=out
    return q

def main():
    p=argparse.ArgumentParser()
    for k in ("traits","previous_cells","previous","annual_size","annual","size_model","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src,oc,prev,az,ann,pm,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(a.traits,a.previous_cells,a.previous,a.annual_size,a.annual,a.size_model,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    if az["version"]!="railstatus-annual-verified-plant-size-v1" or not az["exact_annual_counts_and_2025_moments_replay"] or az["source"]!=src["source"] or any(az["source_sha256"][k]!=hashlib.sha256(p.read_bytes()).hexdigest() for k,p in (("split",a.split),("annual",a.annual),("size",a.traits))) or any(az[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong annual height source")
    if pm["version"]!="railstatus-verified-scaled-plant-size-association-v1" or pm["trait_cells_sha256"]!=hashlib.sha256(a.traits.read_bytes()).hexdigest() or pm["source"]!=src["source"]:raise ValueError("Prior height model changed")
    height_panels={(p["land"],p["coverage"]):p for p in pm["panels"] if p["width"]==200000}
    annual_maps={y:{tuple(r["key"]):r for r in az["grids"][y]} for y in YEARS}
    expected=frozen_candidates(split,ann)
    for y in YEARS:
        if {k:(tuple(v["full"]),tuple(tuple(x) for x in v["halves"])) for k,v in annual_maps[y].items()}!=expected[y]:raise ValueError("Annual support/counts changed")
    if src["version"]!="railstatus-verified-scaled-plant-size-cells-v1" or not src["exact_full_split_target_and_life_replay"] or any(src[k]!=split[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or src["source_SHA256"]!={"LifeForm":LIFE_SHA,"size_archive":SIZE_SHA} or src["prior_cells_sha256"]!=hashlib.sha256(a.previous_cells.read_bytes()).hexdigest():raise ValueError("Wrong verified size source")
    if prev["version"]!="railstatus-disturbance-trait-association-v1" or prev["trait_cells_sha256"]!=src["prior_cells_sha256"] or oc["version"]!="railstatus-disturbance-trait-cells-v1":raise ValueError("Prior model source mismatch")
    old_panels={(p["width"],p["land"],p["coverage"]):p for p in prev["panels"] if p["taxon"]=="plants"}
    if any(prior[k]!=src[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Original census mismatch")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1":raise ValueError("Wrong AGB")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen map/climate/pop mismatch")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(ras)
        for w in (200000,):
            grid=prior["grids"][str(w)];keys={b.key(r) for g in GROUPS for r in grid[g]}
            clim,_=b.source_climate_lookup(land,ras,scale,keys,w);ext,qa=climate_extra(a.climate,ras,clim,w)
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            for g in ('plants',):
                tm={tuple(r["key"]):r for r in src["grids"][str(w)][g]}
                old=extract(split["pairs"][str(w)][g])
                if {k:(tuple(v["full"]),tuple(tuple(x) for x in v["halves"])) for k,v in tm.items()}!=old:raise ValueError("Exact split counts not replayed")
                rows,_,_=b.original_cell_rows([b.convert(v) for v in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                rows,_=b.inject_built(rows,built[str(w)]);rows,_=b.inject_status(rows,status,w)
                for landgate in (.8,.5):
                    rr,_=b.enrich_cells(rows,rd,tr,w,landgate);rr,_=b.join_population(rr,pp,rd)
                    ready=[]
                    for r in rr:
                        k=r["key"]
                        if k not in tm or k not in ext or k not in ag20 or k not in ag24:continue
                        if tuple(tm[k]["full"])!=(r["observations"],r["observers"],r["species"]):raise ValueError("Exact full census mismatch")
                        av,cv=ag20[k],ag24[k]
                        if min(av["valid_area_m2"],cv["valid_area_m2"])<.9*rd[k]["source_GRIP4_total_land_area_m2"]:continue
                        if any(not math.isfinite(z["agb_Mg_ha"]) or z["agb_Mg_ha"]<0 for z in (av,cv)):raise ValueError("Invalid AGB")
                        ready.append(dict(r,agb2020=av["agb_Mg_ha"],agb2024=cv["agb_Mg_ha"],**ext[k]))
                    if len(ready)!=height_panels[landgate,.5]["n_base"]:raise ValueError("Prior geography drift")
                    for k,v in tm.items():compare(annual_maps["2025"][k],v)
                    for gate in (.5,.25,.8):
                        old=height_panels[landgate,gate]["traits"][HEIGHT]
                        z={"width":w,"land":landgate,"coverage":gate,**evaluate_height(ready,annual_maps,gate,old)}
                        panels.append(z)
                        print(json.dumps({"annual_height_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-annual-height-status-validation-v1","source_replay":checks,
         "annual_size_sha256":hashlib.sha256(a.annual_size.read_bytes()).hexdigest(),"prior_size_model_sha256":hashlib.sha256(a.size_model.read_bytes()).hexdigest(),"panels":panels,
         "primary":"200km land80 coverage50, 2023/24 equal-year mean on fixed all-three-year common cells, active-minus-inactive separate log-length slopes with flexible effort, climate and woody AGB; fixed-total active share companion explicitly required for strong status-specific interpretation.",
         "limits":["retrospective earlier years on fixed 2026 railway map, not independent geography or actual closure event",
                   "mean species trait, not individual height or measured wind/mowing",
                   "conditional on known nonwoody trait support, same platform and access",
                   "no more threshold or period searching if unsupported; no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
