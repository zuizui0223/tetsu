#!/usr/bin/env python3
"""Direct common-unit mapped rail and implied road length comparison of height."""
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

ROAD_KM="log_implied_total_land_road_km"
ROAD_VALID="log_implied_valid_land_road_km"
DENSITY="log_road_density"
def road_lengths(r,source):
    land=float(source["source_GRIP4_total_land_area_m2"])
    valid=float(source["source_GRIP4_road_coverage_of_source_land"])
    density=float(source["source_GRIP4_observed_road_density_m_per_km2_land"])
    km=float(source["source_GRIP4_implied_road_km_in_valid_land"])
    if not all(math.isfinite(v) for v in (land,valid,density,km)) or land<=0 or not .95<=valid<=1.0000001 or min(density,km)<0:
        raise ValueError("Invalid source road coverage or implied length")
    whole=density*land/1e9
    if not math.isclose(km,whole*valid,rel_tol=1e-7,abs_tol=1e-7):
        raise ValueError("GRIP4 road density-area identity failed")
    if not math.isclose(r[DENSITY],math.log1p(density),abs_tol=1e-10):
        raise ValueError("Original road density not replayed")
    return dict(r,**{ROAD_KM:math.log1p(whole),ROAD_VALID:math.log1p(km)})

def fit_roads(rows,road_term=ROAD_KM,flex=True,epoch=2024):
    rr=b.with_epoch(b.local_support(rows),epoch)
    fields=GEO+CORE+CLIMATE+(EXTRA if flex else ())+(b.BIO,)
    if fields.count(DENSITY)!=1 or len(fields)!=len(set(fields)):raise ValueError("Original road term missing")
    fields=tuple(road_term if f==DENSITY else f for f in fields)
    terms=(b.ACTIVE,b.INACTIVE,road_term)
    if len(fields)!=len(set(fields)) or not set(terms)<=set(fields):raise ValueError("Invalid road comparison design")
    groups={}
    for r in rr:groups.setdefault((r["continent"],r["block"]),[]).append(r)
    variation={t:sum(len({r[t] for r in v})>=2 for v in groups.values()) for t in terms}
    qa={"n_cells":len(rr),"n_blocks":len(groups),"n_varying_blocks":variation,"continents":dict(Counter(r["continent"] for r in rr))}
    if any(n<20 for n in variation.values()):return dict(qa,status="insufficient exposure variation"),None
    fs={t:linear(rr,fields,"trait_value",t) for t in terms}
    if any(v is None for v in fs.values()):return dict(qa,status="insufficient support or rank"),None
    A,I,R=terms
    fs={"active":fs[A],"inactive":fs[I],"road":fs[R],
        "active_minus_inactive":combine(fs,{A:1.,I:-1.}),
        "active_minus_road":combine(fs,{A:1.,R:-1.}),
        "inactive_minus_road":combine(fs,{I:1.,R:-1.})}
    return dict(qa,status="estimated",effects={k:interval(v) for k,v in fs.items()}),fs

def temporal_roads(by,road_term=ROAD_KM,flex=True,epoch=2024):
    rr=common(by)
    ids=lambda x:{r["key"] for r in x}
    if any(ids(rr[y])!=ids(rr["2025"]) for y in YEARS):raise ValueError("Annual support differs")
    ys={};fs={}
    for y in YEARS:ys[y],fs[y]=fit_roads(rr[y],road_term,flex,epoch)
    out={"n_common_cells":len(rr["2025"]),"years":ys,"earlier_mean":None}
    if all(fs[y] is not None for y in YEARS[:2]):
        out["earlier_mean"]={k:interval(combine({y:fs[y][k] for y in YEARS[:2]},{y:.5 for y in YEARS[:2]})) for k in fs[YEARS[0]]}
    return out

def verify_temporal(new,old):
    if new["n_common_cells"]!=old["n_common_cells"]:raise ValueError("Prior annual common geography changed")
    for y in YEARS:check_effects(new["years"][y],old["years"][y])
    u,v=new["earlier_mean"],old["earlier_mean"]
    if (u is None)!=(v is None):raise ValueError("Prior annual support changed")
    if u is not None:
        for k in v:
            if abs(u[k]["estimate"]-v[k]["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(u[k]["interval95_approx"],v[k]["interval95_approx"]))>1e-9:raise ValueError("Prior earlier coefficient/interval changed")

def evaluate_road_height(ready,maps,gate,old,old_annual):
    for label,flex,share in (("primary",True,False),("linear_effort",False,False),("fixed_total_active_share",True,True)):
        q,_=estimate(attach(ready,maps["2025"],HEIGHT,gate),flex,share);check_effects(q,old[label])
    env={y:ready_year(ready,maps[y]) for y in YEARS}
    by={y:attach(env[y],maps[y],HEIGHT,gate) for y in YEARS}
    original=temporal(by)
    verify_temporal(original,old_annual["primary"])
    q={"native":{y:len(b.local_support(v)) for y,v in by.items()},
       "original_density_replay":original,
       "primary_common_km_scale":temporal_roads(by),
       "observed_valid_land_km":temporal_roads(by,ROAD_VALID),
       "linear_effort":temporal_roads(by,flex=False),
       "AGB2020":temporal_roads(by,epoch=2020),
       "family_centered":temporal_roads({y:attach(env[y],maps[y],HEIGHT,gate,family=True) for y in YEARS}),
       "leave_continent_out":{c:temporal_roads({y:[r for r in by[y] if r["continent"]!=c] for y in YEARS}) for c in sorted({r["continent"] for r in common(by)["2025"]})}}
    # Same cells in all years and both observer sides, paired block scores.
    halves={(y,i):attach(env[y],maps[y],HEIGHT,gate,side=i) for y in YEARS for i in (0,1)}
    keys=set.intersection(*({r["key"] for r in v} for v in halves.values()))
    ff={};hh={}
    for k,v in halves.items():
        hh[k],ff[k]=fit_roads([r for r in v if r["key"] in keys])
    z={"n_common_cells":len(b.local_support([r for r in halves["2025",0] if r["key"] in keys])),
       "earlier_mean":None,"years":{y:{str(i):hh[y,i] for i in (0,1)} for y in YEARS}}
    early=[(y,i) for y in YEARS[:2] for i in (0,1)]
    if all(ff[k] is not None for k in early):
        z["earlier_mean"]={t:interval(combine({k:ff[k][t] for k in early},{k:.25 for k in early})) for t in ff[early[0]]}
    q["observer_halves"]=z
    return q

def main():
    p=argparse.ArgumentParser()
    for k in ("traits","previous_cells","previous","annual_size","annual","size_model","annual_model","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src,oc,prev,az,ann,pm,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(a.traits,a.previous_cells,a.previous,a.annual_size,a.annual,a.size_model,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    am=b.read(a.annual_model)
    if am["version"]!="railstatus-annual-height-status-validation-v1" or am["annual_size_sha256"]!=hashlib.sha256(a.annual_size.read_bytes()).hexdigest() or am["prior_size_model_sha256"]!=hashlib.sha256(a.size_model.read_bytes()).hexdigest():raise ValueError("Prior annual height artifact differs")
    prior_annual={(p["land"],p["coverage"]):p for p in am["panels"]}
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
                        ready.append(road_lengths(dict(r,agb2020=av["agb_Mg_ha"],agb2024=cv["agb_Mg_ha"],**ext[k]),rd[k]))
                    if len(ready)!=height_panels[landgate,.5]["n_base"]:raise ValueError("Prior geography drift")
                    for k,v in tm.items():compare(annual_maps["2025"][k],v)
                    for gate in (.5,):
                        old=height_panels[landgate,gate]["traits"][HEIGHT]
                        z={"width":w,"land":landgate,"coverage":gate,**evaluate_road_height(ready,annual_maps,gate,old,prior_annual[landgate,gate])}
                        panels.append(z)
                        print(json.dumps({"annual_height_road_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-annual-height-road-contrast-v1","source_replay":checks,
         "annual_size_sha256":hashlib.sha256(a.annual_size.read_bytes()).hexdigest(),"prior_size_model_sha256":hashlib.sha256(a.size_model.read_bytes()).hexdigest(),"panels":panels,
         "primary":"200km land80 coverage50, 2023/24 equal-year mean on exact prior common height cells; joint log1p mapped active rail km, inactive rail km and GRIP4 implied total-land road km; direct active-minus-road slope contrast. Active-minus-inactive companion.",
         "limits":["road km inferred from original gridded road density, >=95% valid land coverage and uniform-density extrapolation over missing fraction; not exact road vector lengths",
                   "same units and log1p transform do not imply a physical equal-km replacement or equal causal effect; different maps and reference periods",
                   "exploratory follow-up of prior trait results on same observation platform and coarse geography, not independent replication",
                   "species mean height, not actual local stature or mowing/wind exposure; no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
