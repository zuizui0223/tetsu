#!/usr/bin/env python3
"""Mapped status contrasts for author-scaled plant size traits."""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter,defaultdict
from analyze_world_railstatus_disturbance_traits import (
    b,GEO,CORE,CLIMATE,EXTRA,SHARE,TOTAL,linear,combine,interval,extract,climate_extra,
    GROUPS,WIDTHS,TRAITS,LIFE_SHA,LEP_SHA,SIZE_SHA,SOURCE_SHA,attach,estimate,evaluate)
from scan_world_railstatus_verified_plant_size import AXES,OLD,check_old
def replay(rows,src,old,gate):
    for f in OLD:
        rr=attach(rows,src,f,gate)
        for label,flex,share in (("primary",True,False),("linear_effort",False,False),("fixed_total_active_share",True,True)):
            z,_=estimate(rr,flex,share);q=old["traits"][f][label]
            if z["n_cells"]!=q["n_cells"] or z["status"]!=q["status"]:raise ValueError("Previous mean support changed")
            if z["status"]=="estimated":
                for k,v in z["effects"].items():
                    x=q["effects"][k]
                    if abs(v["estimate"]-x["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],x["interval95_approx"]))>1e-9:raise ValueError("Previous mean coefficients changed")

def life_adjusted(rows,src,f,gate):
    rr=attach(rows,src,f,gate);a={r["key"]:r for r in attach(rows,src,"annual_nonwoody",gate)}
    common=b.with_epoch(b.local_support([dict(r,annual_fraction=a[r["key"]]["trait_value"]) for r in rr if r["key"] in a]),2024)
    out={"n_matched_cells":len(common),"matched_baseline":estimate(common)[0]}
    fields=GEO+CORE+CLIMATE+EXTRA+(b.BIO,"annual_fraction")
    fs={t:linear(common,fields,"trait_value",t) for t in (b.ACTIVE,b.INACTIVE)}
    out["life_adjusted"]=None
    if all(v is not None for v in fs.values()):
        out["life_adjusted"]={k:interval(v) for k,v in {"active":fs[b.ACTIVE],"inactive":fs[b.INACTIVE],"active_minus_inactive":combine(fs,{b.ACTIVE:1.,b.INACTIVE:-1.})}.items()}
    return out
def main():
    p=argparse.ArgumentParser()
    for k in ("traits","previous_cells","previous","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src,oc,prev,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(a.traits,a.previous_cells,a.previous,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
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
        for w in WIDTHS:
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
                    om={tuple(v["key"]):v for v in oc["grids"][str(w)]["plants"]}
                    if set(om)!=set(tm):raise ValueError("Prior plant source support mismatch")
                    for k,r in tm.items():
                        for x,y in zip([r["target_full"],*r["target_halves"]],[om[k]["target_full"],*om[k]["target_halves"]]):
                            if x["counts"]!=y["counts"]:raise ValueError("Target counts changed")
                            check_old(x["traits"],y["traits"])
                    for gate in (.5,.25,.8):
                        pp0=old_panels[w,landgate,gate]
                        if len(ready)!=pp0["n_base"]:raise ValueError("Prior geography drift")
                        replay(ready,tm,pp0,gate)
                        z={"width":w,"land":landgate,"coverage":gate,"taxon":"plants","n_base":len(ready),
                           "traits":{f:{**evaluate(ready,tm,f,gate),"life_adjusted":life_adjusted(ready,tm,f,gate)} for f in AXES}}
                        panels.append(z);print(json.dumps({"verified_plant_size_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-verified-scaled-plant-size-association-v1","source_replay":checks,
         "source":src["source"],"trait_cells_sha256":hashlib.sha256(a.traits.read_bytes()).hexdigest(),"panels":panels,
         "primary":"100km land80 coverage50, at least 20 known nonwoody species; direct active-minus-inactive log rail length coefficient, author-scaled log plant height. Seed mass and stem density secondary; prior frozen environment and flexible effort, known trait coverage, local block FE. Fixed total, observer-half, external family and life-cycle composition matched sensitivity.",
         "limits":["author log-transformed and scaled values; no physical-unit inverse transformation",
                   "species averages, not individual size change, wind resistance or measured seed dispersal",
                   "conditional on known nonwoody traits; missing species not assigned values",
                   "same whole-landscape occurrence source, mapped railway status not closure experiment, no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
