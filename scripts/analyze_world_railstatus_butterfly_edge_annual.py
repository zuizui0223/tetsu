#!/usr/bin/env python3
"""Frozen held-back year retest of the exploratory butterfly edge-affinity signal."""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import extract,interval,combine
from analyze_world_railstatus_coupling_climate import climate_extra
from analyze_world_railstatus_disturbance_traits import attach,estimate
from scan_world_railstatus_butterfly_edge_annual import YEARS,G,TRAIT,WIDTHS,LEP_SHA,SALT,validate_close,full_gate

def compare_years(native,epoch=2024,flex=True,share=False):
    maps={y:{r["key"]:r for r in native[y]} for y in YEARS}
    if any(len(maps[y])!=len(native[y]) for y in YEARS):raise ValueError("Duplicate annual cell")
    keys=set.intersection(*(set(m) for m in maps.values()))
    by={y:b.local_support([maps[y][k] for k in sorted(keys)]) for y in YEARS}
    if any({r["key"] for r in by[y]}!={r["key"] for r in by["2025"]} for y in YEARS):raise ValueError("Annual common geography changed")
    q={"n_native":{y:len(b.local_support(native[y])) for y in YEARS},"n_common":len(by["2025"]),"years":{},"earlier_mean":None,"all_year_mean":None}
    scores={}
    for y in YEARS:q["years"][y],scores[y]=estimate(by[y],epoch=epoch,flex=flex,share=share)
    for label,ys in (("earlier_mean",YEARS[:2]),("all_year_mean",YEARS)):
        if all(scores[y] is not None for y in ys):
            q[label]={t:interval(combine({y:scores[y][t] for y in ys},{y:1/len(ys) for y in ys})) for t in scores[ys[0]]}
    if q["earlier_mean"] is not None and scores["2025"] is not None:
        q["2025_minus_earlier"]={t:interval(combine({y:scores[y][t] for y in YEARS},{"2023":-.5,"2024":-.5,"2025":1.})) for t in scores["2025"]}
    return q,by

def compare_halves(base,maps,gate):
    native={}
    for y in YEARS:
        native[y]={i:attach([r for r in base if r["key"] in maps[y]],maps[y],TRAIT,gate,side=i) for i in (0,1)}
    m={(y,i):{r["key"]:r for r in native[y][i]} for y in YEARS for i in (0,1)}
    keys=set.intersection(*(set(q) for q in m.values()))
    rr={k:b.local_support([q[x] for x in sorted(keys)]) for k,q in m.items()}
    if any({r["key"] for r in q}!={r["key"] for r in rr["2025",0]} for q in rr.values()):raise ValueError("Common observer-half geography changed")
    f={};scores={};out={"n_common":len(rr["2025",0]),"years":{},"earlier_mean":None}
    for y in YEARS:
        for i in (0,1):f[y,i],scores[y,i]=estimate(rr[y,i])
        if all(scores[y,i] is not None for i in (0,1)):
            out["years"][y]={t:interval(combine({i:scores[y,i][t] for i in (0,1)},{0:.5,1:.5})) for t in scores[y,0]}
        else:out["years"][y]=None
    ks=[(y,i) for y in YEARS[:2] for i in (0,1)]
    if all(scores[k] is not None for k in ks):out["earlier_mean"]={t:interval(combine({k:scores[k][t] for k in ks},{k:.25 for k in ks})) for t in scores[ks[0]]}
    return out

def verify_model(a,z):
    if a["n_cells"]!=z["n_cells"] or a["status"]!=z["status"]:raise ValueError("Prior trait model support changed")
    if a["status"]=="estimated":
        for k,v in a["effects"].items():validate_close(v,z["effects"][k])

def evaluate(base,maps,gate,prior):
    native={y:attach([r for r in base if r["key"] in maps[y]],maps[y],TRAIT,gate) for y in YEARS}
    verify_model(estimate(native["2025"])[0],prior["primary"])
    q,common=compare_years(native)
    family={y:attach([r for r in base if r["key"] in maps[y]],maps[y],TRAIT,gate,True) for y in YEARS}
    verify_model(estimate(family["2025"])[0],prior["family_centered"])
    out={"primary":q,"linear_effort":compare_years(native,flex=False)[0],
         "AGB2020":compare_years(native,epoch=2020)[0],
         "family_centered":compare_years(family)[0],
         "fixed_total_active_share":compare_years(native,share=True)[0],
         "observer_half_retest":compare_halves(base,maps,gate)}
    out["earlier_leave_continent_out"]={}
    for c in sorted({r["continent"] for r in common["2025"]}):
        out["earlier_leave_continent_out"][c]=compare_years({y:[r for r in common[y] if r["continent"]!=c] for y in YEARS})[0]["earlier_mean"]
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("annual_edge","traits","prior_model","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();annual,src,prev,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(a.annual_edge,a.traits,a.prior_model,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    if annual["version"]!="railstatus-butterfly-edge-annual-cells-v1" or not annual["exact_2025_full_half_trait_replay"] or annual["LepTraits_SHA256"]!=LEP_SHA or annual["years"]!=list(YEARS) or any(annual["source_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("traits","split")):raise ValueError("Annual edge source invalid")
    if any(len({annual[k],src[k],split[k]})!=1 for k in ("obs_sha256","taxa_sha256","split_seed")) or any(prior[k]!=src[k] for k in ("obs_sha256","taxa_sha256")) or src["split_seed"]!=SALT.decode():raise ValueError("Original source mismatch")
    if prev["version"]!="railstatus-disturbance-trait-association-v1" or prev["trait_cells_sha256"]!=hashlib.sha256(a.traits.read_bytes()).hexdigest() or prev["source"]!=src["source"]:raise ValueError("Wrong previous model")
    prior_panels={(v["width"],v["land"],v["coverage"]):v["traits"][TRAIT] for v in prev["panels"] if v["taxon"]==G}
    if not status["version"].startswith("railstatus-v1") or agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1":raise ValueError("Wrong map or AGB source")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen environmental sources changed")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(ras)
        for w in WIDTHS:
            grid=prior["grids"][str(w)]
            maps={y:{tuple(r["key"]):r for r in annual["grids"][str(w)][y]} for y in YEARS}
            old=extract(split["pairs"][str(w)][G]);oldtr={tuple(r["key"]):r for r in src["grids"][str(w)][G]}
            if {k:(tuple(r["full"]),tuple(tuple(x) for x in r["halves"])) for k,r in maps["2025"].items()}!=old:raise ValueError("Exact 2025 split replay failed")
            for y in YEARS:
                if len(maps[y])!=len(annual["grids"][str(w)][y]) or not set(maps[y])<=set(old) or any(not full_gate(r["full"],r["halves"]) for r in maps[y].values()):raise ValueError("Annual geography/gate changed")
            for k,r in maps["2025"].items():
                for x,z in zip([r["target_full"],*r["target_halves"]],[oldtr[k]["target_full"],*oldtr[k]["target_halves"]]):
                    validate_close(x["counts"],z["counts"]);validate_close(x["traits"][TRAIT],z["traits"][TRAIT])
            keys={b.key(r) for g in ("plants","insects","birds") for r in grid[g]}
            clim,_=b.source_climate_lookup(land,ras,scale,keys,w);ext,_=climate_extra(a.climate,ras,clim,w)
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            rows,_,_=b.original_cell_rows([b.convert(v) for v in grid[G]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
            rows,_=b.inject_built(rows,built[str(w)]);rows,_=b.inject_status(rows,status,w)
            for landgate in (.8,.5):
                rr,_=b.enrich_cells(rows,rd,tr,w,landgate);rr,_=b.join_population(rr,pp,rd)
                ready=[]
                for r in rr:
                    k=r["key"]
                    if k not in old or k not in ext or k not in ag20 or k not in ag24:continue
                    if tuple(oldtr[k]["full"])!=(r["observations"],r["observers"],r["species"]):raise ValueError("Original cell census changed")
                    av,cv=ag20[k],ag24[k]
                    if min(av["valid_area_m2"],cv["valid_area_m2"])<.9*rd[k]["source_GRIP4_total_land_area_m2"]:continue
                    if any(not math.isfinite(z["agb_Mg_ha"]) or z["agb_Mg_ha"]<0 for z in (av,cv)):raise ValueError("Invalid AGB")
                    ready.append(dict(r,agb2020=av["agb_Mg_ha"],agb2024=cv["agb_Mg_ha"],**ext[k]))
                for gate in (.25,.5,.8):
                    q=evaluate(ready,maps,gate,prior_panels[w,landgate,gate])
                    z={"width":w,"land":landgate,"coverage":gate,"trait":TRAIT,**q}
                    panels.append(z);print(json.dumps({"butterfly_edge_annual_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-butterfly-edge-annual-association-v1","source_replay":checks,
         "source_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("annual_edge","traits","prior_model","split","climate")},
         "primary":"100km land80 coverage25, 2023–24 equal-year mean active-minus-inactive edge-affinity contrast in fixed 2023–25 common cells, prior flexible observation/environment model; discovery secondary remains exploratory.",
         "panels":panels,"limits":["earlier years same source and fixed 2025 geography, not independent ecosystems or railway history","coverage25 incomplete, known directional trait subset only; stricter coverage50/80 not relaxed","mapped rail status not actual operation, habitat preference not observed rail use or adaptation","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
