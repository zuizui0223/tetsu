#!/usr/bin/env python3
"""Same-cell direct rail-status contrasts after fixed observation sampling."""
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import GEO,linear,interval,combine
from analyze_world_railstatus_observation_support import make_frame,attach_agb,classify
from analyze_world_railstatus_coupling_climate import climate_extra,EXTRA as CLIMATE
from analyze_world_railstatus_joint_climate_selection import SQUARE_BASE
from analyze_world_rail_status_richness_railstatus import CATS,LOG
from scan_world_railstatus_fixed_observer_richness import candidates,counts,GROUPS,K,MIN_RECORDS,REPS,SEED

BASE=GEO+(b.BIO,)+CLIMATE+tuple("extra_sq_"+k for k in SQUARE_BASE)
MODES=("observer_balanced","record_balanced","full")
STATUS=(b.ACTIVE,b.INACTIVE)

def prepare(rows,data,gs,ext,ag):
    out=[]
    for r in rows:
        k=r["key"]
        if k not in ext or k not in ag or not all(k in data[g] for g in gs):continue
        z=dict(r,**ext[k]);z[b.BIO]=ag[k][2024]
        z["rail_class"]=classify(r)
        for x in SQUARE_BASE:z["extra_sq_"+x]=r[x]**2
        for g in gs:
            q=data[g][k]
            for mode in MODES:
                value=math.log1p(q["full"][2]) if mode=="full" else q[mode]["mean_log1p_species"]
                if not math.isfinite(value) or value<0 or (mode!="full" and not math.log(2)-1e-12<=value<=math.log1p(K)+1e-12):
                    raise ValueError("Invalid fixed effort outcome")
                z[g+"_"+mode]=value
            o,h=[math.log1p(n) for n in q["full"][:2]]
            for term,v in (("records",o),("observers",h)):
                z[g+"_"+term]=v;z[g+"_"+term+"_sq"]=v*v;z[g+"_"+term+"_cu"]=v**3
            z[g+"_effort_interaction"]=o*h
        out.append(z)
    return b.local_support(out)

def fit(rows,gs,mode,extra_effort=False):
    rr=b.local_support(rows);blocks=defaultdict(list)
    for r in rr:blocks[r["continent"],r["block"]].append(r)
    effort=tuple(g+"_"+k for g in gs for k in ("records","observers","records_sq","observers_sq","records_cu","observers_cu","effort_interaction"))
    fields=BASE+(effort if extra_effort else ())
    qa={"n_cells":len(rr),"n_blocks":len(blocks),"continents":dict(Counter(r["continent"] for r in rr)),
        "rail_class_counts":dict(Counter(r["rail_class"] for r in rr)),
        "varying_status_blocks":{k:sum(len({r[k] for r in v})>1 for v in blocks.values()) for k in STATUS}}
    if min(qa["varying_status_blocks"].values())<20:return dict(qa,status="insufficient status support"),None
    fits={}
    for g in gs:
        z={k:linear(rr,fields,g+"_"+mode,k) for k in STATUS}
        if any(v is None for v in z.values()):return dict(qa,status="insufficient support or rank"),None
        fits[g]={"active":z[b.ACTIVE],"inactive":z[b.INACTIVE],"active_minus_inactive":combine(z,{b.ACTIVE:1.,b.INACTIVE:-1.})}
    fits["plant_minus_comparator"]={k:combine({g:fits[g][k] for g in gs},{gs[0]:1.,gs[1]:-1.}) for k in fits[gs[0]]}
    return dict(qa,status="estimated",effects={g:{k:interval(v) for k,v in q.items()} for g,q in fits.items()}),fits

def panel(rows,gs):
    out={"groups":list(gs),"models":{}}
    for lab,mode,adj in (("fixed_observers","observer_balanced",False),("fixed_records","record_balanced",False),
                        ("fixed_observers_effort_adjusted","observer_balanced",True),("full_counts_effort_adjusted","full",True)):
        out["models"][lab]=fit(rows,gs,mode,adj)[0]
    # Sensitivity to how the sample is standardized on identical cells and the same design.
    _,a=fit(rows,gs,"observer_balanced");_,z=fit(rows,gs,"record_balanced")
    if a is not None and z is not None:
        out["observer_minus_record"]={g:{k:interval(combine({"observer":a[g][k],"record":z[g][k]},{"observer":1.,"record":-1.})) for k in a[g]} for g in a}
    out["leave_continent_out"]={c:fit([r for r in rows if r["continent"]!=c],gs,"observer_balanced")[0] for c in sorted({r["continent"] for r in rows})}
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("fixed","prior","plant","rail","status","edges","built","tree","road","population","agb","agb_model","ne","climate","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();sample,prior,plant,rail,status,edge,tree,road,pop,agb,am=map(b.read,(a.fixed,a.prior,a.plant,a.rail,a.status,a.edges,a.tree,a.road,a.population,a.agb,a.agb_model))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    if sample["version"]!="railstatus-fixed-observer-richness-v1" or not sample["exact_candidate_full_counts_replayed"] or (sample["K"],sample["min_full_records"],sample["reps"],sample["seed"])!=(K,MIN_RECORDS,REPS,SEED) or sample["year"]!=prior["year"]:raise ValueError("Wrong fixed effort protocol")
    if sample["source_sha256"]!={"prior":hashlib.sha256(a.prior.read_bytes()).hexdigest(),"observations":prior["obs_sha256"],"taxa":prior["taxa_sha256"]}:raise ValueError("Fixed effort source mismatch")
    if not status["version"].startswith("railstatus-v1") or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Rail census mismatch")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or am["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or am["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass mismatch")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen environmental mismatch")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    cand=candidates(prior);panels=[]
    try:
        scale,_=b.climate_scale(ras)
        for w in b.GRID_WIDTHS:
            data={g:{tuple(r["key"]):r for r in sample["grids"][str(w)][g]} for g in GROUPS}
            for g in GROUPS:
                if len(data[g])!=len(sample["grids"][str(w)][g]) or {k:tuple(v["full"]) for k,v in data[g].items()}!=cand[w][g]:
                    raise ValueError("Fixed effort candidate replay failed")
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            clim,_=b.source_climate_lookup(land,ras,scale,set(rd),w)
            ext,meta=climate_extra(a.climate,ras,clim,w)
            frame=make_frame(clim,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)])
            frame,_=b.inject_built(frame,built[str(w)]);frame,_=b.inject_status(frame,status,w)
            for gate in b.LAND_THRESHOLDS:
                rr,lqa=b.enrich_cells(frame,rd,tr,w,gate);rr,pqa=b.join_population(rr,pp,rd)
                ag=attach_agb(rr,ag20,ag24,rd)
                expected={(100000,.5):13194,(100000,.8):12138,(200000,.5):3263,(200000,.8):2874}
                if len(rr)!=expected[w,gate]:raise ValueError("Environmental frame count changed")
                z={"width_m":w,"land":gate,"n_environmental_frame":len(rr),"climate_source":meta}
                for g in ("birds","insects"):
                    gs=("plants",g);m=prepare(rr,data,gs,ext,ag)
                    z[g+"_comparison"]=panel(m,gs)
                panels.append(z);print(json.dumps({"fixed_observer_rail_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-fixed-observer-richness-rail-contrast-v1","source_replay":checks,
         "source_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("fixed","prior","road","population","status","climate")},
         "primary":"100km land80 plant-vs-bird same cells, plant active-minus-inactive log-km slope contrast on mean log1p observed richness of 20 distinct observers with one record each; local block FE and climate/environment, no raw effort covariates; plant-minus-bird double contrast secondary.",
         "panels":panels,"limits":["observation-standardized recorded richness not true richness, neither observers nor records are random field surveys",
         "observer identities can overlap taxa; not the previous disjoint-observer coupling design",
         "equal effort gate still selects well-sampled cells; no causal rail effect, train dispersal or legacy conclusion",
         "intervals conditional on finite fixed sampling draws; full counts different estimand","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
