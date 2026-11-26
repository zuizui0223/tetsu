#!/usr/bin/env python3
"""Test functional trait breadth, not mean composition, in active vs inactive
railway landscapes. Reuse exact existing incidence-weighted moments.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter,defaultdict
from analyze_world_railstatus_disturbance_traits import (
    b,GEO,linear,interval,combine,extract,climate_extra,CLIMATE,
    GROUPS,WIDTHS,TRAITS,LIFE_SHA,LEP_SHA,SIZE_SHA,SOURCE_SHA,
    flexible,EXTRA,CORE,SHARE,TOTAL,attach,estimate)
AXES={"insects":("log_wingspan_upper",),"birds":("HWI","log_wing_length","log_mass")}
SCALE={"HWI":10.,"log_wingspan_upper":1.,"log_wing_length":1.,"log_mass":1.}
def attach_variance(rows,source,trait,gate,side="full"):
    out=[]
    for r in attach(rows,source,trait,gate,side=side):
        src=source[r["key"]];q=src["target_full"] if side=="full" else src["target_halves"][side]
        t=q["traits"][trait];n=t["known"];v=t["variance"];m=t["mean"]
        if n<2 or v is None or not math.isfinite(v) or v<0 or not math.isfinite(m):raise ValueError("Invalid trait moment")
        # n/(n-1) is the exact expected equal-size simple-random-subsample
        # sample variance, not a claim that observed species are random draws.
        scale=SCALE[trait];mu=m/scale
        out.append(dict(r,trait_value=v*n/(n-1)/(scale*scale),
                       trait_mean=mu,trait_mean_sq=mu*mu,log_known_trait_n=math.log1p(n)))
    return out

def fit(rows,flex=True,share=False,epoch=2024,mean_adjust=False):
    rr=b.with_epoch(b.local_support(rows),epoch)
    ff=GEO+CORE+CLIMATE+(EXTRA if flex else ())+(b.BIO,"log_known_trait_n")
    if mean_adjust:ff+=("trait_mean","trait_mean_sq")
    if share:
        rr=b.local_support([r for r in rr if r["active_normal_km"]+r["inactive_normal_km"]>0])
        ff=tuple(k for k in ff if k not in (b.ACTIVE,b.INACTIVE))+(TOTAL,SHARE)
    if len(set(ff))!=len(ff):raise ValueError("Duplicate design")
    groups=defaultdict(list)
    for r in rr:groups[r["continent"],r["block"]].append(r)
    terms=(SHARE,) if share else (b.ACTIVE,b.INACTIVE)
    var={t:sum(len({r[t] for r in v})>=2 for v in groups.values()) for t in terms}
    qa={"n_cells":len(rr),"n_blocks":len(groups),"n_varying_blocks":var,
        "continents":dict(Counter(r["continent"] for r in rr))}
    if any(n<20 for n in var.values()):return dict(qa,status="insufficient status variation"),None
    fs={t:linear(rr,ff,"trait_value",t) for t in terms}
    if any(v is None for v in fs.values()):return dict(qa,status="insufficient support or rank"),None
    if share:fs={"active_share_at_fixed_total":fs[SHARE]}
    else:fs={"active":fs[b.ACTIVE],"inactive":fs[b.INACTIVE],"active_minus_inactive":combine(fs,{b.ACTIVE:1.,b.INACTIVE:-1.})}
    return dict(qa,status="estimated",effects={k:interval(v) for k,v in fs.items()}),fs

def replay(rows,source,trait,gate,old):
    rr=attach(rows,source,trait,gate)
    for label,flex,share in (("primary",True,False),("linear_effort",False,False),("fixed_total_active_share",True,True)):
        z,_=estimate(rr,flex,share);q=old[label]
        if z["status"]!=q["status"] or z["n_cells"]!=q["n_cells"]:raise ValueError("Prior mean support changed")
        if z["status"]=="estimated":
            for k,v in z["effects"].items():
                x=q["effects"][k]
                if abs(v["estimate"]-x["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],x["interval95_approx"]))>1e-9:raise ValueError("Prior mean coefficient mismatch")

def evaluate(rows,source,trait,gate,old):
    replay(rows,source,trait,gate,old)
    rr=attach_variance(rows,source,trait,gate)
    half={i:attach_variance(rows,source,trait,gate,i) for i in (0,1)}
    keys=set.intersection(*({r["key"] for r in v} for v in half.values()))
    half={i:b.local_support([r for r in v if r["key"] in keys]) for i,v in half.items()}
    if {r["key"] for r in half[0]}!={r["key"] for r in half[1]}:raise ValueError("Half support drift")
    fs={i:fit(v) for i,v in half.items()}
    halves={"halves":{str(i):v[0] for i,v in fs.items()},"paired_mean":None}
    if all(v[1] is not None for v in fs.values()):
        halves["paired_mean"]={k:interval(combine({i:v[1][k] for i,v in fs.items()},{0:.5,1:.5})) for k in fs[0][1]}
    both=[r for r in rr if min(r["active_normal_km"],r["inactive_normal_km"])>0]
    return {"primary":fit(rr)[0],"linear_effort":fit(rr,False)[0],"mean_adjusted":fit(rr,mean_adjust=True)[0],
        "AGB2020":fit(rr,epoch=2020)[0],"fixed_total_active_share":fit(rr,share=True)[0],
        "both_status_positive":fit(both)[0],"both_positive_fixed_total":fit(both,share=True)[0],
        "high_effort":fit([r for r in rr if r["target_observations"]>=200 and r["target_observers"]>=20])[0],
        "observer_half_retest":halves,"leave_continent_out":{c:fit([r for r in rr if r["continent"]!=c])[0] for c in sorted({r["continent"] for r in rr})}}

def main():
    p=argparse.ArgumentParser()
    for k in ("traits","previous","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src,prev,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(a.traits,a.previous,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    if prev["version"]!="railstatus-disturbance-trait-association-v1" or prev["source"]!=src["source"] or prev["trait_cells_sha256"]!=hashlib.sha256(a.traits.read_bytes()).hexdigest():raise ValueError("Previous trait model mismatch")
    previous={(p["width"],p["land"],p["coverage"],p["taxon"]):p for p in prev["panels"]}
    if src["version"]!="railstatus-disturbance-trait-cells-v1" or not src["exact_full_and_split_replay"] or any(src[k]!=split[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or src["source"]["source_SHA256"]!={"LifeForm":LIFE_SHA,"LepTraits":LEP_SHA,"AVONET":SOURCE_SHA,"size_archive":SIZE_SHA}:raise ValueError("Wrong trait source")
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
            for g in AXES:
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
                    for gate in (.5,.25,.8):
                        p0=previous[w,landgate,gate,g]
                        if len(ready)!=p0["n_base"]:raise ValueError("Prior environment support changed")
                        z={"width":w,"land":landgate,"coverage":gate,"taxon":g,"n_base":len(ready),
                           "traits":{f:evaluate(ready,tm,f,gate,p0["traits"][f]) for f in AXES[g]}}
                        panels.append(z);print(json.dumps({"trait_dispersion_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-trait-dispersion-v1","source_replay":checks,"trait_cells_sha256":hashlib.sha256(a.traits.read_bytes()).hexdigest(),
         "previous_mean_model_sha256":hashlib.sha256(a.previous.read_bytes()).hexdigest(),"panels":panels,
         "primary":"100km land80 coverage50, direct active-minus-inactive slope on species-equal HWI sample variance (HWI/10 units) and butterfly log upper-wingspan sample variance, flexible effort, known trait count and coverage, climate, AGB and local block FE. Mean-adjusted variance and fixed-total active share sensitivity; 200km and land50 secondary.",
         "limits":["incidence-weighted observed species sample variance, not individual phenotypic or phylogenetic diversity",
                   "n/(n-1) removes finite-list variance divisor effect, not sampling/detection bias",
                   "HWI and wing size are morphology proxies, not wind tolerance or measured dispersal",
                   "trait missingness and regional habitat composition remain; no family-conditioned variance possible from these aggregates",
                   "mean-adjusted variance is a descriptive sensitivity, not causal conditioning",
                   "whole landscapes, mapped status not randomized closure or actual management, no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
