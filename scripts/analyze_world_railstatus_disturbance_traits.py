#!/usr/bin/env python3
"""Direct active/inactive status contrasts for railway-relevant traits.
Whole landscape, species-equal summaries, explicit known-trait coverage.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import GEO,linear,interval,combine,extract
from analyze_world_railstatus_coupling_climate import climate_extra,EXTRA as CLIMATE
from scan_world_railstatus_disturbance_traits import GROUPS,WIDTHS,TRAITS,LIFE_SHA,LEP_SHA,SIZE_SHA,SOURCE_SHA
from analyze_world_railstatus_cross_taxa_population_local import flexible,EXTRA
CORE=("log_observations","log_observers","coverage")
SHARE="active_share";TOTAL="log_total_mapped_normal_rail"

def attach(rows,source,trait,gate,family=False,side="full"):
    out=[]
    for r in rows:
        z=source[r["key"]]
        q=z["target_full"] if side=="full" else z["target_halves"][side]
        t=q["traits"][trait];n=t["family_known"] if family else t["known"];den=t["denominator"]
        obs,observers,species=q["counts"]
        if den>species or n>den or min(obs,observers,species)<0:raise ValueError("Trait denominator invalid")
        minimum=20 if side=="full" else 10
        if n<minimum or den==0 or n/den<gate or observers<(5 if side=="full" else 3):continue
        v=t["family_excess"] if family else t["mean"]
        if v is None or not math.isfinite(v):raise ValueError("Invalid known trait mean")
        a=r["active_normal_km"];d=r["inactive_normal_km"]
        z=dict(r,trait_value=v,coverage=n/den,known_trait_n=n,target_species=species,
               target_observations=obs,target_observers=observers,
               log_observations=math.log1p(obs),log_observers=math.log1p(observers),
               log_target_species=math.log1p(species),
               **{TOTAL:math.log1p(a+d),SHARE:a/(a+d) if a+d>0 else 0.})
        out.append(flexible(z))
    return out

def estimate(rows,flex=True,share=False,epoch=2024,adjust_biomass=True):
    rr=b.with_epoch(b.local_support(rows),epoch)
    fields=GEO+CORE+CLIMATE+(EXTRA if flex else ())+( (b.BIO,) if adjust_biomass else ())
    if share:
        rr=b.local_support([r for r in rr if r["active_normal_km"]+r["inactive_normal_km"]>0])
        fields=tuple(k for k in fields if k not in (b.ACTIVE,b.INACTIVE))+(TOTAL,SHARE)
    if len(set(fields))!=len(fields):raise ValueError("Duplicate model columns")
    blocks=defaultdict(list)
    for r in rr:blocks[r["continent"],r["block"]].append(r)
    terms=(SHARE,) if share else (b.ACTIVE,b.INACTIVE)
    variation={t:sum(len({r[t] for r in v})>=2 for v in blocks.values()) for t in terms}
    qa={"n_cells":len(rr),"n_blocks":len(blocks),"n_varying_blocks":variation,"continents":dict(Counter(r["continent"] for r in rr))}
    if any(n<20 for n in variation.values()):return dict(qa,status="insufficient status variation"),None
    fs={t:linear(rr,fields,"trait_value",t) for t in terms}
    if any(v is None for v in fs.values()):return dict(qa,status="insufficient support or rank"),None
    if share:fs={"active_share_at_fixed_total":fs[SHARE]}
    else:fs={"active":fs[b.ACTIVE],"inactive":fs[b.INACTIVE],"active_minus_inactive":combine(fs,{b.ACTIVE:1.,b.INACTIVE:-1.})}
    return dict(qa,status="estimated",effects={k:interval(v) for k,v in fs.items()}),fs

def evaluate(rows,source,trait,gate):
    rr=attach(rows,source,trait,gate)
    family=attach(rows,source,trait,gate,True)
    half={i:attach(rows,source,trait,gate,side=i) for i in (0,1)}
    keys=set.intersection(*({r["key"] for r in v} for v in half.values()))
    half={i:b.local_support([r for r in v if r["key"] in keys]) for i,v in half.items()}
    if {r["key"] for r in half[0]}!={r["key"] for r in half[1]}:raise ValueError("Half common support differs")
    h={i:estimate(v) for i,v in half.items()}
    hs={"halves":{str(i):v[0] for i,v in h.items()},"paired_mean":None}
    if all(v[1] is not None for v in h.values()):hs["paired_mean"]={k:interval(combine({i:v[1][k] for i,v in h.items()},{0:.5,1:.5})) for k in h[0][1]}
    both=[r for r in rr if min(r["active_normal_km"],r["inactive_normal_km"])>0]
    return {"n_before_local":len(rr),"primary":estimate(rr)[0],"linear_effort":estimate(rr,False)[0],
            "AGB2020":estimate(rr,epoch=2020)[0],"without_AGB":estimate(rr,adjust_biomass=False)[0],
            "fixed_total_active_share":estimate(rr,share=True)[0],
            "both_status_positive":estimate(both)[0],"both_positive_fixed_total":estimate(both,share=True)[0],
            "high_effort":estimate([r for r in rr if r["target_observations"]>=200 and r["target_observers"]>=20])[0],
            "family_centered":estimate(family)[0],"observer_half_retest":hs,
            "leave_continent_out":{c:estimate([r for r in rr if r["continent"]!=c])[0] for c in sorted({r["continent"] for r in rr})}}

def main():
    p=argparse.ArgumentParser()
    for k in ("traits","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(a.traits,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
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
            for g in GROUPS:
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
                        z={"width":w,"land":landgate,"coverage":gate,"taxon":g,"n_base":len(ready),"traits":{f:evaluate(ready,tm,f,gate) for f in TRAITS[g]}}
                        panels.append(z)
                        compact={f:{k:(v[k] if k=="observer_half_retest" else {"status":v[k]["status"],"n_cells":v[k]["n_cells"],"effects":v[k].get("effects")}) for k in ("primary","linear_effort","fixed_total_active_share","both_status_positive","family_centered","observer_half_retest")} for f,v in z["traits"].items()}
                        print(json.dumps({"disturbance_trait_panel":{"width":w,"land":landgate,"coverage":gate,"taxon":g,"traits":compact}},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-disturbance-trait-association-v1","source_replay":checks,"source":src["source"],"trait_cells_sha256":hashlib.sha256(a.traits.read_bytes()).hexdigest(),"panels":panels,
         "primary":"100km land80 coverage50, at least 20 known trait species, species-equal annual fraction among documented non-woody vascular plants, disturbance-affinity fraction among directional documented butterfly traits, grassland vs forest birds; mapped active-minus-inactive slope, flexible effort plus climate/AGB/geographic controls, local block FE.",
         "limits":["whole-landscape occurrences, not trackside microhabitats; mapped status not measured management or train movement",
                   "non-woody life/size traits conditional on documented growth form; missing trait values not assumed absent",
                   "butterfly denominator taxonomic Papilionoidea, not all insects; upper wingspan is size, not flight performance",
                   "family baseline external pool, not phylogenetic or local species availability control",
                   "species mean traits not within-species adaptation; wind tolerance not measured",
                   "multiple exploratory traits, no universal causal claims; no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
