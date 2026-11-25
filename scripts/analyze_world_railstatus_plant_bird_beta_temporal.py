#!/usr/bin/env python3
"""Retrospective temporal validation of plant-bird community similarity.
Use 2023-24 records on the fixed 2025 candidate geography; pair covariance.
No new observations and no temporal window/threshold selection.
"""
import argparse,json,hashlib
from pathlib import Path
import analyze_world_railstatus_plant_insect_bird_beta as b
import analyze_world_railstatus_observer_split_connectivity as s
from analyze_world_railstatus_observer_annual_connectivity import cells_from_pairs
from scan_world_railstatus_observer_annual_beta import candidate_support,pair_key
from analyze_world_railstatus_plant_bird_disjoint import combine,interval
W=200000
YEARS=("2023","2024","2025")
GROUPS=("plants","birds")
SPECS={"block_Jaccard":("jaccard",False,False,2024),
       "flexible_block_Jaccard":("jaccard",False,True,2024),
       "endpoint_Jaccard":("jaccard",True,False,2024),
       "block_Simpson":("simpson",False,False,2024),
       "block_Jaccard_AGB2020":("jaccard",False,False,2020)}

def estimate(rows,metric,endpoint,flex,epoch):
    r=b.setup(rows,GROUPS,epoch,flex);fields=b.fields(GROUPS,endpoint,flex)
    fs={}
    for i,j in b.ORIENT:
        x=f"plants_{i}_{metric}";y=f"birds_{j}_{metric}"
        z=b.score_fit(r,fields+(x,),y,(x,),endpoint)
        if z is None:return {"status":"insufficient support or rank","n_pairs":len(r)},None
        fs[i,j]=z[x]
    f=combine(fs,{ij:.5 for ij in b.ORIENT})
    return {"status":"estimated","n_pairs":len(r),"effect":interval(f),
            "orientations":{f"{i}{j}":interval(v) for (i,j),v in fs.items()}},f

def validate_previous(rows,old):
    for name,(metric,endpoint,flex,epoch) in SPECS.items():
        prior_key={"block_Jaccard":"primary","flexible_block_Jaccard":"flexible_effort",
                   "endpoint_Jaccard":"endpoint_FE","block_Simpson":"simpson",
                   "block_Jaccard_AGB2020":"AGB2020"}[name]
        new,_=estimate(rows,metric,endpoint,flex,epoch);prev=old[prior_key]
        if new["n_pairs"]!=prev["n_pairs"] or new["status"]!=prev["status"]:raise ValueError("Previous pair support/status changed")
        if new["status"]=="estimated":
            q=prev["effects"]["plant_joint"];v=new["effect"]
            if abs(q["estimate"]-v["estimate"])>1e-9 or max(abs(x-y) for x,y in zip(q["interval95_approx"],v["interval95_approx"]))>1e-9:
                raise ValueError("Previous pair coefficient or interval changed")

def analyze(by):
    def index(rows):
        z={(r["a"],r["b"]):r for r in rows}
        if len(z)!=len(rows):raise ValueError("Repeated pair")
        return z
    maps={y:index(rows) for y,rows in by.items()}
    keys=set.intersection(*(set(v) for v in maps.values()))
    common={y:[maps[y][k] for k in sorted(keys)] for y in YEARS}
    for k in keys:
        if len({(maps[y][k]["continent"],maps[y][k]["block"]) for y in YEARS})!=1:raise ValueError("Annual geography mismatch")
    out={"n_common_pairs":len(keys),"n_native_pairs":{y:len(v) for y,v in by.items()},
         "n_common_blocks":len({(r["continent"],r["block"]) for r in common["2025"]}),
         "n_common_by_continent":{},"models":{}}
    from collections import Counter
    out["n_common_by_continent"]=dict(Counter(r["continent"] for r in common["2025"]))
    for name,spec in SPECS.items():
        native={y:estimate(by[y],*spec)[0] for y in YEARS}
        annual={};fits={}
        for y in YEARS:annual[y],fits[y]=estimate(common[y],*spec)
        z={"native":native,"common":annual}
        for label,ys in (("earlier_2023_2024",YEARS[:2]),("all_year_mean",YEARS)):
            if all(fits[y] is not None for y in ys):
                z[label]=interval(combine({y:fits[y] for y in ys},{y:1/len(ys) for y in ys}))
            else:z[label]=None
        out["models"][name]=z
    # Spatial influence is a diagnostic, not a new model-selection step.
    loo={}
    for c in sorted(out["n_common_by_continent"]):
        q={y:estimate([r for r in common[y] if r["continent"]!=c],"jaccard",False,False,2024)[1] for y in YEARS[:2]}
        loo[c]=interval(combine(q,{y:.5 for y in q})) if all(v is not None for v in q.values()) else None
    out["earlier_block_Jaccard_leave_continent_out"]=loo
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("annual","beta","split","prior","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    annual,previous,src,prior,edge,agb,am,plant,rail,status,tree,road,pop=map(s.read,(a.annual,a.beta,a.split,a.prior,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=s.verify_original_sources(plant,rail);s.verify_replay(prior,plant)
    if previous["version"]!="railstatus-disjoint-community-composition-coupling-v1" or previous["biomass_source_SHA256"]!=agb["source_sha256"] or any(previous["source_SHA256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","edges","agb","agb_model")):raise ValueError("Prior beta source mismatch")
    if annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or annual["width_m"]!=W or annual["qa"]["rows"]!=288721938 or any(annual[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Frozen annual data changed")
    cand,_=candidate_support(src)
    if annual["candidate_pairs"]!={g:len(cand[g]) for g in GROUPS}:raise ValueError("Annual candidate count changed")
    if any(not {pair_key(p) for p in annual["pairs"][y][g]}<=cand[g] for y in YEARS[:2] for g in GROUPS):raise ValueError("Annual geography changed")
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_gate_observations_observers_species"]!=[25,5,10] or src["split_seed"]!="tetsu-railstatus-observer-split-v1":raise ValueError("Wrong frozen split")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")) or any(src["qa"][g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Frozen source mismatch")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or am["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or am["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in s.CATS):raise ValueError("Wrong geometry")
    if tree["source_file_SHA256"]!=s.EXPECTED_SHA or road["source_zip_MD5"]!=s.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong road/tree")
    if s.sha_file(a.ne)!=s.ORIGINAL_NE_SHA or s.sha_file(a.climate)!=s.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=s.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong climate/map/pop")
    built,_=s.read_source_builtup(a.built);land,_=s.source_continents(a.ne);ras,mem,_=s.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=s.climate_scale(ras)
        w=W;grid=prior["grids"][str(w)]
        ps={y:{g:src["pairs"][str(w)][g] if y=="2025" else annual["pairs"][y][g] for g in GROUPS} for y in YEARS}
        keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for g in GROUPS for r in grid[g]}
        keys|={k for y in YEARS for g in GROUPS for p in ps[y][g] for k in pair_key(p)}
        clim,_=s.source_climate_lookup(land,ras,scale,keys,w)
        tr=s.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=s.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=s.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
        ag20=s.make_lookup(agb["grids"]["2020"][str(w)]);ag24=s.make_lookup(agb["grids"]["2024"][str(w)])
        el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(w)][c]} for c in s.CATS}
        raw={y:{g:{pair_key(p):p for p in ps[y][g]} for g in GROUPS} for y in YEARS};geo={}
        for y in YEARS:
            geo[y]={}
            for g in GROUPS:
                data=[s.convert(z) for z in grid[g]] if y=="2025" else cells_from_pairs(ps[y][g])
                rows,_,_=s.original_cell_rows(data,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                rows,_=s.inject_built(rows,built[str(w)]);rows,_=s.inject_status(rows,status,w);geo[y][g]=rows
        oldpanels={p["land"]:p for p in previous["panels"] if p["width_m"]==W}
        if set(oldpanels)!={.5,.8}:raise ValueError("Missing prior 200km land gates")
        for gate in s.LAND_THRESHOLDS:
            by={};counts={}
            for y in YEARS:
                groups={};ag={}
                for g in GROUPS:
                    rr,_=s.enrich_cells(geo[y][g],rd,tr,w,gate);rr,_=s.join_population(rr,pp,rd)
                    full=s.build_pairs([s.adapt_pair(p) for p in ps[y][g]],rr,el)
                    groups[g]=[z for r in full if (z:=s.split_pair(r,raw[y][g][r["a"],r["b"]])) is not None]
                    for k,v in b.attach_agb(rr,ag20,ag24,rd).items():
                        if k in ag and ag[k]!=v:raise ValueError("AGB grid mismatch")
                        ag[k]=v
                counts[y]={g:len(groups[g]) for g in GROUPS}
                if y=="2025" and any(counts[y][g]!=oldpanels[gate]["source_pair_counts"][g] for g in GROUPS):raise ValueError("Exact previous group support not replayed")
                by[y]=b.join(groups,raw[y],ag,GROUPS)
            validate_previous(by["2025"],oldpanels[gate]["plants_birds"])
            z={"land":gate,"source_counts":counts,**analyze(by)}
            panels.append(z);print(json.dumps({"plant_bird_beta_temporal_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-plant-bird-beta-temporal-v1","source_replay":checks,"biomass_source_SHA256":agb["source_sha256"],
         "source_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("annual","beta","split","edges","agb","agb_model")},"panels":panels,
         "primary":"200km land80, fixed common 2023-25 adjacent rail-positive pairs; opposite-observer-side plant and bird Jaccard similarity. 2023/24 equal-year block-FE slope, paired spatial-cluster covariance; flexible effort, endpoint FE, Simpson, 2020 biomass, land50 secondary.",
         "limits":["earlier period retrospective holdback for this composition endpoint, not independent ecosystems",
                   "2025 geography and same static 2026 railway map, biomass 2024 or 2020",
                   "within-taxon shared observer across endpoints adjusted, between taxa opposite observer halves",
                   "Jaccard includes nestedness and richness; not a pure turnover metric",
                   "conditional landscape similarity, no causal interaction, no actual rail movement",
                   "annual insect pair aggregates unavailable; no claim of temporal insect validation","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
