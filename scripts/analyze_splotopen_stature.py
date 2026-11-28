#!/usr/bin/env python3
"""Independent vegetation-plot landscape stature association, frozen transport map."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_disturbance_traits import GEO,CLIMATE,linear,combine,interval
from analyze_global_rail_height_road_contrast import road_lengths,ROAD_KM,ROAD_VALID,DENSITY
from scan_grip4_actual_road_lengths import candidates,combine as combine_vectors
from scan_splotopen_stature import PANELS,W,DOI,URL,LIFE_SHA,SIZE_SHA
from analyze_world_railstatus_coupling_climate import climate_extra
VECTOR="log_actual_grip4_road_km"
EFF=("log_plots","log_species","log_datasets","dominant_dataset_share","median_log_area","mean_year","unspecified_fraction","life_coverage","coverage")
def fit(rows,road_term=ROAD_KM,epoch=2020):
    rr=b.with_epoch(b.local_support(rows),epoch)
    fields=tuple(road_term if k==DENSITY else k for k in GEO)+CLIMATE+(b.BIO,)+EFF
    A,I,R=b.ACTIVE,b.INACTIVE,road_term
    if len(fields)!=len(set(fields)) or not {A,I,R}<=set(fields):raise ValueError("Duplicate/missing model terms")
    groups=defaultdict(list)
    for r in rr:groups[r["continent"],r["block"]].append(r)
    variation={k:sum(len({r[k] for r in v})>=2 for v in groups.values()) for k in (A,I,R)}
    qa={"n_cells":len(rr),"n_blocks":len(groups),"continents":dict(Counter(r["continent"] for r in rr)),"n_varying_blocks":variation}
    if any(n<20 for n in variation.values()):return dict(qa,status="insufficient exposure variation")
    fs={k:linear(rr,fields,"trait_value",k) for k in (A,I,R)}
    if any(v is None for v in fs.values()):return dict(qa,status="insufficient support or rank")
    e={"active":fs[A],"inactive":fs[I],"road":fs[R],"active_minus_inactive":combine(fs,{A:1.,I:-1.}),
       "active_minus_road":combine(fs,{A:1.,R:-1.}),"inactive_minus_road":combine(fs,{I:1.,R:-1.})}
    return dict(qa,status="estimated",effects={k:interval(v) for k,v in e.items()})
def join_source(rows,geo,rail):
    out=[]
    for s in rows:
        k=tuple(s["key"]);g=geo.get(k)
        if g is None:continue
        if s["n_plots"]<5 or s["known_height"]<20 or s["coverage"]<.5 or not s["known_height"]<=s["known_nonwoody"]<=s["species"] or not math.isclose(s["coverage"],s["known_height"]/s["known_nonwoody"],abs_tol=1e-12):raise ValueError("Independent plot gate changed")
        r=dict(s,key=k,continent=g["continent"],block=g["block"],bio1=g["bio1"],bio1_sq=g["bio1"]**2,abs_lat=g["abs_lat"],abs_lat_sq=g["abs_lat"]**2,
               rail_km=rail.get(k,0.),log_rail_km=math.log1p(rail.get(k,0.)),log_plots=math.log1p(s["n_plots"]),log_species=math.log1p(s["species"]),log_datasets=math.log1p(s["n_datasets"]))
        out.append(r)
    return out
def main():
    p=argparse.ArgumentParser()
    for n in ("source","annual_size","road_vectors","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();src,agb,plant,rail,status,tree,road,pop=map(b.read,(a.source,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    if src["version"]!="splotopen-independent-stature-grid-v1" or src["source"]["doi"]!=DOI or src["source"]["url"]!=URL or src["header_audit"]["qa"]["source_plots"]!=95104 or src["width"]!=W or src["traits"]["LifeForm_SHA256"]!=LIFE_SHA or src["traits"]["height_SHA256"]!=SIZE_SHA or set(src["grids"])!=set(PANELS):raise ValueError("Invalid independent source")
    checks=b.verify_original_sources(plant,rail)
    if status["version"]!="railstatus-v1-world-rail-km-by-status-and-class" or status["source_Overture_release"]!=rail["source_Oveture_release"] or agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1":raise ValueError("Wrong map/AGB")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen environmental sources differ")
    prior_keys=candidates(a.annual_size);vectors,receipts=combine_vectors(a.road_vectors,prior_keys)
    keys={tuple(r["key"]) for v in src["grids"].values() for r in v}
    result=[]
    if len(keys)<100:
        obj={"version":"splotopen-independent-stature-association-v1","source":src["source"],"source_support":src["panel_support"],"status":"insufficient trait-supported cells before environmental join","n_cells":len(keys),"panels":[]}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj));print(json.dumps(obj),flush=True);return
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    try:
        scale,_=b.climate_scale(ras);clim,_=b.source_climate_lookup(land,ras,scale,keys,W)
        if len(clim)<100:
            obj={"version":"splotopen-independent-stature-association-v1","source":src["source"],"source_support":src["panel_support"],"status":"insufficient climate land-center support","n_cells":len(clim),"panels":[]}
            a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj));print(json.dumps(obj),flush=True);return
        ext,_=climate_extra(a.climate,ras,clim,W)
    finally:ras.close();mem.close()
    rk={(z["source_equal_area_grid_ix"],z["source_equal_area_grid_iy"]):z["source_real_original_Overture_normal_rail_km_clipped_into_cell"] for z in rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)]}
    tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    ag20=b.make_lookup(agb["grids"]["2020"][str(W)]);ag24=b.make_lookup(agb["grids"]["2024"][str(W)])
    for label in PANELS:
        if len({tuple(r["key"]) for r in src["grids"][label]})!=len(src["grids"][label]):raise ValueError("Duplicate source grid")
        raw=join_source(src["grids"][label],clim,rk)
        raw,_=b.inject_built(raw,built[str(W)]);raw,_=b.inject_status(raw,status,W)
        for lg in (.8,.5):
            rr,_=b.enrich_cells(raw,rd,tr,W,lg);rr,_=b.join_population(rr,pp,rd);ready=[]
            for r in rr:
                k=r["key"]
                if k not in ext or k not in ag20 or k not in ag24:continue
                av,cv=ag20[k],ag24[k]
                if min(av["valid_area_m2"],cv["valid_area_m2"])<.9*rd[k]["source_GRIP4_total_land_area_m2"]:continue
                if any(not math.isfinite(v["agb_Mg_ha"]) or v["agb_Mg_ha"]<0 for v in (av,cv)):raise ValueError("Invalid biomass")
                z=road_lengths(dict(r,agb2020=av["agb_Mg_ha"],agb2024=cv["agb_Mg_ha"],**ext[k]),rd[k])
                ready.append(z)
            overlap=[dict(r,**{VECTOR:math.log1p(vectors[r["key"]])}) for r in ready if r["key"] in prior_keys]
            panel={"period":label,"land":lg,"n_source_cells":len(src["grids"][label]),"n_environment_ready":len(ready),
                   "primary":fit(ready),"AGB2024":fit(ready,epoch=2024),"valid_land_road_km":fit(ready,ROAD_VALID),
                   "original_candidate_overlap_raster":fit(overlap),"original_candidate_overlap_vector":fit(overlap,VECTOR),
                   "leave_continent_out":{c:fit([r for r in ready if r["continent"]!=c]) for c in sorted({r["continent"] for r in ready})}}
            result.append(panel);print(json.dumps({"splot_stature_panel":panel},allow_nan=False),flush=True)
    obj={"version":"splotopen-independent-stature-association-v1","source":src["source"],"source_support":src["panel_support"],"trait_audit":src["traits"],"map_replay":checks,"panels":result,
         "primary":"1990-2015 pooled plot species union, all global eligible 200km cells, land80 coverage50, minimum 5 plots and 20 known nonwoody species, active-minus-road and active-minus-inactive; 2020 AGB; mapped 2026 status",
         "limits":["Plot source independent of photo observations, shared trait sources and contemporary infrastructure geography; taxa name pool differs from iNat intersection",
                   "plots predate 2026 status and 2020 controls; association with present mapped geography only, no historical operation inference",
                   "sPlotOpen favors natural and semi-natural vegetation, not a trackside survey; selected plots, variable area and dataset coverage",
                   "species unions not local abundance; unknown traits not imputed; no causal management inference, no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
