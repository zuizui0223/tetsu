#!/usr/bin/env python3
"""Check whether spatial cross-taxon richness coupling survives additional
climate confounding controls. Same genuine source ZIP, frozen grid centers.
Do not alter the rail exposure, effort or support gates after seeing results.
"""
import argparse,json,math,hashlib,zipfile
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from rasterio.io import MemoryFile
from pyproj import Transformer
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_insect_bird import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,EXPECTED,extract,attach,linear,interval,combine,one,replay_check,
    CROSS,ADDED,extra_attach,fields_for,cv,model)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import cell_center,read_temp
BIO=(4,12,15)
CLIMATE=tuple(f"log_BIO{k}" for k in BIO)
EXTRA=CLIMATE+tuple(k+"_sq" for k in CLIMATE)

def climate_extra(path,reference,keys,width):
    """BIO4 temperature seasonality (source-scaled), BIO12 annual precipitation,
    BIO15 precipitation seasonality; log1p on native nonnegative source values.
    Same cell centers as existing BIO1; no extrapolation or zero imputation."""
    tf=Transformer.from_crs(6933,4326,always_xy=True)
    ks=sorted(keys);xy={k:tf.transform(*cell_center(*k,width)) for k in ks}
    values={k:{} for k in ks};audit={}
    with zipfile.ZipFile(path) as z:
        for n in BIO:
            candidates=[x for x in z.namelist() if x.endswith(f"bio_{n}.tif")]
            if len(candidates)!=1:raise ValueError("Nonunique climate raster")
            raw=z.read(candidates[0])
            with MemoryFile(raw) as mem,mem.open() as ras:
                if ras.count!=1 or ras.crs!=reference.crs or ras.shape!=reference.shape or ras.transform!=reference.transform:raise ValueError("BIO raster geometry differs")
                valid=[]
                for k,(lon,lat) in xy.items():
                    v=read_temp(ras,lat,lon)
                    if v is None:continue
                    if not math.isfinite(v) or not 0<=v<=1e6:raise ValueError("Nonphysical climate value")
                    f=f"log_BIO{n}";x=math.log1p(v);values[k][f]=x;values[k][f+"_sq"]=x*x;valid.append(v)
                if len(valid)<100:raise ValueError("Insufficient authentic climate coverage")
                audit[str(n)]={"member":candidates[0],"sha256":hashlib.sha256(raw).hexdigest(),"valid":len(valid),"range_native":[min(valid),max(valid)],"units":ras.units[0],"scale":ras.scales[0],"offset":ras.offsets[0]}
    return {k:v for k,v in values.items() if set(v)==set(EXTRA)},{"n_centers":len(ks),"n_complete":sum(set(v)==set(EXTRA) for v in values.values()),"raster_metadata":audit}

def attach_climate(rows,climate):
    return b.local_support([dict(r,**climate[r["key"]]) for r in rows if r["key"] in climate])

def fit_model(rows,groups,flex=True,extra=()):
    rr=b.with_epoch(b.local_support(rows),2024)
    fields=fields_for(groups,flex)+extra
    if len(set(fields))!=len(fields):raise ValueError("Duplicate climate terms")
    qa={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr})}
    fits={}
    for i,j in CROSS:
        plant=f"plants_{i}_species";target=f"birds_{j}_species"
        base=fields+(plant,);full=base+tuple(f"{g}_{i}_species" for g in groups)
        z={"plant_without_added_richness":linear(rr,base,target,plant),"plant_with_added_richness":linear(rr,full,target,plant)}
        z.update({g:linear(rr,full,target,f"{g}_{i}_species") for g in groups})
        if any(v is None for v in z.values()):return dict(qa,status="insufficient support or rank"),None
        z["plant_slope_attenuation"]=combine(z,{"plant_without_added_richness":1.,"plant_with_added_richness":-1.})
        if len(groups)==2:z["insects_minus_fungi"]=combine(z,{"insects":1.,"fungi":-1.})
        fits[i,j]=z
    scores={k:combine({ij:z[k] for ij,z in fits.items()},{ij:.5 for ij in CROSS}) for k in fits[CROSS[0]]}
    return dict(qa,status="estimated",effects={k:interval(v) for k,v in scores.items()}),scores

def check(new,old):
    if new["status"]!=old["status"] or new["n_cells"]!=old["n_cells"]:raise ValueError("Previous support/status mismatch")
    if new["status"]=="estimated":
        for k,v in new["effects"].items():
            q=old["effects"][k]
            if abs(v["estimate"]-q["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],q["interval95_approx"]))>1e-9:raise ValueError("Previous climate-baseline effect mismatch")

def evaluate(rows,groups,climate,old=None):
    rr=b.local_support(rows)
    if old is not None:
        for label,flex in (("primary",False),("flexible_effort",True)):
            z,_=fit_model(rr,groups,flex);check(z,old[label])
    matched=attach_climate(rr,climate)
    out={"n_prior_cells":len(rr),"n_matched_cells":len(matched),"models":{}}
    for label,flex in (("linear_effort",False),("flexible_effort",True)):
        a,sa=fit_model(matched,groups,flex);c,sc=fit_model(matched,groups,flex,EXTRA)
        z={"matched_baseline":a,"extra_climate":c}
        if sa is not None and sc is not None:
            z["baseline_minus_extended"]={k:interval(combine({"base":sa[k],"climate":sc[k]},{"base":1.,"climate":-1.})) for k in sa}
            z["holdout_extended"]=cv(b.with_epoch(matched,2024),fields_for(groups,flex)+EXTRA,groups)
        out["models"][label]=z
    out["leave_continent_out"]={c:fit_model([r for r in matched if r["continent"]!=c],groups,True,EXTRA)[0] for c in sorted({r["continent"] for r in matched})}
    return out

def main():
    p=argparse.ArgumentParser()
    for n in ("triad","agb","agb_model","disjoint","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();previous,agb,old,prior,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.triad,a.agb,a.agb_model,a.disjoint,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if previous["version"]!="railstatus-observer-disjoint-plant-insect-bird-biomass-v1" or previous["biomass_source_sha256"]!=agb["source_sha256"] or previous["prior_disjoint_SHA256"]!=hashlib.sha256(a.disjoint.read_bytes()).hexdigest():raise ValueError("Prior insect model mismatch")
    previous_panels={(p["width_m"],p["land"]):p for p in previous["panels"]}
    if prior["version"]!="railstatus-observer-disjoint-plant-bird-biomass-v1" or prior["biomass_source_sha256"]!=agb["source_sha256"] or any(prior["source_artifact_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("agb","agb_model","split","annual")):raise ValueError("Prior disjoint artifact mismatch")
    prior_panels={(p["width_m"],p["land"]):p for p in prior["panels"]}
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
        for w in b.GRID_WIDTHS:
            grid=taxa["grids"][str(w)];yr=("2025",)
            ps={y:{g:src["pairs"][str(w)][g] if y=="2025" else annual["pairs"][y][g] for g in GROUPS} for y in yr}
            maps={y:{g:extract(ps[y][g]) for g in GROUPS} for y in yr}
            prior_counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in GROUPS}
            for g in GROUPS:
                if any(prior_counts[g].get(k)!=v[0] for k,v in maps["2025"][g].items()):raise ValueError("2025 exact source replay failed")
            keys={b.key(z) for g in GROUPS for z in grid[g]}|{k for y in yr for g in GROUPS for k in maps[y][g]}
            clim,_=b.source_climate_lookup(land,raster,scale,keys,w)
            ext,qa=climate_extra(a.climate,raster,clim,w)
            print(json.dumps({"extra_climate_source":{"width_m":w,**qa}}),flush=True)
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
                oldp=prior_panels[w,gate]["2025"]
                replay_check(one(by["2025"])[0],oldp["primary"])
                replay_check(one(by["2025"],flex=True)[0],oldp["flexible_effort"])
                maps_extra={g:extract(src["pairs"][str(w)][g]) for g in ADDED}
                full_extra={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in ADDED}
                if any(full_extra[g].get(k)!=v[0] for g in ADDED for k,v in maps_extra[g].items()):raise ValueError("Added taxa full replay changed")
                r=by["2025"]
                prev=previous_panels[w,gate]
                z={"width_m":w,"land":gate,
                   "plant_bird":evaluate(r,(),ext),
                   "insect":evaluate(extra_attach(r,maps_extra,full_extra,("insects",)),("insects",),ext,prev["insect"]),
                   "both_same_cells":evaluate(extra_attach(r,maps_extra,full_extra,ADDED),ADDED,ext,prev["both_same_cells"])}
                # Exact two-taxon source validation on the original linear/flex disjoint panel
                for flex,label in ((False,"primary"),(True,"flexible_effort")):
                    q,_=fit_model(r,(),flex);prior_effect=oldp[label]["contrasts"]["cross"]
                    if q["n_cells"]!=oldp[label]["n_cells"] or q["status"]!="estimated" or abs(q["effects"]["plant_with_added_richness"]["estimate"]-prior_effect["estimate"])>1e-9:raise ValueError("Original plant-bird effect mismatch")
                panels.append(z);print(json.dumps({"climate_coupling_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    out={"version":"railstatus-observer-disjoint-richness-extra-climate-v1","source_replay":checks,
         "triad_source_sha256":hashlib.sha256(a.triad.read_bytes()).hexdigest(),"worldclim_sha256":b.sha_file(a.climate),"panels":panels,
         "primary":"100km land80 2025; compare plant and insect conditional slopes and plant attenuation on identical cells, flexible effort, additional WorldClim BIO4 BIO12 BIO15 log and squared terms; paired coefficient change and continent deletion",
         "limits":["grid-center climate normals, not whole-grid climate distribution or annual weather",
                   "no causal mediation, no rail-specific effect, no actual prey or bird biomass",
                   "climate controls may themselves affect biological richness; descriptive conditional association",
                   "same source platform and access, not independent field validation","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
