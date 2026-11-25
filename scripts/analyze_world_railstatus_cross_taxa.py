#!/usr/bin/env python3
"""Matched support taxon specificity of rail-status predictive signal.
Same v369 source, v382 covariates, railstatus categories and frozen ridge.
Compare relative held-out MSE gain, signed partial slopes, paired block-cluster
intervals; observational associations only.
"""
import argparse, hashlib, json, math, random
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from analyze_world_rail_status_richness_railstatus import (
    CATS, LOG, ROAD, BASE, LAND_THRESHOLDS, GRID_WIDTHS, CONTINENTS,
    MIN_CONTINENT_RICH_CELLS, N_PLACEBO, band, ho, inject_status,
    read, make_lookup, enrich_cells, read_source_builtup,
    ORIGINAL_NE_SHA, WORLDCLIM_SHA, EXPECTED_SHA, MD5, sha_file,
    source_continents, worldclim_tif, climate_scale, source_climate_lookup,
    original_cell_rows, verify_original_sources)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights

GROUPS=("plants","birds","insects","fungi")
BASELINE=BASE+ROAD
FULL=BASELINE+tuple(LOG[c] for c in CATS)
ACTIVE=LOG["active_normal"]
INACTIVE=LOG["inactive_normal"]
FIELDS=("source_equal_area_grid_ix","source_equal_area_grid_iy")
OLD_GRID="source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"

def key(z): return tuple(z[k] for k in FIELDS)

def convert(z):
    return dict(z,source_year_2025_photo_linked_plants_observation_count=z["observations"],
                source_2025_plant_species_richness=z["species"],
                source_2025_distinct_observers_with_photo_linked_plant_records=z["observers"])

def verify_replay(source,old):
    s=old["source_summary"]
    for new,prior in (("obs_etag","original_source_iNat_obs_S3_ETag"),
                      ("taxa_etag","original_source_iNat_taxa_S3_ETag"),
                      ("obs_sha256","source_iNat_original_2026_10_27_observations_gzip_SHA256_full"),
                      ("taxa_sha256","source_iNat_original_2026_10_27_taxa_gzip_SHA256_full")):
        if not source[new] or source[new]!=s[prior]: raise ValueError("Frozen source replay mismatch "+new)
    for width in GRID_WIDTHS:
        a={key(z):(z["observations"],z["species"],z["observers"]) for z in source["grids"][str(width)]["plants"]}
        b={key(z):(z["source_year_2025_photo_linked_plants_observation_count"],
                   z["source_2025_plant_species_richness"],
                   z["source_2025_distinct_observers_with_photo_linked_plant_records"])
           for z in old["source_grids"][str(width)][OLD_GRID]}
        if a!=b: raise ValueError("Plant counts not exactly identical to v369 at "+str(width))

def restrict(rows,keys):
    out=sorted((r for r in rows if r["key"] in keys),key=lambda r:r["key"])
    if len(out)!=len(keys): raise ValueError("Matched support drift")
    return out

def rel(r):
    v=r["mean_dMSE"]; b=r["base_mean_MSE"]
    return float(v/b) if v is not None and b>0 else None

def signed_cluster(rows,field):
    """FWL within-continent coefficient and paired block-cluster influence.
    Same equal-block weights as v382; interval is conditional on these cells,
    not a causal or independent-continent confidence interval.
    """
    if len(rows)<100: return None,{}
    w=block_weights(rows)
    X=np.array([[r[k] for k in FULL]+[r["log_species"]] for r in rows],float)
    for c in CONTINENTS:
        ix=np.array([i for i,r in enumerate(rows) if r["continent"]==c],int)
        if len(ix): X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None]
    k=FULL.index(field);z=X[:,k];y=X[:,-1]
    other=np.delete(X[:,:-1],k,axis=1)
    rz=z-other@np.linalg.lstsq(other,z,rcond=1e-10)[0]
    ry=y-other@np.linalg.lstsq(other,y,rcond=1e-10)[0]
    den=float(rz@rz)
    if den<=1e-12: return None,{}
    beta=float(rz@ry/den); influence=defaultdict(float)
    for r,v in zip(rows,rz*(ry-beta*rz)/den):
        influence[(r["continent"],r["block"])]+=float(v)
    return beta,dict(influence)

def paired_interval(b1,i1,b2,i2):
    if b1 is None or b2 is None or set(i1)!=set(i2) or len(i1)<2: return None
    n=len(i1); delta=b1-b2
    se=math.sqrt(n/(n-1)*sum((i1[k]-i2[k])**2 for k in i1))
    return {"plant_minus_other":delta,"block_cluster_SE":se,"interval95_approx":[delta-1.96*se,delta+1.96*se],"n_blocks":n}

def contrast(rows,field,conditional=True):
    base=tuple(k for k in FULL if k!=field) if conditional else BASELINE
    result=ho(rows,base,base+(field,),field)
    result["relative_MSE_gain"]=rel(result)
    return result

def permutation_indices(rows,rng):
    groups=defaultdict(list)
    for i,r in enumerate(rows): groups[band(r)].append(i)
    perm=list(range(len(rows)))
    for idx in groups.values():
        src=idx[:];rng.shuffle(src)
        for i,j in zip(idx,src): perm[i]=j
    return perm

def summarize(rows_by_group,label,placebos,seed):
    groups=list(rows_by_group)
    first=rows_by_group[groups[0]]
    if not first: return {"label":label,"n_cells":0,"status":"no common support"}
    for g in groups:
        if [r["key"] for r in rows_by_group[g]]!=[r["key"] for r in first]:
            raise ValueError("Matched cells not identical")
    results={}; signed={}
    for g,rows in rows_by_group.items():
        results[g]={"road_plus":{c:contrast(rows,LOG[c],False) for c in CATS},
                    "conditional":{c:contrast(rows,LOG[c]) for c in CATS},
                    "n_cells_by_continent":dict(Counter(r["continent"] for r in rows))}
        signed[g]={c:signed_cluster(rows,LOG[c]) for c in CATS}
        results[g]["within_continent_full_slopes"]={c:b[0] for c,b in signed[g].items()}
    diffs={}
    if "plants" in groups:
        for g in groups:
            if g=="plants": continue
            d=paired_interval(*signed["plants"]["active_normal"],*signed[g]["active_normal"])
            pg=results["plants"]["conditional"]["active_normal"]["relative_MSE_gain"]
            gg=results[g]["conditional"]["active_normal"]["relative_MSE_gain"]
            diffs[g]={"active_slope_paired_difference":d,
                      "active_relative_gain_difference":pg-gg if pg is not None and gg is not None else None}
    if placebos:
        rng=random.Random(hashlib.sha256(seed.encode()).hexdigest())
        for c in ("active_normal","inactive_normal"):
            field=LOG[c]
            null={g:[] for g in groups}
            for _ in range(N_PLACEBO):
                perm=permutation_indices(first,rng)
                for g,rows in rows_by_group.items():
                    sh=[dict(r,**{field:rows[j][field]}) for r,j in zip(rows,perm)]
                    null[g].append(contrast(sh,field)["relative_MSE_gain"])
            for g in groups:
                actual=results[g]["conditional"][c]["relative_MSE_gain"]
                vals=[v for v in null[g] if v is not None]
                results[g].setdefault("conditional_band_shuffle",{})[c]={
                    "n_null":len(vals),"n_null_ge_real":sum(v>=actual for v in vals) if actual is not None else None,
                    "randomization_p":(1+sum(v>=actual for v in vals))/(1+len(vals)) if actual is not None and len(vals)==N_PLACEBO else None,
                    "null_relative_gains":vals}
            if c=="active_normal" and "plants" in groups:
                for g in diffs:
                    actual=diffs[g]["active_relative_gain_difference"]
                    vals=[a-b for a,b in zip(null["plants"],null[g]) if a is not None and b is not None]
                    diffs[g]["paired_shared_map_null"]={"n_null":len(vals),"n_null_ge_real":sum(v>=actual for v in vals) if actual is not None else None,
                        "note":"shared-map reference for relative-gain contrast; not an exact test of equal nonzero effects"}
    return {"label":label,"n_cells":len(first),"n_blocks":len({(r["continent"],r["block"]) for r in first}),
            "results":results,"plant_vs_other":diffs}

def main():
    p=argparse.ArgumentParser()
    for n in ("taxa","plant","rail","status","built","tree","road","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    source,plant,rail,status,tree,road=map(read,(a.taxa,a.plant,a.rail,a.status,a.tree,a.road))
    verify_replay(source,plant); checks=verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5: raise ValueError("Frozen land source mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA: raise ValueError("Frozen climate/map mismatch")
    if not status["version"].startswith("railstatus-v1"): raise ValueError("Rail status source mismatch")
    from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import inject_built
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne)
    raster,mem,_=worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=climate_scale(raster)
        for width in GRID_WIDTHS:
            grid=source["grids"][str(width)]
            keys={key(z) for g in GROUPS for z in grid[g]}
            clim,_=source_climate_lookup(land,raster,scale,keys,width)
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
            geo={}
            for g in GROUPS:
                rows,_,_=original_cell_rows([convert(z) for z in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)],clim)
                rows,_=inject_built(rows,built[str(width)])
                rows,_=inject_status(rows,status,width)
                geo[g]=rows
            for threshold in LAND_THRESHOLDS:
                rows={g:enrich_cells(geo[g],rd,tr,width,threshold)[0] for g in GROUPS}
                support={g:{r["key"] for r in rows[g]} for g in GROUPS}
                seed=f"railstatus-cross-taxa-{width}-{threshold}"
                # Primary pairwise matched support; all-group and native panels are sensitivity checks.
                fits=[]
                for g in GROUPS[1:]:
                    shared=support["plants"]&support[g]
                    fits.append(summarize({x:restrict(rows[x],shared) for x in ("plants",g)},g+"_matched",True,seed+g))
                common=set.intersection(*support.values())
                fits.append(summarize({g:restrict(rows[g],common) for g in GROUPS},"all_four_matched",False,seed+"all"))
                for g in GROUPS:
                    fits.append(summarize({g:restrict(rows[g],support[g])},g+"_native",False,seed+g+"native"))
                panel={"width_m":width,"min_land_fraction":threshold,"fits":fits}
                panels.append(panel)
                print(json.dumps({"width_m":width,"min_land_fraction":threshold,
                   "fits":[{"label":f["label"],"n_cells":f["n_cells"],
                            "active":{g:{"conditional":r["conditional"]["active_normal"],"slope":r["within_continent_full_slopes"]["active_normal"],"placebo":r.get("conditional_band_shuffle",{}).get("active_normal",{}).get("randomization_p")} for g,r in f.get("results",{}).items()},
                            "plant_vs_other":f.get("plant_vs_other",{})} for f in fits]},allow_nan=False),flush=True)
    finally: raster.close();mem.close()
    result={"version":"railstatus-cross-taxa-v1-analysis","replay":checks,"source_qa":source["qa"],"panels":panels,
            "interpretation":"Plant specificity is suggestive, not a causal discrimination: taxon-specific observer behavior, detectability, and actual rail responses can all differ. Common effects do not establish pure observer bias.",
            "primary":"pairwise exact common cells, 100/200 km, 50% land; 80% land and all-four common cells sensitivity; full status conditional and road-plus; 39 shared band shuffles",
            "limits":["outcomes are recorded richness across whole landscapes","species threshold is the same v369 gate, restricting low-diversity cells","cluster intervals are approximate, conditional on sampled blocks","continent folds share training data","status tagging uneven; no train-frequency measure","relative gain differences not direct equal-effect p-values"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False),encoding="utf8")
if __name__=="__main__":main()
