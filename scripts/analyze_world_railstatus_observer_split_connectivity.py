#!/usr/bin/env python3
"""Paired observer-disjoint beta test, exact common pair support.
Split half outcomes share one symmetric effort/richness adjustment design.
The same-minus-cross outcome gives a paired coefficient difference and
block covariance without comparing separate significance tests.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
from analyze_world_railstatus_direct_boundary_beta import (
    GRID,PAIRS,ACTIVE,INACTIVE,CROSS,CELL_FIELDS,MEAN,DIFF,PRED,build_pairs,fit,
    POP_SHA,join_population,inject_status,CATS,GRID_WIDTHS,LAND_THRESHOLDS,read,
    make_lookup,enrich_cells,read_source_builtup,ORIGINAL_NE_SHA,WORLDCLIM_SHA,
    EXPECTED_SHA,MD5,sha_file,source_continents,worldclim_tif,climate_scale,
    source_climate_lookup,original_cell_rows,verify_original_sources,inject_built,
    continent_holdout)
from analyze_world_railstatus_cross_taxa import GROUPS,convert,verify_replay

EFF=("log_observations","log_observers","log_species")
EFF_PAIR=tuple(prefix+k for prefix in ("mean_","diff_") for k in EFF)
FIXED_PAIR=tuple(k for k in PRED if k not in EFF_PAIR)
ENDPOINT=DIFF+("horizontal",)+CROSS

def adapt_pair(p):
    # Source counts in observation, observer, species order.
    a,b=p["full_a"],p["full_b"];s=p["species_full"]
    if s["union"]!=a[2]+b[2]-s["shared"]:raise ValueError("Full pair union mismatch")
    return {"source_cell_A_original_grid_ix_iy":p["a"],
            "source_cell_B_original_grid_ix_iy":p["b"],
            "source_grid_A_original_plant_species_n":a[2],
            "source_grid_B_original_plant_species_n":b[2],
            "source_shared_plants_species_n":s["shared"],
            "source_union_plants_species_n":s["union"],
            "source_biological_Jaccard_similarity_between_ENTIRE_landscape_cells":s["jaccard"]}

def symmetric_effort(a,b):
    # Common design for both same/cross outcomes. Mean log effort/richness
    # across four endpoint halves; absolute differences averaged over all
    # four possible half-pair combinations. No outcome-specific weighting.
    out={}
    for i,k in enumerate(EFF):
        x=[math.log1p(v[i]) for v in a];y=[math.log1p(v[i]) for v in b]
        out["mean_"+k]=math.fsum(x+y)/4
        out["diff_"+k]=math.fsum(abs(v-u) for v in x for u in y)/4
    return out

def split_pair(base,p):
    if "split_similarity" not in p:return None
    a,b=p["split_counts_a"],p["split_counts_b"]
    if len(a)!=2 or len(b)!=2 or any(not all(v>=g for v,g in zip(x,(25,5,10))) for x in a+b):
        raise ValueError("Split support mismatch")
    r=dict(base,**symmetric_effort(a,b))
    for name,key in (("jaccard","jaccard"),("simpson","simpson_similarity")):
        v={k:p["split_similarity"][k][key] for k in ("00","01","10","11")}
        if not all(x is not None and 0<=x<=1 for x in v.values()):raise ValueError("Invalid split similarity")
        r["same_"+name]=(v["00"]+v["11"])/2
        r["cross_"+name]=(v["01"]+v["10"])/2
        r["same_minus_cross_"+name]=r["same_"+name]-r["cross_"+name]
        for k in v:r[name+"_"+k]=v[k]
    ov=p["observer_full"]["jaccard"]
    if ov is None or not 0<=ov<=1:raise ValueError("Invalid observer overlap")
    r["observer_jaccard"]=ov
    r["same_observer_jaccard"]=sum(p["same_side_observer_overlap"][str(s)]["jaccard"] for s in (0,1))/2
    return r

def both_fe(rows,target,fields=PRED,endpoint_fields=ENDPOINT):
    return {"block_FE":fit(rows,target,False,fields),"endpoint_FE":fit(rows,target,True,endpoint_fields)}

def compare(full,split):
    keys={(r["a"],r["b"]) for r in split}
    if len(keys)!=len(split):raise ValueError("Duplicate eligible pairs")
    common=[r for r in full if (r["a"],r["b"]) in keys]
    if len(common)!=len(split):raise ValueError("Full/split pair support drift")
    out={"n_full_pairs":len(full),"n_common_split_pairs":len(split),"full_native":{},
         "full_common":{},"split_common":{},"observer_diagnostic":{}}
    for name in ("jaccard","simpson"):
        out["full_native"][name]=both_fe(full,name)
        out["full_common"][name]=both_fe(common,name)
        out["split_common"][name]={}
        for target in ("cross_"+name,"same_"+name,"same_minus_cross_"+name):
            out["split_common"][name][target]=both_fe(split,target)
        # Matched design implies paired differences of estimates are exact.
        z=out["split_common"][name]
        for fe in ("block_FE","endpoint_FE"):
            a,b,d=[z[t][fe] for t in ("same_"+name,"cross_"+name,"same_minus_cross_"+name)]
            if all(v["status"]=="estimated" for v in (a,b,d)):
                for f in CROSS:
                    aa,bb,dd=[v["effects"][f] for v in (a,b,d)]
                    if aa is not None and bb is not None and dd is not None and abs(aa["estimate"]-bb["estimate"]-dd["estimate"])>1e-8:
                        raise ValueError("Paired same-minus-cross identity failed")
        if split:
            reduced=tuple(k for k in PRED if k!=ACTIVE)
            out["split_common"][name]["cross_continent_holdout"]=continent_holdout(split,reduced,reduced+(ACTIVE,),"cross_"+name,True,50,300,"disjoint "+name)
    out["observer_diagnostic"]["observer_overlap_on_common"]=both_fe(split,"observer_jaccard")
    out["observer_diagnostic"]["same_side_overlap_on_common"]=both_fe(split,"same_observer_jaccard")
    out["observer_diagnostic"]["disjoint_jaccard_adjusted_overlap"]=both_fe(split,"cross_jaccard",PRED+("observer_jaccard",),ENDPOINT+("observer_jaccard",))
    return out

def taxon_pair(a,b):
    # Same geography and identical joint design for plant, comparator and
    # difference. Each taxon's effort covariates enter separately.
    A={(r["a"],r["b"]):r for r in a};B={(r["a"],r["b"]):r for r in b}
    keys=sorted(set(A)&set(B));rr=[]
    for k in keys:
        p,q=A[k],B[k]
        if any(p[f]!=q[f] for f in FIXED_PAIR):raise ValueError("Common geographical design mismatch")
        r=dict(p,**{"other_"+f:q[f] for f in EFF_PAIR})
        r["other_cross_jaccard"]=q["cross_jaccard"]
        r["plant_minus_other_cross_jaccard"]=p["cross_jaccard"]-q["cross_jaccard"]
        rr.append(r)
    fields=PRED+tuple("other_"+f for f in EFF_PAIR)
    endpoint=ENDPOINT+tuple("other_"+f for f in EFF_PAIR if f.startswith("diff_"))
    return {"n_pairs":len(rr),"effects":{t:both_fe(rr,t,fields,endpoint) for t in
          ("cross_jaccard","other_cross_jaccard","plant_minus_other_cross_jaccard")}}

def main():
    p=argparse.ArgumentParser()
    for n in ("split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    source,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if source["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not source["exact_grid_replay"] or not source["exact_plant_beta_replay"]:raise ValueError("Unverified split source")
    if any(source[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")) or any(source["qa"][g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Split source hash/count mismatch")
    if source["split_gate_observations_observers_species"]!=[25,5,10]:raise ValueError("Split gate drift")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"]:raise ValueError("Invalid edge replay")
    if any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Status km mismatch")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=climate_scale(raster)
        for w in GRID_WIDTHS:
            grid=prior["grids"][str(w)];sp=source["pairs"][str(w)]
            keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for g in GROUPS for r in grid[g]}
            clim,_=source_climate_lookup(land,raster,scale,keys,w)
            tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
            rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
            pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(w)][c]} for c in CATS}
            geo={}
            for g in GROUPS:
                rows,_,_=original_cell_rows([convert(z) for z in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                rows,_=inject_built(rows,built[str(w)]);rows,_=inject_status(rows,status,w);geo[g]=rows
            for threshold in LAND_THRESHOLDS:
                groups={};split_rows={}
                for g in GROUPS:
                    m,_=enrich_cells(geo[g],rd,tr,w,threshold);m,_=join_population(m,pp,rd)
                    full=build_pairs([adapt_pair(p) for p in sp[g]],m,el)
                    lookup={(tuple(p["a"]),tuple(p["b"])):p for p in sp[g]}
                    split=[z for r in full if (z:=split_pair(r,lookup[r["a"],r["b"]])) is not None]
                    split_rows[g]=split;groups[g]=compare(full,split)
                panel={"width_m":w,"min_land_fraction":threshold,"groups":groups,
                       "plant_vs_other":{g:taxon_pair(split_rows["plants"],split_rows[g]) for g in GROUPS if g!="plants"}}
                panels.append(panel)
                print(json.dumps({"observer_split_panel":panel},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-observer-split-connectivity-v1","source_replay":checks,"panels":panels,
        "primary":"100km, 50% land, plants disjoint observer Jaccard effect conditional on symmetric half-specific effort, same support and block FE; endpoint FE, same-half paired difference, 200km and other taxa sensitivity",
        "limits":["posthoc observational falsification not independent replication","disjoint observer IDs do not eliminate shared sampling preferences or site accessibility","half sample gates change support, compare full-common with full-native","simultaneous effort/richness controls are descriptive","whole landscape communities, physical crossing not train service","all intervals approximate; multiple correlated comparisons"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
