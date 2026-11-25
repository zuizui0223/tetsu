#!/usr/bin/env python3
"""Temporal replication of the 2025 bird endpoint-FE connectivity signal.
All-years common support fixed without selecting on observed coefficient.
"""
import argparse,json
from pathlib import Path
from analyze_world_railstatus_observer_split_connectivity import (
    adapt_pair,split_pair,both_fe,taxon_pair,PRED,ENDPOINT,ACTIVE,
    build_pairs,fit,POP_SHA,join_population,inject_status,CATS,
    read,make_lookup,enrich_cells,read_source_builtup,ORIGINAL_NE_SHA,WORLDCLIM_SHA,
    EXPECTED_SHA,MD5,sha_file,source_continents,worldclim_tif,climate_scale,
    source_climate_lookup,original_cell_rows,verify_original_sources,inject_built,
    verify_replay,continent_holdout)
from scan_world_railstatus_observer_annual_beta import candidate_support,pair_key,W,GROUPS
YEARS=("2023","2024","2025")
PRIOR_BIRD={.5:(1041,.02173511226854892),.8:(837,.018960805134287235)}
PRIOR_PLANT={.5:1077,.8:869}

def cells_from_pairs(pairs):
    counts={}
    for p in pairs:
        for side in ("a","b"):
            k=tuple(p[side]);v=tuple(p["full_"+side])
            if k in counts and counts[k]!=v:raise ValueError("Annual cell count inconsistent across pairs")
            counts[k]=v
    return [{"source_equal_area_grid_ix":k[0],"source_equal_area_grid_iy":k[1],
             "source_year_2025_photo_linked_plants_observation_count":v[0],
             "source_2025_distinct_observers_with_photo_linked_plant_records":v[1],
             "source_2025_plant_species_richness":v[2]} for k,v in sorted(counts.items())]

def subset(rows,keys):
    out=[r for r in rows if (r["a"],r["b"]) in keys]
    if len(out)!=len(keys):raise ValueError("Matched pair support changed")
    return out

def model(rows):
    z={}
    for target in ("cross_jaccard","cross_simpson","same_jaccard","same_minus_cross_jaccard","observer_jaccard"):
        z[target]=both_fe(rows,target)
    reduced=tuple(k for k in PRED if k!=ACTIVE)
    z["cv_jaccard"]=continent_holdout(rows,reduced,reduced+(ACTIVE,),"cross_jaccard",True,50,300,"annual disjoint Jaccard") if rows else None
    return z

def main():
    p=argparse.ArgumentParser()
    for n in ("annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    annual,src,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"]:raise ValueError("Wrong 2025 split source")
    if annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or annual["width_m"]!=W or annual["qa"]["rows"]!=288721938 or any(annual[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong annual source")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Prior SHA mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"]:raise ValueError("Invalid edge replay")
    if any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Status km mismatch")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    cand,cells=candidate_support(src)
    if annual["candidate_pairs"]!={g:len(cand[g]) for g in GROUPS}:raise ValueError("Frozen candidates changed")
    source={y:{g:annual["pairs"][y][g] for g in GROUPS} for y in YEARS[:-1]}
    source["2025"]={g:[r for r in src["pairs"][str(W)][g] if "split_similarity" in r] for g in GROUPS}
    for y in YEARS:
        for g in GROUPS:
            keys={pair_key(r) for r in source[y][g]}
            if len(keys)!=len(source[y][g]) or not keys<=cand[g] or any("split_similarity" not in r for r in source[y][g]):raise ValueError("Annual pair candidates changed")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);clim,_=source_climate_lookup(land,raster,scale,set.union(*cells.values()),W)
    finally:raster.close();mem.close()
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(W)][c]} for c in CATS}
    geo={}
    for y in YEARS:
        geo[y]={}
        for g in GROUPS:
            raw=cells_from_pairs(source[y][g])
            if y=="2025":
                pmap={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]):(r["observations"],r["observers"],r["species"]) for r in prior["grids"][str(W)][g]}
                for r in raw:
                    key=(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"])
                    v=(r["source_year_2025_photo_linked_plants_observation_count"],r["source_2025_distinct_observers_with_photo_linked_plant_records"],r["source_2025_plant_species_richness"])
                    if pmap.get(key)!=v:raise ValueError("2025 exact cell replay failed")
            rows,_,_=original_cell_rows(raw,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
            rows,_=inject_built(rows,built[str(W)]);rows,_=inject_status(rows,status,W);geo[y][g]=rows
    panels=[]
    for threshold in (.5,.8):
        rr={}
        for y in YEARS:
            rr[y]={}
            for g in GROUPS:
                m,_=enrich_cells(geo[y][g],rd,tr,W,threshold);m,_=join_population(m,pp,rd)
                full=build_pairs([adapt_pair(p) for p in source[y][g]],m,el)
                lookup={pair_key(p):p for p in source[y][g]}
                rr[y][g]=[split_pair(r,lookup[r["a"],r["b"]]) for r in full]
        old=fit(rr["2025"]["birds"],"cross_jaccard",True,ENDPOINT)
        n,b=PRIOR_BIRD[threshold]
        if len(rr["2025"]["birds"])!=n or len(rr["2025"]["plants"])!=PRIOR_PLANT[threshold] or old["status"]!="estimated" or abs(old["effects"][ACTIVE]["estimate"]-b)>1e-9:raise ValueError("2025 source effect not exactly reproduced")
        keys={g:set.intersection(*({(r["a"],r["b"]) for r in rr[y][g]} for y in YEARS)) for g in GROUPS}
        common=keys["plants"]&keys["birds"]
        z={"land":threshold,"n_all_year_common":{g:len(keys[g]) for g in GROUPS},"n_both_taxa_all_year_common":len(common),
           "annual_native_fixed_2025_candidate":{},"annual_same_pairs":{},"plant_bird_same_pairs":{}}
        for y in YEARS:
            z["annual_native_fixed_2025_candidate"][y]={g:model(rr[y][g]) for g in GROUPS}
            z["annual_same_pairs"][y]={g:model(subset(rr[y][g],keys[g])) for g in GROUPS}
            z["plant_bird_same_pairs"][y]=taxon_pair(subset(rr[y]["plants"],common),subset(rr[y]["birds"],common))
        panels.append(z);print(json.dumps({"annual_observer_beta_panel":z},allow_nan=False),flush=True)
    out={"version":"railstatus-observer-disjoint-bird-temporal-v1","replay":checks,"source_qa":annual["qa"],"panels":panels,
         "primary":"200km land>=50%, bird observer-disjoint Jaccard endpoint FE active crossing on exact all-three-year common pairs; separate 2023 and 2024 retests, 2025 benchmark; plant minus bird on exact both-taxon all-year support, 80% land/Simpson secondary",
         "limits":["2025 geography selection, earlier years are retrospective holdouts, not independent ecosystems","same observers may occur across years but opposite halves are disjoint within each year","no historical train operation, 2026 geometry fixed","whole-landscape photo-recorded species similarity","approximate block intervals and correlated years; no pooled independent p-values"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
