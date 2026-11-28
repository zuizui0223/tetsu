#!/usr/bin/env python3
"""One distinct data-source replication on existing railway-positive geography."""
import argparse,json
from pathlib import Path
from analyze_global_hwi_road_crossings import (road_panel,fit,ACTIVE,FIELDS,END,INACTIVE,CATS,W,read,verify_original_sources,verify_replay,EXPECTED_SHA,MD5,POP_SHA,ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,read_source_builtup,source_continents,worldclim_tif,climate_scale,source_climate_lookup,make_lookup,original_cell_rows,inject_built,inject_status,enrich_cells,join_population,cells_from_pairs,adapt_pair,build_pairs,split_pair,pair_key,SOURCE_SHA,MIN_GUILD_HALF)
from scan_global_grip4_road_crossings import candidates,key_sha,combine
from scan_xenocanto_hwi_replication import PANELS,DOI,MAPPING
def main():
    p=argparse.ArgumentParser()
    for n in ("xc","historical","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","road-crossings","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();xc,hist,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.xc,a.historical,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    common=candidates(a.historical)
    if xc["version"]!="xc-independent-hwi-pairs-v2" or xc["doi"]!=DOI or xc["trait_mapping"]["method"]!=MAPPING or xc["trait_mapping"]["HWI_threshold"]!=21.0 or xc["trait_mapping"]["identical_to_original_iNat_matched_pool"] is not False or xc["trait_sha256"]!=SOURCE_SHA or xc["candidate_sha256"]!=key_sha(common) or xc["n_candidates"]!=len(common) or xc["guild_half_gate"]!=list(MIN_GUILD_HALF):raise ValueError("Wrong source/trait/geography")
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if any(hist[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Prior source mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Invalid geometry")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);ras,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(ras);clim,_=source_climate_lookup(land,ras,scale,{k for pair in common for k in pair},W)
    finally:ras.close();mem.close()
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)]);rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)]);pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(W)][c]} for c in CATS}
    by,receipts=combine(a.road_crossings,common)
    panels=[]
    for period in PANELS:
        gg={pair_key(p):p for p in xc["pairs"][period]}
        if len(gg)!=len(xc["pairs"][period]) or not set(gg)<=common:raise ValueError("Candidate support mismatch")
        raw=[p["full"] for k,p in sorted(gg.items())]
        for landgate in (.5,.8):
            if not raw:
                result={"period":period,"land":landgate,"status":"no full observer-disjoint pairs"}
            else:
                rows,_,_=original_cell_rows(cells_from_pairs(raw),rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
                rows,_=inject_built(rows,built[str(W)]);rows,_=inject_status(rows,status,W);rows,_=enrich_cells(rows,rd,tr,W,landgate);rows,_=join_population(rows,pp,rd)
                pairs=build_pairs([adapt_pair(p) for p in raw],rows,el)
                rr=[split_pair(r,gg[r["a"],r["b"]]["full"]) for r in pairs]
                result={"period":period,"land":landgate,"n_full_pairs_after_geographic_controls":len(rr),"results":road_panel(rr,gg,by)}
            panels.append(result);print(json.dumps({"xc_replication_panel":result},allow_nan=False),flush=True)
    obj={"version":"xc-independent-hwi-road-replication-v2","source_sha256":xc["source_sha256"],"source_audit":xc["audit"],"trait_mapping":xc["trait_mapping"],"prior_source_replay":checks,"panels":panels,"primary":"2018-25 pooled broad uncertainty, land50 habitat x migration standardized, low-minus-high HWI Jaccard, block FE, active-minus-highway direct contrast; active and active-minus-inactive companions","limits":["Different recording platform, same geography and frozen trait table/threshold, expanded AVONET source pool rather than exact prior iNat name intersection; recordists may overlap across platforms","unknown uncertainty is not verified accuracy","recordist names are not guaranteed unique human identities","incidence similarity is not movement or occupancy","retrospective static 2026 map, rail-positive endpoints only","no data from eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False))
if __name__=="__main__":main()
