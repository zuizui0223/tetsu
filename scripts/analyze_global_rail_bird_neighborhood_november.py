"""New source version, unchanged frozen ecological design."""
import argparse, collections, hashlib, json, math
from pathlib import Path
from pyproj import Transformer
from analyze_global_rail_bird_neighborhood_source import (
    s,p,W,as_counts,analyze,selftest)
from audit_railstatus_inat_november_snapshot import Snapshot,taxonomy,traits
def scan(snapshot,canonical):
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    z={k:collections.defaultdict(p.bucket) for k in ("early","late","2025")}
    qa=collections.Counter()
    rr,raw,stream=snapshot.open("observations.csv.gz")
    cols=set(rr.fieldnames or [])
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=cols:
        raise ValueError("observation schema changed")
    if {"geoprivacy","coordinates_obscured","taxon_geoprivacy"}&cols:
        raise ValueError("observation privacy schema changed")
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%p.PROGRESS_EVERY_SOURCE_ROWS==0:
                print(json.dumps({"global_neighborhood_scan":dict(qa)}),flush=True)
            yr=p.period(r["observed_on"])
            if yr is None or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:
                hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"])
                lon,lat=float(r["longitude"]),float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0]!="birds" or not math.isfinite(u) or not 0<=u<=p.MAX_UNCERT_M:continue
            k1=p.source_coarse_cell(lon,lat,tf,100000)
            k=p.source_coarse_cell(lon,lat,tf,W)
            if k1 is None or k is None:continue
            o=int(r["observer_id"]) if r["observer_id"] else None
            p.add(z[yr][k],hit[1],o);qa[yr+"_bird_records"]+=1
            if r["observed_on"].startswith("2025-"):
                p.add(z["2025"][k],hit[1],o);qa["2025_bird_records"]+=1
    finally:stream.close()
    if qa["rows"]<288721938 or raw.bytes!=snapshot.files["observations.csv.gz"]["bytes"] or qa["2025_bird_records"]<3000000:
        raise ValueError("incomplete source")
    snapshot.verify("observations.csv.gz")
    snapshot.verify("taxa.csv.gz")
    return z,dict(qa),raw.hash.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    for k in ("conditional","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        ap.add_argument("--"+k,type=Path)
    a=ap.parse_args();selftest()
    if a.self_test:return
    cond,annual,src,prior,edge,plant,rail,status,tree,road,pop=map(s.read,(a.conditional,a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=s.verify_original_sources(plant,rail);s.verify_replay(prior,plant)
    if cond["version"]!="railstatus-avonet-conditional-hwi-pairs-v1" or not cond["exact_full_pair_replay"] or not cond["exact_hwi_pair_replay"] or cond["source"]["source_SHA256"]!=p.SOURCE_SHA:
        raise ValueError("unverified trait source")
    if any(len({cond[k],annual[k],src[k]})!=1 for k in ("obs_sha256","taxa_sha256","split_seed")) or src["split_seed"]!=p.SALT.decode() or any(prior[k]!=src[k] for k in ("obs_sha256","taxa_sha256")):
        raise ValueError("prior source mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in s.CATS):
        raise ValueError("map geometry mismatch")
    if tree["source_file_SHA256"]!=s.EXPECTED_SHA or road["source_zip_MD5"]!=s.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":
        raise ValueError("environment source mismatch")
    if s.sha_file(a.ne)!=s.ORIGINAL_NE_SHA or s.sha_file(a.climate)!=s.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=s.POP_SHA or pop["relative_population_sum_difference"]>.01:
        raise ValueError("environment source mismatch")
    raw=p.request(p.DERIVATIVE,cap=10_000_000)
    if hashlib.sha256(raw).hexdigest()!=p.SOURCE_SHA:raise ValueError("trait source changed")
    snapshot=Snapshot(p.source_client())
    print(json.dumps({"monthly_source_lock":snapshot.files,"old_archive":snapshot.old_archive}),flush=True)
    canonical,names,txqa,sha=taxonomy(snapshot)
    traitmap,traitqa=traits(p.validate(p.csv_rows(raw)),names)
    print(json.dumps({"new_taxonomy_qa":txqa,"new_trait_qa":traitqa,"taxa_sha256":sha}),flush=True)
    z,scanqa,osha=scan(snapshot,canonical)
    old={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]):(r["observations"],r["observers"],r["species"]) for r in prior["grids"][str(W)]["birds"]}
    new={k:p.full_counts(v) for k,v in z["2025"].items()}
    if sum(v[0] for v in new.values())!=scanqa["2025_bird_records"]:
        raise ValueError("new 2025 internal census inconsistent")
    common=set(new)&set(old)
    drift={"old_cells":len(old),"new_cells":len(new),"common_cells":len(common),
           "identical_common_counts":sum(old[k]==new[k] for k in common),
           "old_records":sum(v[0] for v in old.values()),"new_records":scanqa["2025_bird_records"],
           "old_snapshot_sha256":src["obs_sha256"],"new_snapshot_sha256":osha}
    print(json.dumps({"monthly_census_drift":drift,"eligible_baseline_cells":sum(p.supported(v) for v in z["early"].values()),"scan_qa":scanqa}),flush=True)
    del z["2025"]
    rawrows=as_counts(z["early"])
    keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for r in rawrows}
    built,_=s.read_source_builtup(a.built);land,_=s.source_continents(a.ne);ras,mem,_=s.worldclim_tif(a.climate)
    try:
        scale,_=s.climate_scale(ras);clim,_=s.source_climate_lookup(land,ras,scale,keys,W)
    finally:ras.close();mem.close()
    tr=s.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    rd=s.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    pp=s.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    el={c:{(tuple(v["a"]),tuple(v["b"])):v["crossing_source_segments"] for v in edge["edges"][str(W)][c]} for c in s.CATS}
    rows,_,_=s.original_cell_rows(rawrows,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
    rows,_=s.inject_built(rows,built[str(W)]);rows,_=s.inject_status(rows,status,W)
    panels=[]
    for gate in (.5,.8):
        rr,_=s.enrich_cells(rows,rd,tr,W,gate);rr,_=s.join_population(rr,pp,rd)
        ans=analyze({r["key"]:r for r in rr},z["early"],z["late"],el,traitmap,gate)
        panels.append(ans);print(json.dumps({"bird_neighborhood_record_panel":ans},allow_nan=False),flush=True)
    out={"version":"global-rail-bird-neighborhood-source-Nov2026-v1","frozen_environment_replay":checks,
         "monthly_source":snapshot.files,"old_archive_head":snapshot.old_archive,"taxa_sha256":sha,"taxonomy_qa":txqa,"trait_match_qa":traitqa,"monthly_census_drift":drift,
         "obs_sha256":osha,"trait_sha256":p.SOURCE_SHA,"exact_global_2025_count_replay":False,"new_snapshot_internal_count_replay":True,
         "scan_qa":scanqa,"panels":panels,
         "primary":"200km land50 baseline-supported world frame, target x habitat/migration x number documented source cells and species FE, low-minus-high HWI x active-minus-inactive source exposure, new-record outcome",
         "limits":["all available baseline-supported same-block adjacent neighbors, not all potential source populations or known origins",
          "new record not colonization or actual movement; baseline nonrecord not true absence",
          "source documented in both observer halves; follow target selected on effort only",
          "redetection separate diagnostic, not formal difference in coefficients",
          "2026 map status cannot establish actual 2018-25 operation, status and habitat/access confounding remain",
          "exploratory distinct estimand after unique-source support failure; no threshold or window tuning",
          "November snapshot changes records and taxonomy; do not combine old and new or treat this as independent biological replication",
          "no eBird, individual records, species IDs, observer IDs or fine coordinates exported"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
