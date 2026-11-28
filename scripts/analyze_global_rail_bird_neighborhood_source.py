#!/usr/bin/env python3
"""Global neighborhood-source exposure and subsequent bird records.

Different estimand from unique source per directed edge. Candidate frame is
baseline-supported globally, including cells with no active rail. No origin,
colonization, or movement path is observed. All identities stay in runner.
"""
import argparse
import collections
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import t
from pyproj import Transformer
import analyze_global_rail_bird_new_record_hwi as p

s=p.s
W=200000
MIN_TARGETS=100
MIN_BLOCKS=20
MIN_GUILD=5
CROSS=tuple(s.CROSS)
SOURCE_FIELDS=("bio1","tree_consensus_fullcell","log_road_density","log_population_density",
               "log_active_normal","log_inactive_normal","log_observations","log_observers")
TARGET_FIELDS=("bio1","tree_consensus_fullcell","log_road_density","log_population_density",
               "log_active_normal","log_inactive_normal","log_observations","log_observers",
               "follow_logobs","follow_logobservers")
# Source-mean and environmental distance to target, relative to all
# eligible neighboring source cells. Status exposure is likewise centered.
NAMES=CROSS+tuple("source_mean_"+f for f in SOURCE_FIELDS)+tuple("source_diff_"+f for f in SOURCE_FIELDS[:6])
FIELDS=NAMES+tuple("low_x_"+f for f in NAMES)+tuple("low_x_target_"+f for f in TARGET_FIELDS)
A_NAME="low_x_cross_active_normal"
I_NAME="low_x_cross_inactive_normal"

def scan(client,canonical):
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    z={k:collections.defaultdict(p.bucket) for k in ("early","late","2025")}
    qa=collections.Counter()
    rr,raw,stream=p.open_original_source(client,"observations.csv.gz",p.OBS_ETAG,p.OBS_BYTES)
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
    if qa["rows"]!=288721938 or raw.bytes!=p.OBS_BYTES:
        raise ValueError("incomplete source")
    return z,dict(qa),raw.hash.hexdigest()

def later_effort(v):
    if v is None:return False
    # Follow-up observation effort only, NOT its species richness or
    # any individual candidate species' follow-up outcome.
    return p.full_counts(v)[0]>=50 and p.full_counts(v)[1]>=5 and all(
        a>=25 and b>=5 for a,b,_ in p.split_counts(v))

def as_counts(base):
    out=[]
    for k,v in sorted(base.items()):
        if not p.supported(v):continue
        o,b,n=p.full_counts(v)
        out.append(s.convert({"source_equal_area_grid_ix":k[0],"source_equal_area_grid_iy":k[1],
                              "observations":o,"observers":b,"species":n}))
    return out

def neighbors(dst,geo,edges):
    x=geo[dst]
    ans=[]
    for j in ((dst[0]-1,dst[1]),(dst[0]+1,dst[1]),(dst[0],dst[1]-1),(dst[0],dst[1]+1)):
        if j not in geo:continue
        y=geo[j]
        if x["continent"]!=y["continent"] or x["block"]!=y["block"]:continue
        k=tuple(sorted((dst,j)))
        flags=[float(edges[c].get(k,0)>0) for c in s.CATS]
        for i,c in enumerate(s.CATS):
            if flags[i] and min(x[c+"_km"],y[c+"_km"])<=0:
                raise ValueError("source crossing inconsistent with cell km")
        values=flags+[y[f] for f in SOURCE_FIELDS]+[abs(y[f]-x[f]) for f in SOURCE_FIELDS[:6]]
        ans.append((j,np.asarray(values,dtype=float)))
    return ans

def risk_sets(geo,base,follow,edges,traits):
    by_out={"new_record":[],"redetection":[]};qa=collections.Counter()
    for dst,x in sorted(geo.items()):
        if not later_effort(follow.get(dst)):
            qa["target_follow_effort_insufficient"]+=1;continue
        nb=neighbors(dst,geo,edges)
        if len(nb)<2:
            qa["fewer_than_two_baseline_neighbors"]+=1;continue
        V=np.asarray([v for j,v in nb])
        # At least one conventional railway flag varies between eligible
        # neighbors; no filter on later species or coefficient.
        if not np.any(np.ptp(V[:,:2],axis=0)>0):
            qa["no_neighbor_status_variation"]+=1;continue
        qa["targets_with_source_opportunity"]+=1
        allmean=V.mean(axis=0)
        documented=collections.defaultdict(list)
        for i,(j,v) in enumerate(nb):
            for tid in base[j][4][0]&base[j][4][1]:
                if tid in traits:documented[tid].append(i)
        b=base[dst];f=follow[dst];known=b[4][0]&b[4][1]
        target=dict(x,follow_logobs=math.log1p(p.full_counts(f)[0]),
                    follow_logobservers=math.log1p(p.full_counts(f)[1]))
        target_v=np.asarray([target[k] for k in TARGET_FIELDS],float)
        strata=collections.defaultdict(lambda:{0:[],1:[]})
        for tid,ix in documented.items():
            # All-neighbor occupancy cannot distinguish connected sources.
            if not 0<len(ix)<len(nb):continue
            kind="redetection" if tid in known else "new_record" if tid not in b[1] else None
            if kind is None:continue
            low,hab=traits[tid]
            # Same number of documented source cells in the FE stratum.
            q=V[ix].mean(axis=0)-allmean
            strata[kind,hab,len(ix)][low].append((tid,q))
        for (kind,hab,nsource),guilds in sorted(strata.items()):
            if min(map(len,guilds.values()))<MIN_GUILD:
                qa[kind+"_insufficient_strata"]+=1;continue
            qa[kind+"_eligible_strata"]+=1
            for low,items in guilds.items():
                for tid,q in items:
                    one,two=tid in f[4][0],tid in f[4][1]
                    by_out[kind].append({
                        "tid":tid,"low":low,"target":dst,
                        "stratum":(dst,hab,nsource),"hab":hab,
                        "block":(x["continent"],x["block"]),
                        "x":np.concatenate([q,low*q,low*target_v]),
                        "y":(int(one)+int(two))/2,"both":float(one and two),
                        "any":float(one or two)})
    return by_out,dict(qa)

def fit(rows,outcome="y"):
    targets={r["target"] for r in rows}
    bs={r["block"] for r in rows}
    species={r["tid"] for r in rows}
    qa={"n_cases":len(rows),"n_targets":len(targets),"n_blocks":len(bs),
        "n_species":len(species),"by_HWI":dict(collections.Counter("low" if r["low"] else "high" for r in rows))}
    if len(targets)<MIN_TARGETS or len(bs)<MIN_BLOCKS:
        return dict(qa,status="insufficient geographic support")
    block_targets=collections.defaultdict(set);target_strata=collections.defaultdict(set)
    counts=collections.Counter()
    varies=collections.defaultdict(lambda:[set(),set()])
    for r in rows:
        block_targets[r["block"]].add(r["target"])
        target_strata[r["target"]].add(r["stratum"])
        counts[r["stratum"],r["low"]]+=1
        # Within-stratum exposure support rather than overall status counts.
        for k in (0,1):varies[r["block"],r["stratum"]][k].add(round(float(r["x"][k]),10))
    vb=[{block for (block,g),v in varies.items() if len(v[k])>1} for k in (0,1)]
    qa["blocks_with_within_stratum_source_variation"]=list(map(len,vb))
    if min(map(len,vb))<MIN_BLOCKS:
        return dict(qa,status="insufficient status source variation")
    w=np.array([1/(len(block_targets[r["block"]])*len(target_strata[r["target"]])*2*counts[r["stratum"],r["low"]]) for r in rows])
    w/=w.sum()
    groups=[p.codes(r["stratum"] for r in rows),p.codes(r["tid"] for r in rows)]
    if len(rows)-sum(n for ix,n in groups)-len(FIELDS)<30:
        return dict(qa,status="insufficient residual support")
    mat=np.column_stack((np.asarray([r["x"] for r in rows]),np.asarray([r[outcome] for r in rows])))
    try:mat=p.demean(mat,w,groups)
    except ValueError as e:return dict(qa,status="projection failed",reason=str(e))
    mat*=np.sqrt(w)[:,None]
    A,y=mat[:,:-1],mat[:,-1]
    norm=np.linalg.norm(A,axis=0);keep=norm>1e-10
    fn=[k for k,v in zip(FIELDS,keep) if v]
    if not {A_NAME,I_NAME,"cross_active_normal","cross_inactive_normal"}<=set(fn):
        return dict(qa,status="no residual status contrast")
    Z=A[:,keep]/norm[keep];sv=np.linalg.svd(Z,compute_uv=False)
    if sv[-1]<1e-8*sv[0]:
        return dict(qa,status="rank deficient")
    inv=np.linalg.inv(Z.T@Z);beta=inv@(Z.T@y);res=y-Z@beta
    bcode=p.codes(r["block"] for r in rows)
    scode=p.codes(r["tid"] for r in rows)
    icode=p.codes((r["block"],r["tid"]) for r in rows)
    def result(terms):
        c=np.asarray([terms.get(k,0.) for k in fn])
        est=float(c@(beta/norm[keep]))
        inf=(Z@(inv@(c/norm[keep])))*res
        def variance(code):
            ix,n=code
            sc=np.bincount(ix,weights=inf,minlength=n)
            return float(np.dot(sc,sc)*n/(n-1))
        vb,vs,vi=map(variance,(bcode,scode,icode))
        se=math.sqrt(max(vb,vs,vb+vs-vi,0.))
        h=float(t.ppf(.975,len(bs)-1))*se
        return {"estimate":est,"interval95_approx":[est-h,est+h],
                "block_cluster_SE":math.sqrt(vb),"species_cluster_SE":math.sqrt(vs),
                "conservative_SE":se}
    return dict(qa,status="estimated",effects={
        "high_HWI_active":result({"cross_active_normal":1.}),
        "high_HWI_inactive":result({"cross_inactive_normal":1.}),
        "high_HWI_status_difference":result({"cross_active_normal":1.,"cross_inactive_normal":-1.}),
        "low_HWI_active":result({"cross_active_normal":1.,A_NAME:1.}),
        "low_HWI_inactive":result({"cross_inactive_normal":1.,I_NAME:1.}),
        "low_HWI_status_difference":result({"cross_active_normal":1.,"cross_inactive_normal":-1.,A_NAME:1.,I_NAME:-1.}),
        "low_minus_high_active":result({A_NAME:1.}),
        "low_minus_high_inactive":result({I_NAME:1.}),
        "low_minus_high_active_minus_inactive":result({A_NAME:1.,I_NAME:-1.})})

def analyze(geo,base,follow,edges,traits,gate):
    rr,qa=risk_sets(geo,base,follow,edges,traits)
    # No tuning or replacement. Redetection is a separate diagnostic,
    # not a formal difference against the new-record coefficient.
    ans={"land":gate,"baseline_environment_cells":len(geo),"risk_set_qa":qa,
         "new_record":{"primary":fit(rr["new_record"]),
          "both_observer_sides":fit(rr["new_record"],"both"),
          "resident_only":fit([r for r in rr["new_record"] if r["hab"].endswith("_m1")])},
         "redetection_diagnostic":{"primary":fit(rr["redetection"])}}
    return ans

def selftest():
    assert W==p.W and CROSS[:2]==("cross_active_normal","cross_inactive_normal")
    assert len(set(FIELDS))==len(FIELDS)
    assert fit([])["status"]=="insufficient geographic support"
    a=p.bucket()
    for i in range(100):p.add(a,1,i)
    assert later_effort(a) and not p.supported(a)
    v=np.array([[1.,0.],[0.,1.],[0.,0.]])
    assert np.allclose(v[[0]].mean(axis=0)-v.mean(axis=0),[2/3,-1/3])
    assert np.allclose(v.mean(axis=0)-v.mean(axis=0),[0,0])

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
    client=p.source_client();canonical,_,sha=p.taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("taxonomy changed")
    joined,traitqa=p.conditioned_match(p.validate(p.csv_rows(raw)),p.names(client,canonical,sha))
    if traitqa!=cond["source"]["match"]:raise ValueError("trait membership replay failed")
    traits=p.species_strata(joined)
    if len(traits)!=9855:raise ValueError("trait strata census changed")
    z,scanqa,osha=scan(client,canonical)
    if osha!=src["obs_sha256"] or scanqa["2025_bird_records"]!=prior["qa"]["birds"]:
        raise ValueError("exact original occurrence source failed")
    old={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]):(r["observations"],r["observers"],r["species"]) for r in prior["grids"][str(W)]["birds"]}
    new={k:p.full_counts(v) for k,v in z["2025"].items()}
    if new!=old:raise ValueError("full global 2025 bird grid counts differ")
    print(json.dumps({"full_2025_grid_replay":len(new),"eligible_baseline_cells":sum(p.supported(v) for v in z["early"].values()),"scan_qa":scanqa}),flush=True)
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
        ans=analyze({r["key"]:r for r in rr},z["early"],z["late"],el,traits,gate)
        panels.append(ans);print(json.dumps({"bird_neighborhood_record_panel":ans},allow_nan=False),flush=True)
    out={"version":"global-rail-bird-neighborhood-source-v1","source_replay":checks,
         "obs_sha256":osha,"trait_sha256":p.SOURCE_SHA,"exact_global_2025_count_replay":True,
         "scan_qa":scanqa,"panels":panels,
         "primary":"200km land50 baseline-supported world frame, target x habitat/migration x number documented source cells and species FE, low-minus-high HWI x active-minus-inactive source exposure, new-record outcome",
         "limits":["all available baseline-supported same-block adjacent neighbors, not all potential source populations or known origins",
          "new record not colonization or actual movement; baseline nonrecord not true absence",
          "source documented in both observer halves; follow target selected on effort only",
          "redetection separate diagnostic, not formal difference in coefficients",
          "2026 map status cannot establish actual 2018-25 operation, status and habitat/access confounding remain",
          "exploratory distinct estimand after unique-source support failure; no threshold or window tuning",
          "no eBird, individual records, species IDs, observer IDs or fine coordinates exported"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
