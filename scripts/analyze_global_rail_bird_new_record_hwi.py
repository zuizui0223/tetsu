#!/usr/bin/env python3
"""Global bird prospective new-record test: HWI x mapped railway crossing.
Retrospective equal-length periods, unique documented source among eligible
neighbors, not observed movement or verified colonization. No individual IDs exported.
"""
import argparse, collections, hashlib, json, math, re
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t
from pyproj import Transformer
import analyze_world_railstatus_observer_split_connectivity as s
from analyze_world_railstatus_bird_historical_hwi import cells_from_pairs
from scan_world_railstatus_bird_conditional_hwi import (
    taxonomy, source_client, open_original_source, source_coarse_cell,
    OBS_ETAG, OBS_BYTES, MAX_UNCERT_M, PROGRESS_EVERY_SOURCE_ROWS,
    bucket, add, pair_rows, SALT, SOURCE_SHA, DERIVATIVE, request, names,
    conditioned_match, validate, csv_rows, W, pair_key)
from scan_world_railstatus_observer_split_beta import full_counts, split_counts, eligible, MIN_FULL, MIN_SPLIT

EARLY=tuple(range(2018,2022)); LATE=tuple(range(2022,2026))
MIN_GUILD=5
MIN_BLOCKS=20
MIN_EDGES=100
STATES=tuple(s.CROSS)
ENV=("bio1","tree_consensus_fullcell","log_road_density","log_population_density",
     "log_active_normal","log_inactive_normal")
# Edge-stratum intercept absorbs environmental main effects; retain
# prespecified HWI-dependent environmental differences and other transport crossings.
CONTEXT=tuple(pre+k for k in ENV for pre in ("mean_","diff_"))
EFFORT=("source_base_logobs","target_base_logobs","target_follow_logobs",
        "source_base_logobservers","target_base_logobservers","target_follow_logobservers")
FEATURES=STATES+CONTEXT+EFFORT

def period(date):
    if not isinstance(date,str) or len(date)<10 or date[4]!="-" or date[7]!="-" or not date[:4].isdigit():return None
    y=int(date[:4])
    return "early" if y in EARLY else "late" if y in LATE else None

def scan(client,canonical,cells):
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    z={k:collections.defaultdict(bucket) for k in ("early","late","2025")}
    qa=collections.Counter()
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:
                print(json.dumps({"scan_progress":dict(qa)}),flush=True)
            y=period(r["observed_on"])
            if y is None or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:
                hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"])
                lon,lat=float(r["longitude"]),float(r["latitude"])
            except (ValueError,TypeError):continue
            if hit is None or hit[0]!="birds" or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k not in cells:continue
            o=int(r["observer_id"]) if r["observer_id"] else None
            add(z[y][k],hit[1],o);qa[y+"_bird_records"]+=1
            if r["observed_on"].startswith("2025-"):add(z["2025"][k],hit[1],o)
    finally:stream.close()
    if qa["rows"]!=288721938 or raw.bytes!=OBS_BYTES:raise ValueError("incomplete original source")
    return z,dict(qa),raw.hash.hexdigest()

def supported(v):
    return v is not None and eligible(full_counts(v),MIN_FULL) and all(eligible(c,MIN_SPLIT) for c in split_counts(v))

def species_strata(traits):
    out={}
    for tid,gg in traits.items():
        low="low_HWI" in gg
        if not (low or "high_HWI" in gg):continue
        strata=[x for x in gg if re.fullmatch(r"h[123]_m[123]_(low|high)",x)]
        if len(strata)==1:
            g=strata[0].rsplit("_",1)[0]
            if strata[0].endswith("_low")!=low:raise ValueError("HWI stratum mismatch")
            out[tid]=(int(low),g)
    return out

def candidates(edges,base,follow,traits):
    """Outcome-free source risk set; require >=2 eligible source neighbors.
    A documented source is unique only among eligible mapped-frame neighbors.
    Baseline presence in both source halves, absence in ALL target baseline
    records, absence from every other eligible source neighbor's baseline records.
    """
    nb=collections.defaultdict(list)
    for e in edges:
        a,b=e["a"],e["b"]
        if not all(supported(z.get(k)) for z in (base,follow) for k in (a,b)):continue
        nb[a].append((b,e));nb[b].append((a,e))
    cases=[];qa=collections.Counter()
    for dst,neighbors in sorted(nb.items()):
        if len(neighbors)<2:continue
        qa["eligible_targets"]+=1
        nsource=collections.Counter()
        for src,e in neighbors:nsource.update(base[src][1])
        missing=base[dst][1];later=follow[dst]
        for src,e in neighbors:
            group=collections.defaultdict(list)
            eb=base[src];known=eb[4][0]&eb[4][1]
            for tid in known:
                if tid in missing or nsource[tid]!=1 or tid not in traits:continue
                low,stratum=traits[tid]
                group[(stratum,low)].append(tid)
            for stratum in sorted({g for g,l in group}):
                lo,hi=group[stratum,1],group[stratum,0]
                if min(len(lo),len(hi))<MIN_GUILD:
                    qa["edge_strata_insufficient_guilds"]+=1;continue
                ec=(src,dst,stratum)
                effort={}
                for label,v in (("source_base",eb),("target_base",base[dst]),("target_follow",later)):
                    obs,observers,_=full_counts(v)
                    effort[label+"_logobs"]=math.log1p(obs);effort[label+"_logobservers"]=math.log1p(observers)
                f=[float(e[k]) for k in STATES+CONTEXT]+[effort[k] for k in EFFORT]
                for low,ids in ((1,lo),(0,hi)):
                    for tid in ids:
                        a0=tid in later[4][0];a1=tid in later[4][1]
                        # Each later target observer side is disjoint from one
                        # independently documented baseline source side.
                        cases.append({"tid":tid,"low":low,"stratum":stratum,"edge":(src,dst),
                                      "edge_stratum":ec,"src":src,"dst":dst,
                                      "block":(e["continent"],e["block"]),"x":[low*v for v in f],
                                      "y":(int(a0)+int(a1))/2,"both":float(a0 and a1),
                                      "any":float(a0 or a1),"side0":float(a0),"side1":float(a1)})
                qa["eligible_edge_strata"]+=1
    return cases,dict(qa)

def codes(values):
    m={}
    return np.asarray([m.setdefault(x,len(m)) for x in values],dtype=np.int64),len(m)

def demean(A,w,groups):
    """Alternating weighted projections onto fixed-effect orthogonal spaces."""
    A=A.copy()
    for _ in range(1000):
        change=0.
        for ix,n in groups:
            den=np.bincount(ix,weights=w,minlength=n)
            for j in range(A.shape[1]):
                mu=np.bincount(ix,weights=w*A[:,j],minlength=n)/den
                change=max(change,float(np.max(abs(mu))))
                A[:,j]-=mu[ix]
        if change<1e-9:return A
    raise ValueError("fixed-effect projections did not converge")

def fit(rows,endpoint=False,target="y"):
    bs={r["block"] for r in rows};edges={r["edge"] for r in rows}
    qa={"n_cases":len(rows),"n_directed_edges":len(edges),"n_blocks":len(bs),
        "n_species":len({r["tid"] for r in rows}),
        "by_HWI":dict(collections.Counter("low" if r["low"] else "high" for r in rows))}
    if len(edges)<MIN_EDGES or len(bs)<MIN_BLOCKS:return dict(qa,status="insufficient support")
    # One block weight, one directed edge weight inside block; one stratum
    # weight inside edge; equal weights for low/high within the stratum.
    block_edges=collections.defaultdict(set);edge_strata=collections.defaultdict(set)
    count=collections.Counter()
    variation=collections.defaultdict(lambda:[set(),set()])
    for r in rows:
        block_edges[r["block"]].add(r["edge"])
        edge_strata[r["edge"]].add(r["stratum"])
        count[r["edge_stratum"],r["low"]]+=1
        if r["low"]:
            variation[r["block"]][0].add(r["x"][0]);variation[r["block"]][1].add(r["x"][1])
    nv=[sum(len(v[i])>1 for v in variation.values()) for i in (0,1)]
    qa["varying_blocks_active_inactive"]=nv
    if min(nv)<MIN_BLOCKS:return dict(qa,status="insufficient status variation")
    w=np.asarray([1/(len(block_edges[r["block"]])*len(edge_strata[r["edge"]])*2*count[r["edge_stratum"],r["low"]]) for r in rows],float)
    w/=w.sum()
    groups=[codes(r["edge_stratum"] for r in rows),codes(r["tid"] for r in rows)]
    if endpoint:
        groups += [codes((r["src"],r["low"]) for r in rows),codes((r["dst"],r["low"]) for r in rows)]
    mat=np.asarray([r["x"]+[r[target]] for r in rows],float)
    try:mat=demean(mat,w,groups)
    except ValueError as e:return dict(qa,status="projection failed",reason=str(e))
    mat*=np.sqrt(w)[:,None]
    A,y=mat[:,:-1],mat[:,-1]
    norm=np.linalg.norm(A,axis=0);keep=norm>1e-10
    fn=[k for k,v in zip(FEATURES,keep) if v]
    if not set(STATES[:2])<=set(fn):return dict(qa,status="no residual status contrast")
    Z=A[:,keep]/norm[keep]
    sv=np.linalg.svd(Z,compute_uv=False)
    if sv[-1]<1e-8*sv[0]:return dict(qa,status="rank deficient")
    inv=np.linalg.inv(Z.T@Z);beta=inv@(Z.T@y);res=y-Z@beta
    # Direct low-HWI x (active-inactive) coefficient, with block and
    # species clustering. Report both one-way and conservative max variance
    # including two-way intersection; no pseudo-replication at species-cell level.
    def result(c):
        est=float(c@(beta/norm[keep]))
        influence=(Z@(inv@(c/norm[keep])))*res
        def vv(values):
            ix,n=codes(values)
            if n<2:return 0.
            sc=np.bincount(ix,weights=influence,minlength=n)
            return float(np.dot(sc,sc)*n/(n-1))
        vb=vv(r["block"] for r in rows);vs=vv(r["tid"] for r in rows)
        vi=vv((r["block"],r["tid"]) for r in rows)
        se=math.sqrt(max(vb,vs,vb+vs-vi,0))
        crit=float(student_t.ppf(.975,len(bs)-1))
        return {"estimate":est,"interval95_approx":[est-crit*se,est+crit*se],
                "block_cluster_se":math.sqrt(vb),"species_cluster_se":math.sqrt(vs),
                "conservative_cluster_se":se}
    a=np.array([float(k==STATES[0]) for k in fn]);i=np.array([float(k==STATES[1]) for k in fn])
    return dict(qa,status="estimated",effects={"low_HWI_x_active":result(a),
        "low_HWI_x_inactive":result(i),"low_HWI_x_active_minus_inactive":result(a-i)})

def analyze(edges,z,traits):
    out={}
    for label,first,second in (("forward","early","late"),("reverse_diagnostic","late","early")):
        rows,qa=candidates(edges,z[first],z[second],traits)
        summary={"candidate_qa":qa,"n_cases":len(rows),"primary":fit(rows),
                 "endpoint_HWI_FE":fit(rows,endpoint=True),"both_observer_sides":fit(rows,target="both"),
                 "either_observer_side":fit(rows,target="any"),
                 "resident_only":fit([r for r in rows if r["stratum"].endswith("_m1")])}
        # Do not select a positive annual window or threshold after seeing results.
        out[label]=summary
    return out

def selftest():
    assert W==200000
    assert STATES[:2]==("cross_active_normal","cross_inactive_normal")
    assert period("2018-01-01")=="early" and period("2022-01-01")=="late" and period("2017-01-01") is None
    assert species_strata({1:("low_HWI","h2_m1_low"),2:("high_HWI","h2_m1_high")})=={1:(1,"h2_m1"),2:(0,"h2_m1")}
    v=np.array([[1.,2.],[3.,4.],[5.,6.],[7.,8.]])
    r=demean(v,np.ones(4),[codes([0,0,1,1])])
    assert np.max(abs(r-np.array([[-1.,-1.],[1.,1.],[-1.,-1.],[1.,1.]])))<1e-12
    assert fit([])["status"]=="insufficient support"

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true")
    for n in ("conditional","annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path)
    a=p.parse_args();selftest()
    if a.self_test:return
    cond,annual,src,prior,edge,plant,rail,status,tree,road,pop=map(s.read,(a.conditional,a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=s.verify_original_sources(plant,rail);s.verify_replay(prior,plant)
    if cond["version"]!="railstatus-avonet-conditional-hwi-pairs-v1" or not cond["exact_full_pair_replay"] or not cond["exact_hwi_pair_replay"] or cond["source"]["source_SHA256"]!=SOURCE_SHA:
        raise ValueError("unverified HWI source")
    if any(len({cond[k],annual[k],src[k]})!=1 for k in ("obs_sha256","taxa_sha256","split_seed")) or src["split_seed"]!=SALT.decode():
        raise ValueError("prior source mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in s.CATS):
        raise ValueError("map geometry mismatch")
    if tree["source_file_SHA256"]!=s.EXPECTED_SHA or road["source_zip_MD5"]!=s.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":
        raise ValueError("environment source mismatch")
    if s.sha_file(a.ne)!=s.ORIGINAL_NE_SHA or s.sha_file(a.climate)!=s.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=s.POP_SHA or pop["relative_population_sum_difference"]>.01:
        raise ValueError("environment source mismatch")
    previous={y:{pair_key(p):p for p in annual["pairs"][y]["birds"]} for y in ("2023","2024")}
    previous["2025"]={pair_key(p):p for p in src["pairs"][str(W)]["birds"] if "split_similarity" in p}
    common=set.intersection(*(set(v) for v in previous.values()))
    if len(common)!=1645 or any({pair_key(p) for p in cond["pairs"][y]}!=common for y in previous):
        raise ValueError("frozen candidate geography changed")
    cells={k for pair in common for k in pair}
    raw=request(DERIVATIVE,cap=10_000_000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError("trait source changed")
    client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=src["taxa_sha256"]:raise ValueError("taxonomy changed")
    joined,traitqa=conditioned_match(validate(csv_rows(raw)),names(client,canonical,sha))
    if traitqa!=cond["source"]["match"]:raise ValueError("trait membership replay failed")
    traits=species_strata(joined)
    print(json.dumps({"frozen_frame":{"pairs":len(common),"cells":len(cells),"trait_species":len(traits),"windows":[EARLY,LATE]}}),flush=True)
    z,scanqa,osha=scan(client,canonical,cells)
    if osha!=src["obs_sha256"]:raise ValueError("observation SHA mismatch")
    replay={pair_key(p):p for p in pair_rows(z["2025"]) if pair_key(p) in common and "split_similarity" in p}
    if set(replay)!=common or any(json.loads(json.dumps(p))!=previous["2025"][k] for k,p in replay.items()):
        raise ValueError("exact original 2025 counts/overlap replay failed")
    built,_=s.read_source_builtup(a.built);land,_=s.source_continents(a.ne);ras,mem,_=s.worldclim_tif(a.climate)
    try:
        scale,_=s.climate_scale(ras);clim,_=s.source_climate_lookup(land,ras,scale,cells,W)
    finally:ras.close();mem.close()
    tr=s.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    rd=s.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    pp=s.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    el={c:{(tuple(v["a"]),tuple(v["b"])):v["crossing_source_segments"] for v in edge["edges"][str(W)][c]} for c in s.CATS}
    raw=[replay[k] for k in sorted(common)];panels=[]
    for gate in (.5,.8):
        rows,_,_=s.original_cell_rows(cells_from_pairs(raw),rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
        rows,_=s.inject_built(rows,built[str(W)]);rows,_=s.inject_status(rows,status,W)
        rows,_=s.enrich_cells(rows,rd,tr,W,gate);rows,_=s.join_population(rows,pp,rd)
        ee=s.build_pairs([s.adapt_pair(p) for p in raw],rows,el)
        if len(ee)!={.5:819,.8:649}[gate]:raise ValueError("original geography replay failed")
        ans={"land":gate,"n_frozen_edges":len(ee),"results":analyze(ee,z,traits)}
        panels.append(ans);print(json.dumps({"bird_temporal_new_record_panel":ans},allow_nan=False),flush=True)
    out={"version":"global-rail-HWI-prospective-record-v1","source_replay":checks,"scan_qa":scanqa,
        "source_SHA256":SOURCE_SHA,"obs_sha256":osha,"panels":panels,
        "primary":"200km land50, 2018-21 documented source and baseline target nonrecord, 2022-25 target new record, low-HWI x active-minus-inactive; edge x habitat/migration strata and species FE",
        "limits":["new record not colonization, detection nonrecord not true absence, no actual movement tracking",
            "unique documented source among eligible neighbors only; not a known origin or all possible sources",
            "both-endpoint observation support in both periods and 2025 geography selected retrospectively",
            "fixed 2026 railway status map cannot establish historical active operation or closure",
            "observer halves, same edge habitat/migration strata and species effects reduce but do not eliminate detection, geographic range and habitat confounding",
            "reverse-time diagnostic is not a causal placebo; no further threshold/window tuning; no raw IDs or coordinates exported; no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
