#!/usr/bin/env python3
"""Posthoc mechanism probe: mapped active versus inactive rail association
with GloNAF-documented naturalized vs unlisted plant record richness.
Reuse exact v397 three-year coarse counts, no new raw records.
"""
import argparse,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from analyze_world_railstatus_population_local import (
    FULL,POP,ACTIVE,INACTIVE,POP_SHA,join_population,inject_status,
    enrich_cells,make_lookup,read,MD5,EXPECTED_SHA,block_weights)
from analyze_global_naturalized_count_decomposition_v389 import source_rows,count_outcomes
from analyze_world_railstatus_cross_taxa_population_local import flexible,EXTRA

W=100000
YEARS=("2023","2024","2025")
OUT=("log1p_documented","log1p_unlisted","log1p_count_contrast")
FIELDS=FULL+(POP,)
TARGETS={"active":{ACTIVE:1.},"inactive":{INACTIVE:1.},"active_minus_inactive":{ACTIVE:1.,INACTIVE:-1.}}
def strata_key(r,local):
    return (r["region"],r["continent"],r["block"]) if local else (r["region"],)
def supported_local(rows):
    d=defaultdict(list)
    for r in rows:d[strata_key(r,True)].append(r)
    return [r for v in d.values() if len(v)>=2 for r in v]
def fit(rows,fields=FIELDS,local=True):
    gs=defaultdict(list);bs=defaultdict(list)
    for i,r in enumerate(rows):
        gs[strata_key(r,local)].append(i);bs[(r["continent"],r["block"])].append(i)
    qa={"n":len(rows),"n_blocks":len(bs),"n_strata":len(gs),"n_regions":len({r["region"] for r in rows}),
        "continents":dict(Counter(r["continent"] for r in rows)),
        "n_blocks_active_varies":sum(len({rows[i][ACTIVE] for i in v})>1 for v in bs.values()),
        "n_blocks_inactive_varies":sum(len({rows[i][INACTIVE] for i in v})>1 for v in bs.values())}
    if len(bs)<20 or len(rows)-len(gs)-len(fields)<30:return dict(qa,status="insufficient support")
    w=block_weights(rows)
    X=np.array([[r[f] for f in fields+OUT] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite source")
    for ix in gs.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None];A=X[:,:len(fields)];Y=X[:,len(fields):]
    scale=np.sqrt((A*A).sum(axis=0));ok=scale>1e-10
    names=[f for f,k in zip(fields,ok) if k]
    if not {ACTIVE,INACTIVE}<=set(names):return dict(qa,status="no status variation")
    Z=A[:,ok]/scale[ok];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return dict(qa,status="rank deficient")
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@Y);e=Y-Z@b
    coef=b/scale[ok,None]
    scores=np.array([(inv@(Z[ix].T@e[ix]))/scale[ok,None] for ix in bs.values()])
    effects={}
    for t,terms in TARGETS.items():
        c=np.array([terms.get(f,0.) for f in names])
        beta=c@coef;inf=np.einsum("i,bio->bo",c,scores)
        if abs(beta[0]-beta[1]-beta[2])>1e-9:raise ValueError("Count contrast identity failed")
        se=np.sqrt(len(bs)/(len(bs)-1)*(inf**2).sum(axis=0))
        effects[t]={o:{"estimate":float(v),"interval95_approx":[float(v-1.96*z),float(v+1.96*z)]} for o,v,z in zip(OUT,beta,se)}
    return dict(qa,status="estimated",effects=effects,condition_number=float(s[0]/s[-1]))
def period_model(rows,period):
    local=supported_local(rows)
    high=supported_local([r for r in rows if r["observations"]>=(150 if period=="2023-2025" else 200) and r["observers"]>=(15 if period=="2023-2025" else 20)])
    both=supported_local([r for r in rows if r["active_normal_km"]>0 and r["inactive_normal_km"]>0])
    return {"local":fit(local),"local_flexible_effort":fit(local,FIELDS+EXTRA),
            "region_FE":fit(rows,local=False),"high_effort":fit(high),"both_status_positive":fit(both),
            "leave_one_continent_out":{c:fit(supported_local([r for r in local if r["continent"]!=c])) for c in sorted({r["continent"] for r in local})}}
def normalize(raw):
    r=dict(raw,key=tuple(raw["key"]),block=tuple(raw["block"]))
    y=count_outcomes(r)[:3]
    for o,v in zip(OUT,y):
        if abs(r[o]-v)>1e-12:raise ValueError("v397 count contrast replay mismatch")
    return r
def main():
    p=argparse.ArgumentParser()
    for k in ("annual","v383","status","road","tree","population","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    annual,old,status,road,tree,pop=map(read,(a.annual,a.v383,a.status,a.road,a.tree,a.population))
    replay=annual["source_replay"]
    if annual["version"]!="v397-original-2023-2025-coarse-aggregates" or replay["all_rows"]!=288721938 or not replay["original_2025_all_cell_counts_exact"] or replay["source_observation_gzip_sha256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279":raise ValueError("Wrong annual source")
    if old["version"]!="v383-exact-v379d-coarse-cell-aggregate" or old["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:raise ValueError("Wrong v383 source")
    if not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong status source")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or tree["source_file_SHA256"]!=EXPECTED_SHA:raise ValueError("Wrong road/tree source")
    if pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong population source")
    rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    prior={z["region_guard_degrees"]:source_rows(z) for z in old["panels"] if z["width_m"]==W}
    src={(z["period"],z["region_guard_degrees"]):[normalize(r) for r in z["rows"]] for z in annual["panels"]}
    if len(src)!=8:raise ValueError("Incomplete annual panels")
    for pad in (.3,0.):
        aa={r["key"]:r for r in src["2025",pad]};bb={r["key"]:r for r in prior[pad]}
        if set(aa)!=set(bb):raise ValueError("2025 support mismatch")
        for k in aa:
            for f in ("region","block","continent","observations","observers","species","vascular_richness","documented_richness"):
                if aa[k][f]!=bb[k][f]:raise ValueError("2025 replay mismatch "+f)
    panels=[]
    for pad in (.3,0.):
        for land in (.5,.8):
            sets={}
            for period in (*YEARS,"2023-2025"):
                rows,_=inject_status(src[period,pad],status,W)
                rows,_=enrich_cells(rows,rd,tr,W,land);rows,_=join_population(rows,pp,rd)
                sets[period]=[flexible(r) for r in rows]
            common=set.intersection(*({r["key"] for r in sets[y]} for y in YEARS))
            common_local={r["key"] for r in supported_local([r for r in sets["2025"] if r["key"] in common])}
            z={"pad":pad,"land":land,"n_common_local":len(common_local),
               "native_period_support":{y:period_model(rr,y) for y,rr in sets.items()},
               "same_cells_annual":{y:fit([r for r in sets[y] if r["key"] in common_local]) for y in YEARS}}
            panels.append(z);print(json.dumps({"naturalized_status_panel":z},allow_nan=False),flush=True)
    result={"version":"railstatus-naturalized-multiyear-v1","source_replay":replay,"panels":panels,
       "primary":"100 km, strict 0.3 degree GloNAF region guard, 50% land, 2023-25 pooled local region x block FE, active minus inactive coefficient on log1p documented minus log1p unlisted count; same-cell annual sensitivity",
       "limitations":["GloNAF documented versus unlisted is not alien versus native","2026 status applied to 2023-25 records, not opening/closure histories","pooled species/observer unions are not independent replicates","posthoc mechanism probe, not causal evidence of train transport","regional checklist coverage and abandoned tags uneven","approximate cluster intervals; residual spatial dependence"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
