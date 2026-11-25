#!/usr/bin/env python3
"""Status-specific cold-climate assembly on the verified v403 finite-pool
aggregates, not a new iNat scan. Direct interaction contrast, same cells
for both reference definitions and occupancy/family nulls.
"""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_world_railstatus_population_local import (
    BASE,ROAD,POP,LOG,CATS,inject_status,block_weights,read)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import RICH_BUILT
W=100000
MODES=("exclude_3x3_blocks","exclude_continent")
LEVELS=("occupancy","family_occupancy","genus_occupancy")
TAGS=("T_variance","width_C","warm_offset_C")
OUT=tuple(t+"_"+k for k in TAGS for t in ("observed","expected","excess"))
RAIL=tuple(LOG[c] for c in CATS)
ACTIVE=LOG["active_normal"]
INACTIVE=LOG["inactive_normal"]
BASE_FIELDS=BASE+ROAD+(POP,"bio4_raw","log1p_bio12","bio15_raw","trait_coverage")
MODIFIERS=("log_road_density",RICH_BUILT,POP,"tree_consensus_fullcell",
           "land_fraction","log_observations","log_observers","log1p_bio12",
           "bio4_raw","bio15_raw","trait_coverage")
FIELDS=BASE_FIELDS+RAIL+tuple(k+"_cold" for k in RAIL+MODIFIERS)
TARGETS={
    "active_cold":{ACTIVE+"_cold":1.},
    "inactive_cold":{INACTIVE+"_cold":1.},
    "active_minus_inactive_cold":{ACTIVE+"_cold":1.,INACTIVE+"_cold":-1.},
    "active_at_5C":{ACTIVE:1.,ACTIVE+"_cold":1.},
    "inactive_at_5C":{INACTIVE:1.,INACTIVE+"_cold":1.},
    "active_minus_inactive_at_5C":{ACTIVE:1.,ACTIVE+"_cold":1.,INACTIVE:-1.,INACTIVE+"_cold":-1.},
    "active_at_15C":{ACTIVE:1.},"inactive_at_15C":{INACTIVE:1.},
}
EXPECTED={(.5,False):1323,(.5,True):742,(.8,False):1116,(.8,True):618}

def prepare(r):
    z=dict(r,key=tuple(r["key"]),block=tuple(r["block"]))
    cold=(15.-z["bio1"])/10.
    for k in RAIL+MODIFIERS:z[k+"_cold"]=z[k]*cold
    for k in TAGS:
        if abs(z["observed_"+k]-z["expected_"+k]-z["excess_"+k])>1e-9:
            raise ValueError("Source null decomposition not additive")
    return z

def fit(rows,fields=FIELDS):
    blocks=defaultdict(list)
    for i,r in enumerate(rows):blocks[(r["continent"],r["block"])].append(i)
    q={"n_cells":len(rows),"n_blocks":len(blocks),
       "cells_by_continent":dict(Counter(r["continent"] for r in rows)),
       "n_both_positive":sum(r.get("active_normal_km",0)>0 and r.get("inactive_normal_km",0)>0 for r in rows)}
    if len(rows)<100 or len(blocks)<20 or len(rows)-len(blocks)-len(fields)<30:
        return dict(q,status="insufficient support")
    q["temperature_quantiles_C"]=list(map(float,np.quantile([r["bio1"] for r in rows],[0,.05,.5,.95,1])))
    w=block_weights(rows)
    X=np.array([[r[f] for f in fields+OUT] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite source")
    for ix in blocks.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None]
    A=X[:,:len(fields)];Y=X[:,len(fields):];scale=np.sqrt((A*A).sum(axis=0))
    keep=scale>1e-10
    names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE,ACTIVE+"_cold",INACTIVE+"_cold"}<=set(names):
        return dict(q,status="insufficient status variation")
    Z=A[:,keep]/scale[keep]
    singular=np.linalg.svd(Z,compute_uv=False)
    if singular[-1]<1e-8*singular[0]:
        return dict(q,status="rank deficient",condition_number=float(singular[0]/max(singular[-1],1e-20)))
    inv=np.linalg.inv(Z.T@Z);b=inv@(Z.T@Y);e=Y-Z@b
    coef=b/scale[keep,None]
    scores=np.array([(inv@(Z[ix].T@e[ix]))/scale[keep,None] for ix in blocks.values()])
    effects={}
    for label,terms in TARGETS.items():
        c=np.array([terms.get(k,0.) for k in names])
        est=c@coef;inf=np.einsum("i,bio->bo",c,scores)
        se=np.sqrt(len(blocks)/(len(blocks)-1)*(inf**2).sum(axis=0))
        for j in range(0,len(OUT),3):
            if abs(est[j]-est[j+1]-est[j+2])>1e-8:raise ValueError("Fitted null decomposition failed")
        effects[label]={o:{"estimate":float(v),"interval95_approx":[float(v-1.96*z),float(v+1.96*z)]}
                        for o,v,z in zip(OUT,est,se)}
    return dict(q,status="estimated",condition_number=float(singular[0]/singular[-1]),
                dropped_constant_fields=[f for f in fields if f not in names],effects=effects)

def summary(rows):
    return {"all_status":fit(rows),"both_positive":fit([r for r in rows if r["active_normal_km"]>0 and r["inactive_normal_km"]>0]),
            "leave_one_continent_out":{c:fit([r for r in rows if r["continent"]!=c]) for c in sorted({r["continent"] for r in rows})}}

def run(src,niche,prior,status):
    q=src["source_QA"]
    if src["version"]!="v403-original-finite-pool-null-coarse" or niche["version"]!="v399-crossfit-niche-coarse-aggregate" or prior["version"]!="v407-original-common-cell-taxonomic-null-decomposition":
        raise ValueError("Wrong original coarse artifacts")
    if q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:
        raise ValueError("Wrong original observation replay")
    if prior["source_QA"]!=q or niche["source_QA"]["source_observation_SHA256"]!=q["source_observation_SHA256"]:
        raise ValueError("Original niche source mismatch")
    if not status["version"].startswith("railstatus-v1") or status["source_Overture_release"]!="2026-11-18.0" or abs(status["total_km"]["active_normal"]-2237043.481)>.01:
        raise ValueError("Wrong original status source")
    raw=src["rows"]
    if len(raw)!=1323*2*3:raise ValueError("Original v403 count changed")
    raw=[dict(r,key=tuple(r["key"]),block=tuple(r["block"])) for r in raw]
    rows,worst=inject_status(raw,status,W)
    rows=[prepare(r) for r in rows]
    land80={tuple(r["key"]) for p in niche["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    previous={(p["land"],p["support"],p["reference"]):p for p in prior["panels"]}
    panels=[]
    for land in (.5,.8):
        for high in (False,True):
            support="high_effort" if high else "all"
            maps={}
            for mode in MODES:
                for level in LEVELS:
                    base=[r for r in rows if r["reference_exclusion"]==mode and r["null_level"]==level and
                          (land==.5 or r["key"] in land80) and
                          (not high or (r["observations"]>=200 and r["observers"]>=20))]
                    if len(base)!=EXPECTED[land,high]:raise ValueError("v403 support denominator mismatch")
                    eligible=[r for r in base if r["pool_cells"]>=3 and r["movable_fraction"]>=.20]
                    maps[mode,level]={r["key"]:r for r in eligible}
                    if len(maps[mode,level])!=previous[land,support,mode]["original_support"][level]:
                        raise ValueError("v407 eligible support mismatch")
            for levels in (LEVELS[:2],LEVELS):
                keys=set.intersection(*(set(maps[m,l]) for m in MODES for l in levels))
                if not keys:raise ValueError("No common source niche cells")
                for k in keys:
                    a=maps[MODES[0],levels[0]][k]
                    for m in MODES:
                        for l in levels:
                            b=maps[m,l][k]
                            if any(abs(a[f]-b[f])>1e-10 for f in FIELDS) or (a["continent"],a["block"])!=(b["continent"],b["block"]):
                                raise ValueError("Covariates differ on identical source geography")
                z={"land":land,"support":support,"levels":levels,"n_common":len(keys),"models":[]}
                for mode in MODES:
                    # Replay common outcome equality for the observed assemblage.
                    for k in keys:
                        a=maps[mode,levels[0]][k]
                        for level in levels[1:]:
                            if any(abs(a["observed_"+t]-maps[mode,level][k]["observed_"+t])>1e-10 for t in TAGS):
                                raise ValueError("Observed niche changed by null")
                    for level in levels:
                        r=[maps[mode,level][k] for k in sorted(keys)]
                        z["models"].append({"reference":mode,"null":level,**summary(r)})
                panels.append(z)
                print(json.dumps({"status_niche_panel":z},allow_nan=False),flush=True)
    return {"version":"railstatus-niche-pool-v1","source_QA":q,"active_km_max_replay_error":worst,"panels":panels,
        "primary":"100km, land>=50%, all effort, occupancy+family common cells across both references, occupancy null excess among-species temperature variance, active minus inactive cold interaction; family null and alternate reference sensitivity",
        "limits":["posthoc mechanism probe based on prior climate results, not independent replication","2018-22 realized occupancy macroclimate, not tolerance or measured ballast temperature","regional null pool includes focal cells, not independent available flora","mapped active candidate status not observed operations","landscape records, not trackside assemblages","approximate cluster intervals; correlated blocks and multiple comparisons","conditional on widespread taxa and shared reference support"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("pool","niche","prior","status","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();r=run(*map(read,(a.pool,a.niche,a.prior,a.status)))
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,allow_nan=False))
if __name__=="__main__":main()
