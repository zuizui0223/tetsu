#!/usr/bin/env python3
"""v402: geographic held-out prediction of niche width and convergence."""
import argparse,json,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_crossfit_climate_niche_v399 import FIELDS,TARGET,OUTCOMES,regression
from analyze_global_niche_common_support_v400 import (
    prepare,accepted,high,MODES,EXPECTED,INTERACTION_FIELDS,COVERAGE_FIELDS,COLD)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    fit_predict,metrics,CONTINENTS,RIDGE)
SEED="tetsu-v402-geographic-holdout-20261125"
TARGETS=OUTCOMES[:-1]
MIN_TRAIN=300
MIN_TEST=100
MIN_TEST_BLOCKS=10
MIN_CONTINENTS=4
N_FOLDS=5
BASE=tuple(k for k in INTERACTION_FIELDS if k not in (TARGET,COLD))
MODELS={"environment":BASE,"rail_main":BASE+(TARGET,),
        "rail_main_and_cold":BASE+(TARGET,COLD)}
def block(r):return (r["continent"],r["block"])
def fold(k):
    c,(x,y)=k
    h=hashlib.sha256(f"{SEED}:{c}:{x}:{y}".encode()).digest()
    return int.from_bytes(h[:8],"big")%N_FOLDS
def separated(a,b):
    # Remove training blocks within the held-out 3x3 block neighbourhood;
    # even when continents differ near a boundary.
    return max(abs(a[1][0]-b[1][0]),abs(a[1][1]-b[1][1]))>1
def partitions(rows,mode):
    if mode=="continent":
        for c in CONTINENTS:
            te=[r for r in rows if r["continent"]==c]
            tr=[r for r in rows if r["continent"]!=c]
            yield c,tr,te
    elif mode=="buffered_block":
        for f in range(N_FOLDS):
            te=[r for r in rows if fold(block(r))==f]
            tb={block(r) for r in te}
            # Build a spatial buffer before fitting. No response-driven splits.
            excluded={(x+i,y+j) for _,(x,y) in tb for i in (-1,0,1) for j in (-1,0,1)}
            tr=[r for r in rows if r["block"] not in excluded]
            yield str(f),tr,te
    else:raise ValueError("unknown holdout")
def evaluate(rows,mode,models):
    report=[];per_target=defaultdict(list)
    for label,tr,te in partitions(rows,mode):
        nblocks=len({block(r) for r in te});trainblocks=len({block(r) for r in tr})
        if {block(r) for r in tr}&{block(r) for r in te}:raise ValueError("Block overlap")
        if mode=="buffered_block" and not all(separated(a,b) for a in {block(r) for r in tr} for b in {block(r) for r in te}):raise ValueError("Buffer overlap")
        rec={"fold":label,"n_train":len(tr),"n_test":len(te),"train_blocks":trainblocks,"test_blocks":nblocks}
        if len(tr)<max(MIN_TRAIN,4*max(map(len,models.values()))) or len(te)<MIN_TEST or nblocks<MIN_TEST_BLOCKS:
            report.append(dict(rec,status="insufficient heldout/training support"));continue
        rr={}
        for t in TARGETS:
            mse={};mae={}
            for k,fields in models.items():
                pred=fit_predict(tr,te,fields,t,t=="below_reference_q05_fraction")
                mse[k],mae[k]=metrics(te,pred,t)
            d={"rail_main_delta_MSE":mse["environment"]-mse["rail_main"],
               "rail_joint_delta_MSE":mse["environment"]-mse["rail_main_and_cold"],
               "interaction_delta_MSE":mse["rail_main"]-mse["rail_main_and_cold"]}
            rr[t]={"mse":mse,"mae":mae,**d}
            per_target[t].append(d)
        report.append(dict(rec,status="estimated",outcomes=rr))
    n=sum(r["status"]=="estimated" for r in report)
    summary={}
    for t,v in per_target.items():
        summary[t]={k:{"mean":float(np.mean([x[k] for x in v])),
                      "positive_folds":sum(x[k]>0 for x in v)} for k in v[0]}
    return {"n_valid_folds":n,"min_valid_folds":MIN_CONTINENTS if mode=="continent" else N_FOLDS,
        "support_pass":n>=(MIN_CONTINENTS if mode=="continent" else N_FOLDS),
        "folds":report,"fold_equal_summary":summary}
def run(src):
    q=src["source_QA"]
    if src["version"]!="v399-crossfit-niche-coarse-aggregate" or q["all_source_rows"]!=288721938 or q["eligible_by_year"]["2025"]!=7827936 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original source")
    panels={(p["reference_exclusion"],p["min_land_fraction"]):[prepare(r) for r in p["rows"]] for p in src["panels"]}
    if set(panels)!=set(EXPECTED):raise ValueError("Missing original panel")
    out=[]
    for gate in (.5,.8):
        for i,effort in enumerate(("all","high_effort")):
            rr={m:[r for r in panels[m,gate] if accepted(r) and (not i or high(r))] for m in MODES}
            common=set.intersection(*[{r["key"] for r in a} for a in rr.values()])
            if len(common)!={(.5,0):1326,(.5,1):742,(.8,0):1119,(.8,1):618}[gate,i]:raise ValueError("Exact v400 support mismatch")
            for mode in MODES:
                n,d,h=EXPECTED[mode,gate][i];r=rr[mode]
                old=regression(r,"block",FIELDS,(TARGET,),OUTCOMES,False)
                if len(r)!=n or abs(old["coef"][TARGET]["temp_width_C"]-d)>1e-9:raise ValueError("Exact v399 source coefficient mismatch")
                r=[x for x in r if x["key"] in common]
                p={"land":gate,"support":effort,"reference_exclusion":mode,"n":len(r),
                   "cells_by_continent":dict(Counter(x["continent"] for x in r)),"fits":[]}
                for label,fields in (("competing",MODELS),("coverage_sensitivity",{k:v+("trait_coverage","coverage_cold") for k,v in MODELS.items()})):
                    for holdout in ("continent","buffered_block"):
                        p["fits"].append({"model":label,"holdout":holdout,**evaluate(r,holdout,fields)})
                out.append(p)
                print(json.dumps({"v402_panel":dict(p,fits=[{k:v for k,v in f.items() if k!="folds"}|{"fold_support":[{k:v for k,v in x.items() if k!="outcomes"} for x in f["folds"]]} for f in p["fits"]])},allow_nan=False),flush=True)
    return {"version":"v402-original-niche-geographic-heldout-prediction","source_QA":q,"panels":out,
            "ridge_lambda":RIDGE,"models":MODELS,"seed":SEED,
            "limits":["held-out regression geography is not independent niche reference reconstruction; species may recur across folds",
             "strict 3x3 block buffers can remove most training data; insufficient folds remain unestimated",
             "out-of-continent prediction limited by geographically uneven taxon coverage",
             "predictive gain is unsigned and not causation, no local microclimate inference",
             "same observational snapshot, no model tuning on held-out outcomes"]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();o=run(json.loads(a.source.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
