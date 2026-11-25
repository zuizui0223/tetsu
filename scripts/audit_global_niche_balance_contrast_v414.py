#!/usr/bin/env python3
"""v414 paired block contrast of balanced versus original annual niche estimates."""
import argparse,json
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from analyze_global_niche_common_support_v400 import MODES,prepare,high,INTERACTION_FIELDS,COVERAGE_FIELDS,RESP
from analyze_global_niche_annual_change_v410 import YEARS,make_change,STAT,stat
from analyze_global_crossfit_climate_niche_v399 import group_indices,BOOT
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
SEED=41420261125
def contrast(rows0,rows1,fields,old0,old1,seed):
    n=len(rows0)
    if len(rows1)!=n or [r["key"] for r in rows0]!=[r["key"] for r in rows1]:raise ValueError("Different matched geography")
    g,b=group_indices(rows0,"block")
    for r,s in zip(rows0,rows1):
        if (r["continent"],r["block"])!=(s["continent"],s["block"]):raise ValueError("Different blocks")
    out={"n_cells":n,"n_blocks":len(b),"cells_by_continent":dict(Counter(r["continent"] for r in rows0))}
    if n<100 or len(b)<20 or n-len(g)-len(fields)<30:
        if any(m["status"]=="estimated" for m in (old0,old1)):raise ValueError("Prior support gate changed")
        return dict(out,status="insufficient common-year support")
    X=[np.asarray([[r[k] for k in fields+RESP] for r in rows],float) for rows in (rows0,rows1)]
    def est(w):
        a,b=(stat(x,g,w,fields) for x in X)
        if a is None or b is None:return None
        return np.stack((a,b,b-a))
    c=est(block_weights(rows0))
    if c is None:raise ValueError("Prior supported estimates lost")
    for k,old in enumerate((old0,old1)):
        if old["status"]!="estimated" or old["n_cells"]!=n or old["n_blocks"]!=len(b):raise ValueError("Prior support changed")
        for i,s in enumerate(STAT):
            for j,t in enumerate(RESP):
                if abs(c[k,i,j]-old["coefficients"][s][t])>1e-9:raise ValueError("Prior coefficient replay mismatch")
    bycont=defaultdict(list)
    for k in sorted(b):bycont[k[0]].append(k)
    rng=np.random.default_rng(seed);draws=[]
    for _ in range(BOOT):
        w=np.zeros(n)
        for bs in bycont.values():
            mult=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for k,m in zip(bs,mult):w[b[k]]=m/len(b[k])
        z=est(w)
        if z is not None:draws.append(z)
    names=("original","balanced","balanced_minus_original")
    out["status"]="estimated"
    out["coefficients"]={k:{s:dict(zip(RESP,map(float,c[h,i]))) for i,s in enumerate(STAT)} for h,k in enumerate(names)}
    out["intervals"]={k:{s:{t:{"valid":len(draws),
        "low":float(np.quantile(np.asarray(draws)[:,h,i,j],.025)) if len(draws)>=.95*BOOT else None,
        "high":float(np.quantile(np.asarray(draws)[:,h,i,j],.975)) if len(draws)>=.95*BOOT else None}
        for j,t in enumerate(RESP)} for i,s in enumerate(STAT)} for h,k in enumerate(names)}
    return out
def run(src,bal,niche,old):
    q=src["source_QA"]
    if src["version"]!="v405-original-annual-niche-coarse" or bal["version"]!="v412-original-balanced-annual-coarse" or niche["version"]!="v399-crossfit-niche-coarse-aggregate" or old["source_run"]!=45575269539 or q!=bal["source_QA"] or q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279":raise ValueError("Wrong original source")
    # JSON dict equality checks entire original source QA including annual counts.
    orig={tuple(r["key"]):r["annual"] for r in src["cells"]}
    data={tuple(r["key"]):r["balanced"] for r in bal["cells"]}
    if len(orig)!=1323 or set(data)!={k for k,r in orig.items() if all(r[y] is not None for y in YEARS)} or len(data)!=857:raise ValueError("Wrong original common support")
    land80={tuple(r["key"]) for p in niche["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    out=[]
    for p in old["panels"]:
        land=p["land"];h=p["support"]=="high_effort";mode=p["balance"]
        base=sorted(k for k in data if (land==.5 or k in land80) and (not h or all(high(orig[k][y][0]) for y in YEARS)))
        keys=[k for k in base if data[k][mode]["annual"] is not None and (not h or data[k][mode]["all_draws_high"])]
        if len(base)!=p["original_n"] or len(keys)!=p["n_retained"]:raise ValueError("Prior retained cells changed")
        panel={"land":land,"support":p["support"],"balance":mode,"n_retained":len(keys),"fits":[]}
        for j,ref in enumerate(MODES):
            f=p["fits"][j]
            if f["reference"]!=ref:raise ValueError("Prior ref mismatch")
            rr=[{y:[prepare(orig[k][y][j] if i==0 else data[k][mode]["annual"][y][j]) for k in keys] for y in YEARS} for i in range(2)]
            by={(m["variant"],m["model"]):m["change"] for m in f["models"]}
            for label,fields in (("competing",INTERACTION_FIELDS),("coverage_sensitivity",COVERAGE_FIELDS)):
                a,fs=make_change(rr[0],fields,"2025_minus_earlier_mean")
                b,gs=make_change(rr[1],fields,"2025_minus_earlier_mean")
                if fs!=gs:raise ValueError("Different model fields")
                z=contrast(a,b,fs,by["original_same_cells",label],by["balanced",label],SEED+int(land*10)+h)
                panel["fits"].append({"reference":ref,"model":label,**z})
        out.append(panel);print(json.dumps({"v414_panel":panel},allow_nan=False),flush=True)
    return {"version":"v414-original-v412-paired-observation-balance-contrast","source_QA":q,"panels":out,
      "limits":["same cell support and paired block draws, but estimates use differing observation subsets and covariates",
        "fixed record subsamples and reference draws; not full detection correction or uncertainty propagation",
        "difference is descriptive balance sensitivity, not causal observation bias decomposition",
        "posthoc comparison, no threshold relaxation; nonzero intervals not multiplicity adjusted",
        "failure to reject a difference is not proof of equivalence"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","balanced","niche","verified","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text());o=run(read(a.source),read(a.balanced),read(a.niche),read(a.verified));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
