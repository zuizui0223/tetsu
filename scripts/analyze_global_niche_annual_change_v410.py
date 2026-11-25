#!/usr/bin/env python3
"""v410: fixed-geography change in observed niche metrics, not causal rail effect."""
import argparse,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_crossfit_climate_niche_v399 import group_indices,slopes,BOOT,TARGET
from analyze_global_niche_common_support_v400 import COLD,RESP,MODES,COVERAGE_FIELDS,INTERACTION_FIELDS,prepare,high
from analyze_global_niche_temporal_v405 import STAT,stat
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
YEARS=("2023","2024","2025")
SEED=41020261125
DYNAMIC=("log_observations","log_observers","log_species","trait_coverage")
def is_dynamic(k):
    return k in DYNAMIC or k in tuple(x+"_cold" for x in DYNAMIC) or k=="coverage_cold"
def make_change(rows,fields,contrast):
    """Baseline (2023–24 mean or 2023) and annual change covariates.
    Outcomes and effort/coverage changes computed inside the same cell.
    """
    out=[]
    for rr in zip(*(rows[y] for y in YEARS)):
        a,b,c=rr
        if len({r["key"] for r in rr})!=1:raise ValueError("Cell mismatch")
        prev=(a,b) if contrast=="2025_minus_earlier_mean" else (a,)
        current=c if contrast=="2025_minus_earlier_mean" else b
        mean=lambda k:float(np.mean([r[k] for r in prev]))
        z=dict(current)
        for k in fields:
            if is_dynamic(k):
                z[k]=mean(k)
                z["delta_"+k]=current[k]-z[k]
            elif max(r[k] for r in rr)-min(r[k] for r in rr)>1e-10:
                raise ValueError("Stable covariate changed: "+k)
        for k in RESP:z[k]=current[k]-mean(k)
        out.append(z)
    return out,fields+tuple("delta_"+k for k in fields if is_dynamic(k))
def fit(rows,fields,seed):
    g,b=group_indices(rows,"block")
    n=len(rows);o={"n_cells":n,"n_blocks":len(b),"cells_by_continent":dict(Counter(r["continent"] for r in rows))}
    if n<100 or len(b)<20 or n-len(g)-len(fields)<30:return dict(o,status="insufficient common-year support")
    X=np.asarray([[r[k] for k in fields+RESP] for r in rows],float);w=block_weights(rows)
    a=stat(X,g,w,fields)
    if a is None:return dict(o,status="insufficient residual exposure")
    bc=defaultdict(list)
    for k in sorted(b):bc[k[0]].append(k)
    rng=np.random.default_rng(seed);ds=[]
    for _ in range(BOOT):
        bw=np.zeros(n)
        for bs in bc.values():
            m=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for k,v in zip(bs,m):bw[b[k]]=v/len(b[k])
        s=stat(X,g,bw,fields)
        if s is not None:ds.append(s)
    o["status"]="estimated"
    o["coefficients"]={s:dict(zip(RESP,map(float,a[i]))) for i,s in enumerate(STAT)}
    o["intervals"]={s:{t:{"valid":len(ds),"low":float(np.quantile(np.asarray(ds)[:,i,j],.025)) if len(ds)>=.95*BOOT else None,
                           "high":float(np.quantile(np.asarray(ds)[:,i,j],.975)) if len(ds)>=.95*BOOT else None}
                         for j,t in enumerate(RESP)} for i,s in enumerate(STAT)}
    return o
def replay(rows,fields,expected):
    g,_=group_indices(rows["2025"],"block");w=block_weights(rows["2025"])
    for y in YEARS:
        X=np.asarray([[r[k] for k in fields+RESP] for r in rows[y]],float)
        c=stat(X,g,w,fields)
        if c is None:raise ValueError("Prior fit not reproduced")
        for i,s in enumerate(STAT):
            for j,t in enumerate(RESP):
                if abs(c[i,j]-expected[y][s][t])>1e-9:
                    raise ValueError("Exact v405 replay mismatch")
def run(src,niche,old):
    q=src["source_QA"]
    if src["version"]!="v405-original-annual-niche-coarse" or niche["version"]!="v399-crossfit-niche-coarse-aggregate" or old["source_run"]!=45568291544 or q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original source")
    if len(src["cells"])!=1323:raise ValueError("Original geography mismatch")
    land80={tuple(r["key"]) for p in niche["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    data={tuple(r["key"]):r["annual"] for r in src["cells"]}
    out=[]
    for p in old["panels"]:
        land=p["land"];h=p["support"]=="high_effort"
        keys=sorted(k for k,v in data.items() if (land==.5 or k in land80) and all(v[y] is not None and (not h or high(v[y][0])) for y in YEARS))
        if len(keys)!=p["common_three_year_cells"]:raise ValueError("Exact v405 sample count mismatch")
        for j,mode in enumerate(MODES):
            rr={y:[prepare(data[k][y][j]) for k in keys] for y in YEARS}
            f=p["fits"][j]
            if f["reference"]!=mode:raise ValueError("Reference order mismatch")
            z={"land":land,"support":p["support"],"reference":mode,"n":len(keys),"fits":[]}
            for label,fields in (("competing",INTERACTION_FIELDS),("coverage_sensitivity",COVERAGE_FIELDS)):
                replay(rr,fields,f["models"][label]["coefficients"])
                for contrast in ("2025_minus_earlier_mean","2024_minus_2023"):
                    a,fs=make_change(rr,fields,contrast)
                    m=fit(a,fs,SEED+int(land*10)+h)
                    non_eu=[r for r in a if r["continent"]!="Europe"]
                    eu=[r for r in a if r["continent"]=="Europe"]
                    z["fits"].append({"model":label,"contrast":contrast,"all":m,
                      "without_Europe":fit(non_eu,fs,SEED+int(land*10)+h),
                      "Europe_only":fit(eu,fs,SEED+int(land*10)+h)})
            out.append(z);print(json.dumps({"v410_panel":z},allow_nan=False),flush=True)
    return {"version":"v410-original-annual-within-cell-niche-change","source_QA":q,"panels":out,
      "limits":["fixed mapped rail and fixed reference climate, not exposure or actual climate change",
      "observed composition change can reflect observer/season/taxon selection; not proof of ecological change",
      "2025-selected geography, repeated observations not independent ecosystems",
      "baseline/change effort and coverage adjustment is sensitivity not detection correction",
      "2024 minus 2023 is a descriptive preceding-year comparison, not independent placebo",
      "posthoc same-data analysis, multiple correlated outcomes, no causal inference"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","verified","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.source),read(a.niche),read(a.verified));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
