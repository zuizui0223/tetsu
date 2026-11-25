#!/usr/bin/env python3
"""v408: support-selection audit of the three-level v407 decomposition."""
import argparse,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_global_niche_pool_null_v403 import LEVELS,MIN_MOVABLE,EXPECTED
from analyze_global_niche_common_support_v400 import MODES,prepare,high,COVERAGE_FIELDS,COLD
from analyze_global_crossfit_climate_niche_v399 import regression,TARGET
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
FEATURES=("bio1","abs_lat","log_road_density","log1p_1000_times_original_2020_builtup_fraction",
          "log_population_density","tree_consensus_fullcell","land_fraction","log_observations",
          "log_observers","log1p_bio12","bio4_raw","bio15_raw","log_rail_km","trait_coverage",
          "known_niche_species","pool_cells","pool_taxa")
TARGETS=("observed_T_variance","observed_width_C","observed_warm_offset_C")
def support(r):return r["pool_cells"]>=3 and r["movable_fraction"]>=MIN_MOVABLE
def describe(rows):
    if not rows:return {"n":0}
    return {"n":len(rows),"n_blocks":len({(r["continent"],r["block"]) for r in rows}),
          "continents":dict(Counter(r["continent"] for r in rows)),
          "mean":{k:float(np.mean([r[k] for r in rows])) for k in FEATURES},
          "quantiles":{k:list(map(float,np.quantile([r[k] for r in rows],[.05,.25,.5,.75,.95]))) for k in FEATURES}}
def compare(a,b):
    if not a or not b:return {}
    X=np.array([[r[k] for k in FEATURES] for r in a]);Y=np.array([[r[k] for k in FEATURES] for r in b])
    s=np.sqrt((X.var(axis=0)+Y.var(axis=0))/2)
    return {k:float(d/v) if v>1e-12 else None for k,d,v in zip(FEATURES,X.mean(axis=0)-Y.mean(axis=0),s)}
def run(src,orig,verified):
    q=src["source_QA"]
    if src["version"]!="v403-original-finite-pool-null-coarse" or orig["version"]!="v399-crossfit-niche-coarse-aggregate" or q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong source")
    if verified["source_run_id"]!=45570583651:raise ValueError("Wrong v407 source")
    check={(p["land"],p["support"],p["reference"],len(p["levels"])):p for p in verified["panels"]}
    if len(check)!=16:raise ValueError("Missing prior panels")
    land80={tuple(r["key"]) for p in orig["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    rows=[prepare(r) for r in src["rows"]]
    if len(rows)!=1323*2*3:raise ValueError("Wrong source count")
    out=[]
    for land in (.5,.8):
        for h in (False,True):
            effort="high_effort" if h else "all"
            for mode in MODES:
                maps={}
                for lev in LEVELS:
                    r=[x for x in rows if x["reference_exclusion"]==mode and x["null_level"]==lev and (land==.5 or x["key"] in land80) and (not h or high(x))]
                    if len(r)!=EXPECTED[land,h]:raise ValueError("Wrong original panel")
                    maps[lev]={x["key"]:x for x in r if support(x)}
                base=set(maps[LEVELS[0]])&set(maps[LEVELS[1]])
                incl=base&set(maps[LEVELS[2]])
                keep=[maps[LEVELS[0]][k] for k in sorted(incl)]
                drop=[maps[LEVELS[0]][k] for k in sorted(base-incl)]
                full=[maps[LEVELS[0]][k] for k in sorted(base)]
                p={"land":land,"support":effort,"reference":mode,"base_n":len(full),
                   "included":describe(keep),"excluded":describe(drop),
                   "retention_by_continent":{},"included_minus_excluded_SMD":compare(keep,drop)}
                b=Counter(r["continent"] for r in full);a=Counter(r["continent"] for r in keep)
                p["retention_by_continent"]={k:{"base":n,"included":a[k],"fraction":a[k]/n} for k,n in b.items()}
                p["models"]={}
                for name,r,old in (("base",full,2),("genus_eligible",keep,3),("genus_ineligible",drop,None)):
                    m=regression(r,"block",COVERAGE_FIELDS,(TARGET,COLD),TARGETS,False)
                    if old is not None:
                        c=check[land,effort,mode,old]
                        if len(r)!=c["n"] or abs(m["coef"][COLD]["observed_T_variance"]-c["coefficients"]["observed_T_variance"])>1e-9:
                            raise ValueError("Exact v407 observed slope replay failed")
                    p["models"][name]=m
                out.append(p);print(json.dumps({"v408_panel":p},allow_nan=False),flush=True)
    return {"version":"v408-original-genus-support-selection-audit","panels":out,"source_QA":q,
       "limits":["descriptive support selection, not correction for missing taxa or causal selection",
          "SMD uses unweighted cell moments, not independent spatial sample statistics",
          "observed slope subsets are not independent replications; no new outcome tests",
          "conditional label support does not establish actual botanical availability"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","verified","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.source),read(a.niche),read(a.verified));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
