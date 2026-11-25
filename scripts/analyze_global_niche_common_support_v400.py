#!/usr/bin/env python3
"""v400: same-cell niche-reference support and competing cold interactions."""
import argparse,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_global_crossfit_climate_niche_v399 import (
    FIELDS,TARGET,OUTCOMES,MIN_KNOWN,MIN_COVERAGE,BOOT,group_indices,slopes,regression)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import RICH_BUILT
SEED=40020261125
COLD="rail_cold"
MODIFIERS=("log_road_density",RICH_BUILT,"log_population_density","tree_consensus_fullcell",
           "land_fraction","log_observations","log_observers","log1p_bio12","bio4_raw","bio15_raw")
INTERACTION_FIELDS=FIELDS+(COLD,)+tuple(k+"_cold" for k in MODIFIERS)
COVERAGE_FIELDS=INTERACTION_FIELDS+("trait_coverage","coverage_cold")
RESP=("temp_width_C","among_species_temperature_sd_C","warm_offset_C","cold_tail_offset_C",
      "below_reference_q05_fraction","dry_offset_log_precip","precip_width_log")
MODES=("exclude_3x3_blocks","exclude_continent")
EXPECTED={
 ("exclude_3x3_blocks",.5):[(1954,-.04675956367826274,-.02753656126494103),(1193,-.03212981165487984,-.036984554324859316)],
 ("exclude_3x3_blocks",.8):[(1672,-.048992546896861125,-.036529581140456166),(1010,-.04362103851556237,-.0570844660670073)],
 ("exclude_continent",.5):[(1333,-.025646015648890377,-.028186779356425117),(742,-.01950400746352303,-.018420702472028663)],
 ("exclude_continent",.8):[(1125,-.023207725848818455,-.029694137431485548),(618,-.011515494070069276,-.0470235887475799)]}
def accepted(r):
    return r["known_niche_species"]>=MIN_KNOWN and r["trait_coverage"]>=MIN_COVERAGE
def high(r):return r["observations"]>=200 and r["observers"]>=20
def prepare(r):
    c=(15-r["bio1"])/10.
    return dict(r,key=tuple(r["key"]),block=tuple(r["block"]),rail_cold=r[TARGET]*c,
                coverage_cold=r["trait_coverage"]*c,**{k+"_cold":r[k]*c for k in MODIFIERS})
def summarize(v,outcomes):
    return {k:{"valid":len(v),"low":float(np.quantile(np.asarray(v)[:,j],.025)) if len(v)>=.95*BOOT else None,
            "high":float(np.quantile(np.asarray(v)[:,j],.975)) if len(v)>=.95*BOOT else None} for j,k in enumerate(outcomes)}
def model(rows,fields,seed):
    g,b=group_indices(rows,"block")
    o={"n_cells":len(rows),"n_blocks":len(b),"cells_by_continent":dict(Counter(r["continent"] for r in rows))}
    if len(rows)<100 or len(b)<20 or len(rows)-len(g)-len(fields)<30:
        return dict(o,status="insufficient support")
    X=np.asarray([[r[k] for k in fields+RESP] for r in rows],float);w=block_weights(rows)
    targets=(TARGET,COLD)
    def est(X,g,w):
        a=slopes(X,g,w,fields,targets)
        if a is None:return None
        z={k:np.asarray(v) for k,v in a.items()}
        for t in (5,15,25):z["rail_slope_at_"+str(t)+"C"]=z[TARGET]+((15-t)/10.)*z[COLD]
        return z
    real=est(X,g,w)
    if real is None:return dict(o,status="insufficient residual variation")
    o.update(status="estimated",coefficients={k:dict(zip(RESP,map(float,v))) for k,v in real.items()})
    bc=defaultdict(list)
    for k in sorted(b):bc[k[0]].append(k)
    rng=np.random.default_rng(seed);draws={k:[] for k in real}
    for _ in range(BOOT):
        bw=np.zeros(len(rows))
        for bs in bc.values():
            m=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for k,n in zip(bs,m):bw[b[k]]=n/len(b[k])
        z=est(X,g,bw)
        if z is not None:
            for k,v in z.items():draws[k].append(v)
    o["block_percentile_95_descriptive"]={k:summarize(v,RESP) for k,v in draws.items()}
    for label,units,keyfunc in (("leave_one_block_out",sorted(b),lambda r:(r["continent"],r["block"])),
                                ("leave_one_continent_out",sorted(bc),lambda r:r["continent"])):
        v={k:[] for k in real}
        for u in units:
            sub=[r for r in rows if keyfunc(r)!=u];gg,bb=group_indices(sub,"block")
            if len(sub)<100 or len(bb)<20 or len(sub)-len(gg)-len(fields)<30:continue
            q=est(np.asarray([[r[k] for k in fields+RESP] for r in sub],float),gg,block_weights(sub))
            if q is not None:
                for k,a in q.items():v[k].append(a)
        o[label]={"attempted":len(units),"ranges":{k:{t:{"valid":len(a),
            "min":float(np.min(np.asarray(a)[:,j])) if a else None,"max":float(np.max(np.asarray(a)[:,j])) if a else None,
            "negative":int(np.sum(np.asarray(a)[:,j]<0)) if a else 0} for j,t in enumerate(RESP)} for k,a in v.items()}}
    o["observed_local_temperature_quantiles_C"]=list(map(float,np.quantile([r["bio1"] for r in rows],[0,.05,.5,.95,1])))
    return o
def run(src):
    if src["version"]!="v399-crossfit-niche-coarse-aggregate":raise ValueError("Wrong aggregate version")
    q=src["source_QA"]
    if q["all_source_rows"]!=288721938 or q["eligible_by_year"]["2025"]!=7827936 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original replay")
    panels={(p["reference_exclusion"],p["min_land_fraction"]):[prepare(r) for r in p["rows"]] for p in src["panels"]}
    if set(panels)!=set(EXPECTED):raise ValueError("Missing original panels")
    out=[]
    for gate in (.5,.8):
        for i,effort in enumerate(("all","high_effort")):
            rr={k:[r for r in panels[k,gate] if accepted(r) and (i==0 or high(r))] for k in MODES}
            maps={k:{r["key"]:r for r in rr[k]} for k in MODES}
            common=set(maps[MODES[0]])&set(maps[MODES[1]])
            if any(len(maps[k])!=len(rr[k]) for k in MODES):raise ValueError("Duplicate cells")
            for k in common:
                a,b=maps[MODES[0]][k],maps[MODES[1]][k]
                if any(a[f]!=b[f] for f in FIELDS+("continent","block","species","observations","observers")):raise ValueError("Same geography covariates differ")
            p={"min_land_fraction":gate,"support":effort,"common_n":len(common),"fits":[]}
            for mode in MODES:
                r=rr[mode];n,d,h=EXPECTED[mode,gate][i]
                prior=regression(r,"block",FIELDS,(TARGET,),OUTCOMES,False)
                inter=regression(r,"block",FIELDS+(COLD,),(TARGET,COLD),OUTCOMES,False)
                if len(r)!=n or abs(prior["coef"][TARGET]["temp_width_C"]-d)>1e-9 or abs(inter["coef"][COLD]["among_species_temperature_sd_C"]-h)>1e-9:raise ValueError("Exact v399 coefficient not reproduced")
                r=[x for x in r if x["key"] in common]
                seed=SEED+int(gate*10)+i
                m={"reference_exclusion":mode,"original_n":n,"same_cells":len(r),
                    "same_cell_main":regression(r,"block",FIELDS,(TARGET,),RESP,True),
                    "same_cell_original_interaction":model(r,FIELDS+(COLD,),seed),
                    "same_cell_competing_interactions":model(r,INTERACTION_FIELDS,seed),
                    "same_cell_competing_plus_coverage":model(r,COVERAGE_FIELDS,seed)}
                p["fits"].append(m)
            out.append(p)
            summary=[]
            for m in p["fits"]:
                v={"reference_exclusion":m["reference_exclusion"],"original_n":m["original_n"],"same_cells":m["same_cells"]}
                for k in ("same_cell_main","same_cell_original_interaction","same_cell_competing_interactions","same_cell_competing_plus_coverage"):
                    a=m[k];v[k]={"status":a["status"],"coef":a.get("coef",a.get("coefficients")),"ci":a.get("block_percentile_95_descriptive"),
                       "block_deletion":a.get("leave_one_block_out"),"continent_deletion":a.get("leave_one_continent_out")}
                summary.append(v)
            print(json.dumps({"v400_panel":dict(p,fits=summary)},allow_nan=False),flush=True)
    return {"version":"v400-original-common-niche-support-and-competing-climate-interactions",
       "panels":out,"source_QA":q,"modifiers":MODIFIERS,"bootstrap_draws":BOOT,
       "limits":["same evaluated cells but reference-known species sets may differ",
                 "coverage adjustment is sensitivity, may condition on selection; does not recover unknown taxa",
                 "competitor interactions adjust measured context, not causality",
                 "marginal slopes at fixed 5/15/25C may extrapolate; report observed support",
                 "same source, exploratory follow-up to v399, not independent confirmation",
                 "realized macroclimatic occupancy, not physiological tolerance or railway-side plants"]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();o=run(json.loads(a.source.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
