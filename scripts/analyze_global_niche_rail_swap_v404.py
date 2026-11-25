#!/usr/bin/env python3
"""v404: local conditional railway exposure swaps on shared, rarefied niche cells."""
import argparse,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_crossfit_climate_niche_v399 import TARGET,group_indices,slopes
from analyze_global_niche_common_support_v400 import MODES,prepare,high,COVERAGE_FIELDS,COLD,RESP
from analyze_global_naturalized_joint_population_local_v394 import MATCH
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
SEED=40420261125
DRAWS=199
CALIPER=.75
FEATURES=MATCH+("bio4_raw","log1p_bio12","bio15_raw")
EXPECTED={
 (.5,False):(1323,(-.021188242841640774,-.02977653502712443),(-.021763360918019763,-.016243670898003725)),
 (.5,True):(742,(-.02696033802708125,-.034587764782342774),(-.02152474089170227,-.05629486458290198)),
 (.8,False):(1116,(-.026261373863920728,-.03528462063037689),(-.02281961035138914,-.026779815305515806)),
 (.8,True):(618,(-.050053421950253134,-.03563783216494552),(-.041904818601527144,-.05300146675424044))}
def bkey(r):return (r["continent"],r["block"])
def pairs(rows):
    X=np.asarray([[r[k] for k in FEATURES] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite matching features")
    C=defaultdict(list);B=defaultdict(list)
    for i,r in enumerate(rows):C[r["continent"]].append(i);B[bkey(r)].append(i)
    Z=X.copy()
    for ix in C.values():
        s=X[ix].std(axis=0);s[s<1e-12]=1.;Z[ix]=(X[ix]-X[ix].mean(axis=0))/s
    out=[]
    for ix in B.values():
        cand=[]
        for j,i in enumerate(ix):
            for k in ix[j+1:]:
                d=np.abs(Z[i]-Z[k])
                if d.max()<=CALIPER:cand.append((float(d@d),i,k))
        used=set()
        for _,i,k in sorted(cand):
            if i in used or k in used:continue
            used.update((i,k));out.append((i,k))
    return out
def swapped(X,ps,rng,j,k,cold):
    out=X.copy()
    for i,h in ps:
        if rng.integers(0,2):
            out[i,j],out[h,j]=X[h,j],X[i,j]
    out[:,k]=out[:,j]*cold
    if not np.array_equal(np.sort(X[:,j]),np.sort(out[:,j])):raise ValueError("Exposure multiset changed")
    return out
def statistics(c):
    a,b=np.asarray(c[TARGET]),np.asarray(c[COLD])
    return {"rail_cold":dict(zip(RESP,map(float,b))),"rail_at_5C":dict(zip(RESP,map(float,a+b)))}
def analyze(rows,ps,expected,seed):
    fields=COVERAGE_FIELDS
    g,b=group_indices(rows,"block");w=block_weights(rows)
    X=np.asarray([[r[k] for k in fields+RESP] for r in rows],float)
    actual=statistics(slopes(X,g,w,fields,(TARGET,COLD)))
    if abs(actual["rail_cold"]["among_species_temperature_sd_C"]-expected[0])>1e-9 or abs(actual["rail_at_5C"]["temp_width_C"]-expected[1])>1e-9:raise ValueError("Exact v401 original fit mismatch")
    mov=[(i,j) for i,j in ps if rows[i][TARGET]!=rows[j][TARGET]]
    qa={"n_cells":len(rows),"n_pairs":len(ps),"n_movable_cells":2*len(mov),
        "movable_fraction":2*len(mov)/len(rows),"movable_blocks":len({bkey(rows[i]) for i,j in mov}),
        "movable_cells_by_continent":dict(Counter(rows[i]["continent"] for p in mov for i in p))}
    out={"matching":qa,"actual":actual}
    if qa["movable_fraction"]<.20 or qa["movable_blocks"]<20:
        return dict(out,status="insufficient fixed-caliper exposure mobility; no rank")
    rng=np.random.default_rng(seed);vals=defaultdict(list);changes=[]
    j=fields.index(TARGET);k=fields.index(COLD);cold=np.array([(15-r["bio1"])/10 for r in rows])
    for _ in range(DRAWS):
        Z=swapped(X,ps,rng,j,k,cold)
        changes.append(int(np.count_nonzero(Z[:,j]!=X[:,j])))
        c=slopes(Z,g,w,fields,(TARGET,COLD))
        if c is None:continue
        for stat,dd in statistics(c).items():
            for t,v in dd.items():vals[stat,t].append(v)
    out["status"]="technical conditional swaps, not randomization inference"
    out["mean_changed_cells"]=float(np.mean(changes))
    out["rank"]={s:{t:{"valid":len(vals[s,t]),"actual":v,"fake_mean":float(np.mean(vals[s,t])) if vals[s,t] else None,
        "fake_low_high":[float(min(vals[s,t])),float(max(vals[s,t]))] if vals[s,t] else None,
        "lower_tail_plus_one":(1+sum(x<=v+1e-12 for x in vals[s,t]))/(len(vals[s,t])+1) if len(vals[s,t])>=.95*DRAWS else None,
        "upper_tail_plus_one":(1+sum(x>=v-1e-12 for x in vals[s,t]))/(len(vals[s,t])+1) if len(vals[s,t])>=.95*DRAWS else None}
        for t,v in dd.items()} for s,dd in actual.items()}
    return out
def run(src,orig):
    q=src["source_QA"]
    if src["version"]!="v401-original-shared-taxon-rarefaction-coarse" or orig["version"]!="v399-crossfit-niche-coarse-aggregate" or q["all_source_rows"]!=288721938 or q["eligible_by_year"]["2025"]!=7827936 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original source")
    land80={tuple(r["key"]) for p in orig["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    rows=[r for r in src["cells"] if r["variants"]]
    if len(rows)!=1323 or len({tuple(r["key"]) for r in rows})!=len(rows):raise ValueError("Original shared support mismatch")
    out=[]
    for land in (.5,.8):
        for h in (False,True):
            sub=sorted((r for r in rows if (land==.5 or tuple(r["key"]) in land80) and (not h or high(r["variants"][0]["shared_rarefied"]))),key=lambda r:tuple(r["key"]))
            n,*e=EXPECTED[land,h]
            if len(sub)!=n:raise ValueError("Original v401 support mismatch")
            rr=[[prepare(r["variants"][j]["shared_rarefied"]) for r in sub] for j in range(2)]
            for a,b in zip(*rr):
                if a["key"]!=b["key"] or bkey(a)!=bkey(b) or any(a[k]!=b[k] for k in FEATURES+(TARGET,"trait_coverage")):raise ValueError("Common geography and covariates differ")
            ps=pairs(rr[0]);seed=SEED+int(land*10)+h
            panel={"land":land,"support":"high_effort" if h else "all","n":n,"fits":[]}
            for j,mode in enumerate(MODES):
                panel["fits"].append({"reference":mode,**analyze(rr[j],ps,e[j],seed)})
            out.append(panel);print(json.dumps({"v404_panel":panel},allow_nan=False),flush=True)
    return {"version":"v404-original-shared-rarefied-conditional-rail-swaps","source_QA":q,"panels":out,
        "features":FEATURES,"caliper":CALIPER,"draws":DRAWS,
        "limits":["technical exchange reference not a causal p-value or realistic counterfactual rail geometry",
        "pairs formed without niche outcomes or railway exposure, within 1000km blocks, many cells may remain fixed",
        "actual and fake rail-cold terms recomputed using local temperature",
        "fixed reference sample-size variants share species, not independent ecological replications",
        "no new original observations, no physiological or causal railway-habitat inference"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();o=run(json.loads(a.source.read_text()),json.loads(a.niche.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
