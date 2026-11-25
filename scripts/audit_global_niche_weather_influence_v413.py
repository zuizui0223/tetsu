#!/usr/bin/env python3
"""v413 influence audit of the actual rainfall-adjusted annual change result."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
from analyze_global_niche_precip_change_v411 import add_weather,PREC,DRY_INT,MODIFIERS
from analyze_global_niche_annual_change_v410 import make_change,YEARS,RESP
from analyze_global_niche_common_support_v400 import MODES,prepare,high,INTERACTION_FIELDS,COVERAGE_FIELDS,COLD
from analyze_global_crossfit_climate_niche_v399 import TARGET,group_indices,slopes
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
TGT=(TARGET,COLD,PREC,DRY_INT)
D="dry_offset_log_precip"
def evaluate(rows,fields):
    g,b=group_indices(rows,"block")
    n=len(rows)
    o={"n_cells":n,"n_blocks":len(b)}
    if n<100 or len(b)<20 or n-len(g)-len(fields)<30:return dict(o,status="insufficient support")
    X=np.asarray([[r[k] for k in fields+RESP] for r in rows],float)
    c=slopes(X,g,block_weights(rows),fields,TGT)
    if c is None:return dict(o,status="insufficient residual variation")
    return dict(o,status="estimated",coef={k:dict(zip(RESP,v)) for k,v in c.items()})
def summary(f):
    c=f["coef"]
    return {"dry_rail_15C_zero_decline":c[TARGET][D],
            "dry_rail_5C_zero_decline":c[TARGET][D]+c[COLD][D],
            "dry_rail_x_precip_decline":c[DRY_INT][D]}
def audit(rows,fields,old):
    original=evaluate(rows,fields)
    if original["status"]!="estimated" or old["status"]!="estimated":raise ValueError("Original fit unsupported")
    if original["n_cells"]!=old["n_cells"] or original["n_blocks"]!=old["n_blocks"]:raise ValueError("Original support changed")
    for k in TGT:
        for t in RESP:
            if abs(original["coef"][k][t]-old["coef"][k][t])>1e-9:raise ValueError("Exact v411 coefficient replay failed")
    full=summary(original)
    def remove(keys,continent=False):
        out=[]
        for k in keys:
            kept=[r for r in rows if (r["continent"] if continent else (r["continent"],r["block"]))!=k]
            f=evaluate(kept,fields)
            z={"omitted":k,"n_cells":f["n_cells"],"n_blocks":f["n_blocks"],"status":f["status"]}
            if f["status"]=="estimated":z["statistics"]=summary(f)
            out.append(z)
        return out
    b=remove(sorted({(r["continent"],r["block"]) for r in rows}))
    c=remove(sorted({r["continent"] for r in rows}),True)
    v={}
    for key,value in full.items():
        x=[(r["omitted"],r["statistics"][key]) for r in b if "statistics" in r]
        if not x:v[key]={"valid":0};continue
        maxv=max(x,key=lambda p:abs(p[1]-value))
        v[key]={"valid":len(x),"minimum":min(t for _,t in x),"maximum":max(t for _,t in x),
            "positive":sum(t>0 for _,t in x),"negative":sum(t<0 for _,t in x),
            "largest_change_block":maxv[0],"largest_absolute_change":abs(maxv[1]-value)}
    return {"full":full,"block_deletion_summary":v,"block_deletions":b,"continent_deletions":c}
def run(src,niche,weather_src,old):
    q=src["source_QA"]
    if src["version"]!="v405-original-annual-niche-coarse" or niche["version"]!="v399-crossfit-niche-coarse-aggregate" or old["source_run"]!=45573280418 or q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original source")
    if weather_src["source_QA"]!=old["qa"]:raise ValueError("Original CPC provenance changed")
    data={tuple(r["key"]):r["annual"] for r in src["cells"]}
    if len(data)!=1323 or len(weather_src["annual_cell_precipitation_mm"])!=1323:raise ValueError("Original cell support changed")
    weather={y:{tuple(r["key"]):r[y] for r in weather_src["annual_cell_precipitation_mm"] if r[y] is not None} for y in YEARS}
    if [len(weather[y]) for y in YEARS]!=[1289,1317,1318]:raise ValueError("Original CPC support changed")
    land80={tuple(r["key"]) for p in niche["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    out=[]
    for p in old["panels"]:
        land=p["land"];h=p["support"]=="high_effort";j=MODES.index(p["reference"])
        keys=sorted(k for k,v in data.items() if (land==.5 or k in land80) and all(v[y] is not None and (not h or high(v[y][0])) for y in YEARS))
        if len(keys)!=p["n_original"]:raise ValueError("Original annual support changed")
        rr={y:[prepare(data[k][y][j]) for k in keys] for y in YEARS}
        a={"land":land,"support":p["support"],"reference":p["reference"],"fits":[]}
        for f in p["fits"]:
            fields=INTERACTION_FIELDS if f["model"]=="competing" else COVERAGE_FIELDS
            x,fs=make_change(rr,fields,f["contrast"])
            x=add_weather(x,weather,f["contrast"])
            if len(x)!=f["n_complete_precip"]:raise ValueError("Original weather support changed")
            ff=fs+("cpc_baseline_log_precip",PREC)+tuple(k+"_x_"+PREC for k in MODIFIERS)
            z=audit(x,ff,f["fits"]["precip_competing_interactions"])
            a["fits"].append({"model":f["model"],"contrast":f["contrast"],**z})
        out.append(a)
        print(json.dumps({"v413_panel":{"land":land,"support":p["support"],"reference":p["reference"],"fits":[{k:v for k,v in f.items() if k!="block_deletions"} for f in a["fits"]]}},allow_nan=False),flush=True)
    return {"version":"v413-v411-actual-weather-adjusted-influence-audit","source_QA":q,"precip_QA":weather_src["source_QA"],"panels":out,
       "limits":["deletion slopes are same-data influence diagnostics, not independent replications or significance tests",
       "no threshold relaxation: continent exclusions can fail existing 100-cell/20-block gate",
       "no new climate, observational or causal claim; annual precipitation may postdate observations"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","weather","verified","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text());o=run(read(a.source),read(a.niche),read(a.weather),read(a.verified));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
