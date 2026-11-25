#!/usr/bin/env python3
"""v411 NOAA CPC original annual precipitation change on frozen v410 cell support."""
import argparse,json,calendar,hashlib,math
from pathlib import Path
from collections import Counter
import numpy as np
from netCDF4 import Dataset,num2date
from pyproj import Transformer
from analyze_global_niche_annual_change_v410 import make_change,YEARS,RESP,STAT,stat
from analyze_global_niche_common_support_v400 import MODES,prepare,high,INTERACTION_FIELDS,COVERAGE_FIELDS,COLD
from analyze_global_crossfit_climate_niche_v399 import TARGET,group_indices,slopes,regression,W
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights,cell_center
from generate_global_equal_area_climate_land_grid_v351 import sha_file
BASE_URL="https://psl.noaa.gov/thredds/fileServer/Datasets/cpc_global_precip/"
MODIFIERS=(TARGET,"log_road_density","log1p_1000_times_original_2020_builtup_fraction","log_population_density","tree_consensus_fullcell")
PREC="precip_decline_log"
DRY_INT=TARGET+"_x_"+PREC
MIN_VALID=1.
def aggregate_annual(path,year,keys):
    with Dataset(path) as ds:
        p=ds.variables["precip"];lat=np.asarray(ds.variables["lat"][:],float);lon=np.asarray(ds.variables["lon"][:],float)
        if p.dimensions!=("time","lat","lon") or len(lat)!=360 or len(lon)!=720 or not np.allclose(np.abs(np.diff(lat)),.5) or not np.allclose(np.diff(lon),.5):raise ValueError("Unexpected original CPC grid")
        unit=str(p.units).lower().strip()
        if unit not in ("mm","mm/day","mm day-1","mm/d","mm day^-1","mm d-1","mm/day "):raise ValueError("Unexpected original CPC precip units: "+unit)
        tv=ds.variables["time"];dates=num2date(tv[:],tv.units,calendar=getattr(tv,"calendar","standard"))
        expected=365+calendar.isleap(year)
        days=[(d.year,d.month,d.day) for d in dates]
        if len(days)!=expected or len(set(days))!=expected or any(y!=year for y,m,d in days):raise ValueError("Incomplete original calendar")
        tf=Transformer.from_crs(6933,4326,always_xy=True)
        coords=[tf.transform(*cell_center(k[0],k[1],W)) for k in keys]
        iy=np.array([np.argmin(np.abs(lat-y)) for x,y in coords]);ix=np.array([np.argmin(np.abs(((lon-x+180)%360)-180)) for x,y in coords])
        if any(abs(lat[i]-y)>.251 or abs(((lon[j]-x+180)%360)-180)>.251 for (x,y),i,j in zip(coords,iy,ix)):raise ValueError("CPC nearest-cell mismatch")
        sums=np.zeros(len(keys));n=np.zeros(len(keys),int)
        for start in range(0,expected,16):
            a=np.ma.asarray(p[start:min(start+16,expected),:,:],dtype=float).filled(np.nan)
            v=a[:,iy,ix];valid=np.isfinite(v)&(v>=0)&(v<10000)
            sums+=np.where(valid,v,0).sum(axis=0);n+=valid.sum(axis=0)
        out={k:float(s) for k,s,c in zip(keys,sums,n) if c==expected}
        return out,{"year":year,"url":BASE_URL+path.name,"sha256":sha_file(path),"size_bytes":path.stat().st_size,
                    "calendar_days":expected,"valid_all_days_cells":len(out),"requested_cells":len(keys),"units":unit,
                    "valid_days_min":int(n.min()),"valid_days_max":int(n.max()),"source_title":str(getattr(ds,"title",""))}
def add_weather(a,weather,contrast):
    r=[]
    for x in a:
        k=x["key"]
        if any(k not in weather[y] for y in YEARS):continue
        z=dict(x)
        prev=(weather["2023"][k],weather["2024"][k]) if contrast=="2025_minus_earlier_mean" else (weather["2023"][k],)
        cur=weather["2025" if contrast=="2025_minus_earlier_mean" else "2024"][k]
        z[PREC]=float(np.mean(np.log1p(prev)))-math.log1p(cur)
        z["cpc_baseline_log_precip"]=float(np.mean(np.log1p(prev)))
        for m in MODIFIERS:z[m+"_x_"+PREC]=z[m]*z[PREC]
        r.append(z)
    return r
def run(src,niche,old,weather,qa):
    q=src["source_QA"]
    if src["version"]!="v405-original-annual-niche-coarse" or niche["version"]!="v399-crossfit-niche-coarse-aggregate" or old["source_run"]!=45572691489 or q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:raise ValueError("Wrong original source")
    from analyze_global_niche_annual_change_v410 import make_change
    land80={tuple(r["key"]) for p in niche["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    data={tuple(r["key"]):r["annual"] for r in src["cells"]}
    if len(data)!=1323:raise ValueError("Wrong v405 cells")
    out=[]
    for p in old["panels"]:
        land=p["land"];h=p["support"]=="high_effort";j=MODES.index(p["reference"])
        keys=sorted(k for k,v in data.items() if (land==.5 or k in land80) and all(v[y] is not None and (not h or high(v[y][0])) for y in YEARS))
        if len(keys)!=p["n"]:raise ValueError("Original v410 support changed")
        rr={y:[prepare(data[k][y][j]) for k in keys] for y in YEARS}
        panel={"land":land,"support":p["support"],"reference":p["reference"],"n_original":len(keys),"fits":[]}
        for f in p["fits"]:
            fields=INTERACTION_FIELDS if f["model"]=="competing" else COVERAGE_FIELDS
            a,fs=make_change(rr,fields,f["contrast"])
            g,_=group_indices(a,"block");X=np.asarray([[r[k] for k in fs+RESP] for r in a],float)
            c=stat(X,g,block_weights(a),fs)
            if c is None or any(abs(c[i,j]-f["all"][s][t])>1e-9 for i,s in enumerate(STAT) for j,t in enumerate(RESP)):raise ValueError("Exact v410 prior fit not reproduced")
            a=add_weather(a,weather,f["contrast"])
            ff=fs+("cpc_baseline_log_precip",PREC)
            # No uncalibrated single-year climate proxy invented for missing cells.
            m={"model":f["model"],"contrast":f["contrast"],"n_complete_precip":len(a),
               "precip_decline_quantiles":list(map(float,np.quantile([r[PREC] for r in a],[0,.05,.5,.95,1]))) if a else None,
               "n_positive_decline":sum(r[PREC]>0 for r in a),"fits":{}}
            for label,fields in (("original_common_weather",fs),("annual_precip_main",ff),("precip_competing_interactions",ff+tuple(k+"_x_"+PREC for k in MODIFIERS))):
                targets=(TARGET,COLD,PREC,DRY_INT) if label=="precip_competing_interactions" else (TARGET,COLD,PREC) if label=="annual_precip_main" else (TARGET,COLD)
                m["fits"][label]=regression(a,"block",fields,targets,RESP,True)
            panel["fits"].append(m)
        out.append(panel);print(json.dumps({"v411_panel":panel},allow_nan=False),flush=True)
    return {"version":"v411-original-CPC-precipitation-change-niche","source_QA":q,"precipitation_source_QA":qa,"panels":out,
       "limits":["CPC gauge-interpolated annual rainfall at nearest 0.5-degree grid center; not actual rail microclimate, PET, drought stress or growing-season moisture",
       "annual totals may include post-observation rain; contemporaneous association, not temporal causal ordering",
       "fixed 2025-selected observational geography, same three years and historical realized niche reference",
       "posthoc mechanism sensitivity with competing interactions; no independent confirmation or causal interpretation",
       "source downloads recorded by SHA and size, all calendar days valid per cell; no missing-day imputation"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","verified","climate_dir","out","climate_out"):p.add_argument("--"+k.replace("_","-"),dest=k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    src=read(a.source);keys=sorted(tuple(r["key"]) for r in src["cells"]);weather={};qa=[]
    for y in YEARS:
        weather[y],z=aggregate_annual(a.climate_dir/("precip."+y+".nc"),int(y),keys)
        qa.append(z);print(json.dumps({"v411_original_CPC":z}),flush=True)
    a.climate_out.parent.mkdir(parents=True,exist_ok=True)
    a.climate_out.write_text(json.dumps({"source_QA":qa,"annual_cell_precipitation_mm":[{"key":k,**{y:weather[y].get(k) for y in YEARS}} for k in keys]},allow_nan=False))
    o=run(src,read(a.niche),read(a.verified),weather,qa);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
