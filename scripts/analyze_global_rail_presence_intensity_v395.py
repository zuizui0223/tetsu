#!/usr/bin/env python3
"""v395: separate mapped rail presence and within-positive rail intensity."""
import argparse,json,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_naturalized_count_decomposition_v389 import source_rows,OUTCOMES,BOOT,EXPECTED
from analyze_global_naturalized_quality_effort_v390 import spatial_key,high
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
from analyze_global_naturalized_population_v392 import add_population
from analyze_global_naturalized_local_swap_v393 import supported_local
from analyze_global_naturalized_joint_population_local_v394 import (
    POP_SHA,FULL,LOCAL_EXPECT,groups,mat,supported)
from analyze_global_naturalized_count_decomposition_v389 import conditional
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
BASE=tuple(k for k in FULL if k!="log_rail_km")
SEED=39520261125
REF_KM=100.
PRIOR={
 (100000,.3,.5):(.003640234533337464,.04401791006162525),
 (100000,.3,.8):(.006899656442287413,.04347667097520801),
 (100000,0.,.5):(.02955371004598179,.035506534438446764),
 (100000,0.,.8):(.030669427782566573,.03307246733237476)}
def add_exposure(row):
    km=float(row["rail_km"])
    if not math.isfinite(km) or km<0 or abs(row["log_rail_km"]-math.log1p(km))>1e-10:
        raise ValueError("Invalid original exposure")
    return dict(row,rail_present=float(km>0),
                rail_positive_centered=math.log1p(km)-math.log1p(REF_KM) if km>0 else 0.)
def support_summary(rows):
    g,b=groups(rows)
    return {"n_cells":len(rows),"n_regions":len({r["region"] for r in rows}),
        "n_region_block_strata":len(g),"n_blocks":len(b),
        "n_zero":sum(r["rail_present"]==0 for r in rows),
        "n_positive":sum(r["rail_present"]==1 for r in rows),
        "cells_by_continent":dict(Counter(r["continent"] for r in rows))}
def mixed(rows):
    g=defaultdict(list)
    for r in rows:g[spatial_key(r)].append(r)
    keys={k for k,rr in g.items() if {r["rail_present"] for r in rr}=={0.,1.}}
    return [r for r in rows if spatial_key(r) in keys]
def coefficients(X,g,w,fields,targets):
    if not np.isfinite(w).all() or w.sum()<=0 or np.any(w<0):return None
    w=w/w.sum();Z=X.copy()
    for ix in g:
        if w[ix].sum()>0:Z[ix]-=np.average(Z[ix],axis=0,weights=w[ix])
    Z*=np.sqrt(w)[:,None]
    n=len(fields);A=Z[:,:n];Y=Z[:,n:]
    if not np.isfinite(Z).all():raise ValueError("Nonfinite regression matrix")
    result={}
    for k in targets:
        j=fields.index(k);x=A[:,j];other=np.delete(A,j,axis=1)
        scale=np.sqrt((other**2).sum(axis=0));ok=scale>1e-12
        other=other[:,ok]/scale[ok]
        rx=x-other@np.linalg.lstsq(other,x,rcond=1e-10)[0]
        var=float(rx@rx)
        if var<=1e-10:return None
        ry=Y-other@np.linalg.lstsq(other,Y,rcond=1e-10)[0]
        b=rx@ry/var
        if abs(b[0]-b[1]-b[2])>1e-9:raise ValueError("Count contrast slopes do not add")
        result[k]=dict(zip(OUTCOMES,map(float,b)))
    return result
def model(rows,fields,targets,require_mixed,seed):
    out=support_summary(rows)
    if not supported(rows,fields) or (require_mixed and min(out["n_zero"],out["n_positive"])<30):
        out["status"]="insufficient prespecified support";return out
    g,blocks=groups(rows);X=mat(rows,fields)
    real=coefficients(X,g,block_weights(rows),fields,targets)
    if real is None:out["status"]="insufficient residual exposure variation";return out
    out["status"]="estimated";out["coefficients"]=real
    by=defaultdict(list)
    for b in sorted(blocks):by[b[0]].append(b)
    rng=np.random.default_rng(seed)
    samples={k:{t:[] for t in OUTCOMES} for k in targets}
    for _ in range(BOOT):
        bw=np.zeros(len(rows))
        for bs in by.values():
            mult=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for b,m in zip(bs,mult):
                ix=blocks[b];bw[ix]=m/len(ix)
        a=coefficients(X,g,bw,fields,targets)
        if a is None:continue
        for k in targets:
            for t in OUTCOMES:samples[k][t].append(a[k][t])
    out["block_percentile_95_descriptive"]={k:{t:{"valid":len(v),
        "low":float(np.quantile(v,.025)) if len(v)>=.95*BOOT else None,
        "high":float(np.quantile(v,.975)) if len(v)>=.95*BOOT else None}
        for t,v in ss.items()} for k,ss in samples.items()}
    for label,units,keyfunc in (("leave_one_block_out",sorted(blocks),lambda x:(x["continent"],x["block"])),
                                ("leave_one_continent_out",sorted({r["continent"] for r in rows}),lambda x:x["continent"])):
        v={k:{t:[] for t in OUTCOMES} for k in targets}
        for u in units:
            r=[x for x in rows if keyfunc(x)!=u]
            if not supported(r,fields):continue
            if require_mixed and min(sum(x["rail_present"]==0 for x in r),sum(x["rail_present"]==1 for x in r))<30:continue
            gg,_=groups(r);c=coefficients(mat(r,fields),gg,block_weights(r),fields,targets)
            if c is None:continue
            for k in targets:
                for t in OUTCOMES:v[k][t].append(c[k][t])
        out[label]={"n_attempted":len(units),"ranges":{k:{t:{"valid":len(a),"min":min(a) if a else None,
            "max":max(a) if a else None,"positive":sum(x>0 for x in a)} for t,a in ss.items()} for k,ss in v.items()}}
    return out
def panel(rows,key,i):
    local=[add_exposure(r) for r in supported_local(rows)]
    n,old=LOCAL_EXPECT[key][i]
    if len(local)!=n:raise ValueError("v394 local denominator changed")
    o={"original_local":support_summary(local),"models":{}}
    if old is None:o["status"]="v394 local support insufficient";return o
    g,_=groups(local)
    base=conditional(mat(local,FULL),g,block_weights(local),FULL)
    if base is None or abs(base["log_rail_km"]["log1p_documented"]-PRIOR[key][i])>1e-9:
        raise ValueError("Exact v394 original coefficient not reproduced")
    o["original_joint_log_rail_coefficients"]=base
    s=mixed(local)
    # Intensity estimand uses only rail-positive cells and rechecks within-stratum support.
    p=supported_local([r for r in local if r["rail_present"]==1])
    seed=SEED+key[0]+int(key[1]*100)+int(key[2]*10)+i
    o["models"]["binary_presence_only"]=model(s,BASE+("rail_present",),("rail_present",),True,seed)
    o["models"]["presence_and_intensity_at_100km_reference"]=model(s,BASE+("rail_present","rail_positive_centered"),("rail_present","rail_positive_centered"),True,seed)
    o["models"]["positive_cells_intensity_only"]=model(p,FULL,("log_rail_km",),False,seed)
    o["status"]="completed"
    return o
def run(cells,roads,tree,pop):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:raise ValueError("Wrong source cells")
    if roads["source_zip_MD5"]!=MD5 or roads["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong roads")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree")
    if pop["version"]!="v391-original-JRC-GHS-POP-2020-1km-global-counts" or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong population")
    out=[]
    for p in cells["panels"]:
        w,pad=p["width_m"],p["region_guard_degrees"];r=source_rows(p)
        rd=make_lookup(roads["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
        for land in LAND_THRESHOLDS:
            k=w,pad,land
            m,_=enrich_cells(r,rd,tr,w,land);z,_=add_population(m,pp,rd)
            if len(m)!=EXPECTED[k][0] or len(z)!=len(m):raise ValueError("Original v394 source support changed")
            x={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":land,"fits":[]}
            for i,mode in enumerate(("all","high_effort")):
                rr=eligible_region_rows(z if i==0 else [r for r in z if high(r)])
                x["fits"].append({"support":mode,"test":panel(rr,k,i)})
            out.append(x);print(json.dumps({"v395_panel":x},allow_nan=False),flush=True)
    return {"version":"v395-original-rail-presence-versus-positive-length",
        "panels":out,"population_SHA256":POP_SHA,"bootstrap_draws":BOOT,"positive_length_reference_km":REF_KM,
        "limits":["zero denotes no mapped candidate rail length, not certified railway absence",
          "binary comparison among mixed strata; joint presence coefficient contrasts 100 km versus zero, not a literal onset effect",
          "positive-only intensity conditions on mapped rail presence and may change support",
          "same 2025 observed richness records, not temporal change or causal transport",
          "unlisted is not native; no independent replication or standardized detection"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","population","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    r=run(read(a.cells),read(a.road),read(a.tree),read(a.population));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,allow_nan=False))
if __name__=="__main__":main()
