#!/usr/bin/env python3
"""v389: numerator/denominator decomposition of v383 regional documented fraction."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_global_rail_naturalized_road_tree_v383 import (
    ORIGINAL_REPLAY,eligible_region_rows,regression,ADJ,ENV)
from analyze_global_rail_documented_naturalized_v379 import within_region
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights,RICH_BASE
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import RICH_BUILT
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
FIELDS=RICH_BASE+(RICH_BUILT,)+ENV+("log_rail_km",)
OUTCOMES=("log1p_documented","log1p_unlisted","log1p_count_contrast","known_fraction_pp")
TARGETS=("log_rail_km","log_road_density")
BOOT=399
SEED=38920261125
EXPECTED={
 (100000,.3,.5):(723,393,.4428467087806793),
 (100000,.3,.8):(705,379,.4622420610746872),
 (100000,0.,.5):(1178,641,.5269058952435043),
 (100000,0.,.8):(1153,619,.510326182091794),
 (200000,.3,.5):(130,57,-.015630918405012708),
 (200000,.3,.8):(127,52,.18188035747435277),
 (200000,0.,.5):(240,87,.2971586451745747),
 (200000,0.,.8):(235,81,.4248486362584915)}
def count_outcomes(row):
    n=row["vascular_richness"];a=row["documented_richness"]
    if any(type(x)!=int for x in (n,a)) or not 0<=a<=n<=row["species"] or n<20:
        raise ValueError("Invalid original aggregate documented/vascular counts")
    if abs(row["known_fraction"]-a/n)>1e-12:raise ValueError("Original fraction mismatch")
    b=n-a
    x,y=math.log1p(a),math.log1p(b)
    return (x,y,x-y,100*a/n)
def source_rows(panel):
    r=[]
    for raw in panel["rows"]:
        x=dict(raw,key=tuple(raw["key"]),block=tuple(raw["block"]))
        y=count_outcomes(x);x.update(dict(zip(OUTCOMES,y)))
        r.append(x)
    if len({x["key"] for x in r})!=len(r):raise ValueError("Duplicate original cell")
    k=(panel["width_m"],panel["region_guard_degrees"])
    n,b=ORIGINAL_REPLAY[k]
    if len(r)!=n or abs(within_region(r)["coefficient_pp_per_log1p_rail_km"]-b)>1e-9:
        raise ValueError("Original v379d/v383 cell aggregate replay mismatch")
    return r
def conditional(X,groups,w,fields):
    if not np.isfinite(w).all() or w.sum()<=0 or np.any(w<0):return None
    w=w/w.sum();Z=X.copy()
    for ix in groups:
        if w[ix].sum()>0:Z[ix]-=np.average(Z[ix],axis=0,weights=w[ix])
    Z*=np.sqrt(w)[:,None]
    if not np.isfinite(Z).all():raise ValueError("Nonfinite original model values")
    n=len(fields);result={}
    for k in TARGETS:
        j=fields.index(k);x=Z[:,j];Y=Z[:,n:]
        other=np.delete(Z[:,:n],j,axis=1)
        scale=np.sqrt((other**2).sum(axis=0));ok=scale>1e-12
        other=other[:,ok]/scale[ok]
        rx=x-other@np.linalg.lstsq(other,x,rcond=1e-10)[0]
        ry=Y-other@np.linalg.lstsq(other,Y,rcond=1e-10)[0]
        var=float(rx@rx)
        if var<=1e-10:return None
        b=rx@ry/var
        if abs(b[0]-b[1]-b[2])>1e-9:raise ValueError("Count component log-contrast coefficient not additive")
        result[k]=dict(zip(OUTCOMES,map(float,b)))
    return result
def fit(rows,fields,boot=True):
    groups=defaultdict(list);blocks=defaultdict(list)
    for i,r in enumerate(rows):
        groups[r["region"]].append(i);blocks[(r["continent"],r["block"])].append(i)
    out={"n_cells":len(rows),"n_regions":len(groups),"n_1000km_blocks":len(blocks),
         "cells_by_continent":dict(Counter(r["continent"] for r in rows)),"coefficients":None}
    if len(groups)<2 or len(rows)<30 or len(blocks)<5:return out
    X=np.array([[r[k] for k in fields+OUTCOMES] for r in rows],float)
    ids=[np.array(v,dtype=int) for v in groups.values()]
    c=conditional(X,ids,block_weights(rows),fields);out["coefficients"]=c
    if not c or not boot:return out
    bycont=defaultdict(list)
    for b in sorted(blocks):bycont[b[0]].append(b)
    rng=np.random.default_rng(SEED)
    samples={k:{t:[] for t in OUTCOMES} for k in TARGETS}
    for i in range(BOOT):
        w=np.zeros(len(rows))
        for bs in bycont.values():
            mult=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for b,m in zip(bs,mult):
                ix=blocks[b];w[ix]=m/len(ix)
        z=conditional(X,ids,w,fields)
        if z is None:continue
        for k in TARGETS:
            for t in OUTCOMES:samples[k][t].append(z[k][t])
    out["block_bootstrap_percentile_95_descriptive"]={k:{t:{
        "valid":len(v),"low":float(np.quantile(v,.025)) if len(v)>=.95*BOOT else None,
        "high":float(np.quantile(v,.975)) if len(v)>=.95*BOOT else None}
        for t,v in ss.items()} for k,ss in samples.items()}
    return out
def run(cells,road,tree):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:
        raise ValueError("Not verified v383 original aggregate")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or not road["version"].startswith("v381b"):
        raise ValueError("Not source v381b road")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:
        raise ValueError("Not source v380b tree")
    panels=[]
    for p in cells["panels"]:
        width,pad=p["width_m"],p["region_guard_degrees"];rr=source_rows(p)
        rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
        for gate in LAND_THRESHOLDS:
            matched,qa=enrich_cells(rr,rd,tr,width,gate)
            whole=eligible_region_rows(matched)
            n,m,b=EXPECTED[width,pad,gate]
            old=regression(whole,ADJ,("log_rail_km",),False)
            if len(matched)!=n or len(whole)!=m or abs(old["coefficients_pp"]["log_rail_km"]-b)>1e-9:
                raise ValueError("Exact v383 adjusted sample/estimate mismatch")
            out={"width_m":width,"region_guard_degrees":pad,"min_land_fraction":gate,
                 "matched":n,"v383_original_adjusted_rail_fraction_pp":b,"fits":[]}
            for mode in ("all","high_effort"):
                s=matched if mode=="all" else [r for r in matched if r["observations"]>=200 and r["observers"]>=20]
                s=eligible_region_rows(s)
                # Conditioning on the observed total can mechanically alter
                # component slopes. Primary count/ratio model omits total
                # richness; original conditional fraction is retained separately.
                base=fit(s,FIELDS)
                conditional_total=regression(s,ADJ,("log_rail_km","log_road_density"),True)
                out["fits"].append({"support":mode,"component_and_fraction_no_total_richness_control":base,
                    "fraction_with_original_total_richness_control":conditional_total})
            panels.append(out)
            print(json.dumps({"v389_panel":out},allow_nan=False),flush=True)
    return {"version":"v389-v383-verified-documented-vs-unlisted-count-decomposition",
        "source_replay":cells["source_replay"],"panels":panels,"bootstrap_draws":BOOT,"seed":SEED,
        "limits":["GloNAF documented naturalized counts vs nonlisted counts, NOT alien versus native",
        "observed species richness only; log1p count contrast is not exact log odds when counts small",
        "conditional on observed regional support, no standardized detectability or causal transport",
        "high effort subset changes geography, report actual composition; not independent replication",
        "block percentile stability is descriptive, not proof of external generalization"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.cells),read(a.road),read(a.tree));a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
