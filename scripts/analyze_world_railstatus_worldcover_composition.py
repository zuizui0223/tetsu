#!/usr/bin/env python3
"""Separate vegetation amount from vegetation composition in existing paired
WorldCover fractions. Exact original result replay; no new imagery or sites.

G=V*P where V=tree+shrub+grass cover, P=grass/V. Symmetric decomposition
of active-minus-inactive at near and background separately, then subtract.
Not a causal mediation model. All four areas must have V>=.20.
"""
import argparse,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_world_railstatus_local_worldcover import CLASS,evaluate,interval,MIN_PAIRS,MIN_BLOCKS
V_FIELDS=("tree","shrub","grass")
MIN_V=.20
SIDES=("active","inactive")
ZONES=("near","background")
def components(active,inactive):
    """Exact product difference. V and P already defined on valid support."""
    va,pa=active;vi,pi=inactive
    amount=(va-vi)*(pa+pi)/2
    allocation=(pa-pi)*(va+vi)/2
    actual=va*pa-vi*pi
    if abs(actual-amount-allocation)>1e-12:raise ValueError("Product decomposition identity")
    return actual,amount,allocation,pa-pi,va-vi

def fraction(row,side,zone):
    q=row[side][zone]["fractions"]
    if set(q)!=set(CLASS.values()) or any(not isinstance(v,(float,int)) or not math.isfinite(v) or not 0<=v<=1 for v in q.values()) or abs(sum(q.values())-1)>1e-9:raise ValueError("Invalid original class fractions")
    return q

def decompose(row):
    q={(s,z):fraction(row,s,z) for s in SIDES for z in ZONES}
    vv={k:sum(v[f] for f in V_FIELDS) for k,v in q.items()}
    if min(vv.values())<MIN_V:return None
    pp={k:q[k]["grass"]/v for k,v in vv.items()}
    d={}
    for z in ZONES:
        d[z]=components(*[(vv[s,z],pp[s,z]) for s in SIDES])
    rr=[x-y for x,y in zip(d["near"],d["background"])]
    if abs(rr[0]-rr[1]-rr[2])>1e-12:raise ValueError("Spatial decomposition identity")
    return dict(zip(("grass_total","vegetation_amount_component","within_vegetation_allocation_component","grass_share_in_tree_shrub_grass","tree_shrub_grass_total"),rr))

def stats(rows,targets):
    n=len(rows);nb=len({tuple(r["block1000"]) for r in rows})
    qa={"n_pairs":n,"n_blocks":nb}
    if n<MIN_PAIRS or nb<MIN_BLOCKS:return dict(qa,status="insufficient support")
    effects={k:interval([(tuple(r["block1000"]),r["value"][k]) for r in rows]) for k in targets}
    return dict(qa,status="estimated",effects=effects)

def raw(row):
    q={(s,z):fraction(row,s,z) for s in SIDES for z in ZONES}
    def did(f):
        return sum(w*sum(q[s,z][k] for k in f) for s,z,w in (("active","near",1),("inactive","near",-1),("active","background",-1),("inactive","background",1)))
    return {"built":did(("built",)),"grass":did(("grass",)),"tree_shrub_grass_total":did(V_FIELDS),"all_vegetation_classes":did(("tree","shrub","grass","crop","wetland","mangrove","moss_lichen")),"bare":did(("bare",))}

def check_effects(replayed,old):
    if replayed["status"]!=old["status"] or replayed["n_complete_pairs"]!=old["n_complete_pairs"] or replayed["n_blocks"]!=old["n_blocks"]:raise ValueError("Original support not replayed")
    if replayed["status"]!="estimated":raise ValueError("Prior model not estimated")
    for name,z in replayed["effects"].items():
        for k,v in z.items():
            p=old["effects"][name][k]
            if abs(v["estimate"]-p["estimate"])>1e-10 or max(abs(a-b) for a,b in zip(v["interval95_approx"],p["interval95_approx"]))>1e-10:raise ValueError("Prior result changed")
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--historical",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    src=json.loads(a.source.read_text());hist=json.loads(a.historical.read_text())
    if src["version"]!="railstatus-local-worldcover-2021-v1" or src["source"]!="ESA WorldCover 2021 v200" or src["geometry"]["draws"]!=8192 or src["geometry"]["n_selected"]!=126 or len(src["pairs"])!=126 or len({r["draw"] for r in src["pairs"]})!=126 or len(src["tile_receipts"])!=103:raise ValueError("Wrong frozen WorldCover source")
    if hist["version"]!="railstatus-worldcover-historical-map-v1" or hist["blocked"]!="Historical API HTTP 502" or hist["support"]["complete_two_date_pairs"]!=0:raise ValueError("Prior history audit changed; review before fitting")
    check_effects(evaluate(src["pairs"]),src["primary"])
    complete=[r for r in src["pairs"] if all(r[k]["status"]=="complete" for k in SIDES)]
    if len(complete)!=123 or len({tuple(r["block1000"]) for r in complete})!=45:raise ValueError("Prior complete sample differs")
    all_rows=[dict(r,value=raw(r)) for r in complete]
    decomposed=[]
    for r in complete:
        v=decompose(r)
        if v is not None:decomposed.append(dict(r,value=v))
    keys={r["draw"] for r in decomposed}
    targets=("grass_total","vegetation_amount_component","within_vegetation_allocation_component","grass_share_in_tree_shrub_grass","tree_shrub_grass_total")
    matched=stats(decomposed,targets)
    if matched["status"]=="estimated":
        e=matched["effects"]
        if abs(e["grass_total"]["estimate"]-e["vegetation_amount_component"]["estimate"]-e["within_vegetation_allocation_component"]["estimate"])>1e-12:raise ValueError("Weighted decomposition identity")
    out={"version":"railstatus-worldcover-composition-decomposition-v1","source_sha256":hashlib.sha256(a.source.read_bytes()).hexdigest(),"historical_sha256":hashlib.sha256(a.historical.read_bytes()).hexdigest(),"history_audit_status":"blocked at preflight (HTTP 502); no historical confirmation","exact_original_result_replay":True,
         "primary":matched,"unrestricted_raw":stats(all_rows,tuple(all_rows[0]["value"])),
         "same_support_raw":stats([r for r in all_rows if r["draw"] in keys],tuple(all_rows[0]["value"])),
         "n_excluded_low_vegetation":len(complete)-len(decomposed),"min_tree_shrub_grass_fraction_each_area":MIN_V,
         "limits":["conditional exploratory decomposition after prior satellite results; not independent confirmation",
                   "tree+shrub+grass is a fixed class subtotal, not all vegetation or biomass",
                   "no imputation when vegetation denominator small; selected subset differs from full sample",
                   "fractional cover composition closure; components are arithmetic not causal",
                   "spatial near-background contrast, not temporal difference-in-differences or before-after closure",
                   "2021 satellite imagery vs 2026 map status; historical validation blocked, no measured mowing or wind, no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({"worldcover_composition_result":out},allow_nan=False),flush=True)
if __name__=="__main__":main()
