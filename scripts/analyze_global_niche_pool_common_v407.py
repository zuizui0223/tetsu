#!/usr/bin/env python3
"""v407: common-cell decomposition of taxonomy-conditioned niche nulls."""
import argparse,json
from pathlib import Path
import numpy as np
from analyze_global_niche_pool_null_v403 import LEVELS,RESP,MIN_MOVABLE,EXPECTED
from analyze_global_niche_common_support_v400 import MODES,COVERAGE_FIELDS,INTERACTION_FIELDS,COLD,prepare,high
from analyze_global_crossfit_climate_niche_v399 import regression,TARGET
TAGS=("T_variance","width_C","warm_offset_C")
def attach(maps,levels):
    keys=set.intersection(*(set(maps[l]) for l in levels))
    out=[]
    for k in sorted(keys):
        first=maps[levels[0]][k];r=dict(first)
        for lev in levels:
            a=maps[lev][k]
            if any(a[f]!=first[f] for f in COVERAGE_FIELDS+("continent","block","species","observations","observers")):
                raise ValueError("Common geography covariates changed")
            for tag in TAGS:
                if abs(a["observed_"+tag]-first["observed_"+tag])>1e-10:
                    raise ValueError("Observed community outcome changed across nulls")
                r[lev+"_expected_"+tag]=a["expected_"+tag]
                r[lev+"_excess_"+tag]=a["excess_"+tag]
        for before,after in zip(levels,levels[1:]):
            for tag in TAGS:
                r[after+"_minus_"+before+"_"+tag]=r[after+"_expected_"+tag]-r[before+"_expected_"+tag]
        out.append(r)
    return out
def response_names(levels):
    return tuple("observed_"+tag for tag in TAGS)+tuple(lev+"_"+term+"_"+tag for lev in levels for term in ("expected","excess") for tag in TAGS)+tuple(after+"_minus_"+before+"_"+tag for before,after in zip(levels,levels[1:]) for tag in TAGS)
def verify_additive(m,levels):
    for c in m.get("coef",{}).values():
        for tag in TAGS:
            for lev in levels:
                if abs(c["observed_"+tag]-c[lev+"_expected_"+tag]-c[lev+"_excess_"+tag])>1e-8:
                    raise ValueError("Observed expectation residual decomposition failed")
            for before,after in zip(levels,levels[1:]):
                if abs(c[after+"_minus_"+before+"_"+tag]-(c[before+"_excess_"+tag]-c[after+"_excess_"+tag]))>1e-8:
                    raise ValueError("Between null contrast decomposition failed")
def run(src,old,expected):
    q=src["source_QA"]
    if src["version"]!="v403-original-finite-pool-null-coarse" or old["version"]!="v399-crossfit-niche-coarse-aggregate" or q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279" or not q["exact_2025_all_cell_effort_richness_replay"]:
        raise ValueError("Wrong original source")
    land80={tuple(r["key"]) for p in old["panels"] if p["min_land_fraction"]==.8 for r in p["rows"]}
    rows=[prepare(r) for r in src["rows"]]
    if len(rows)!=1323*2*3:raise ValueError("Original v403 aggregate denominator changed")
    ee={(r["land"],r["support"],r["ref"],r["level"]):(r["n"],r["coef"]) for r in expected}
    if len(ee)!=24:raise ValueError("Incomplete verified prior coefficients")
    out=[]
    for land in (.5,.8):
        for h in (False,True):
            support="high_effort" if h else "all"
            for mode in MODES:
                maps={}
                for lev in LEVELS:
                    base=[r for r in rows if r["reference_exclusion"]==mode and r["null_level"]==lev and (land==.5 or r["key"] in land80) and (not h or high(r))]
                    if len(base)!=EXPECTED[land,h]:raise ValueError("Original v403 support changed")
                    rr=[r for r in base if r["pool_cells"]>=3 and r["movable_fraction"]>=MIN_MOVABLE]
                    n,v=ee[land,support,mode,lev]
                    prior=regression(rr,"block",COVERAGE_FIELDS,(TARGET,COLD),RESP,False)
                    if len(rr)!=n or abs(prior["coef"][COLD]["excess_T_variance"]-v)>1e-9:raise ValueError("Original v403 exact replay failed")
                    maps[lev]={r["key"]:r for r in rr}
                    if len(maps[lev])!=n:raise ValueError("Duplicate cells")
                p={"land":land,"support":support,"reference":mode,"original_support":{k:len(v) for k,v in maps.items()},"fits":[]}
                for levels in (LEVELS[:2],LEVELS):
                    a=attach(maps,levels);names=response_names(levels)
                    f={"levels":levels,"n":len(a),"fits":{}}
                    for label,fields in (("competing",INTERACTION_FIELDS),("coverage_sensitivity",COVERAGE_FIELDS)):
                        m=regression(a,"block",fields,(TARGET,COLD),names,True)
                        verify_additive(m,levels);f["fits"][label]=m
                    p["fits"].append(f)
                out.append(p);print(json.dumps({"v407_panel":p},allow_nan=False),flush=True)
    return {"version":"v407-original-common-cell-taxonomic-null-decomposition","source_QA":q,"panels":out,
            "limits":["descriptive decomposition on observed support, not causal mediation or independent phylogenetic inference",
            "different support between two-level and three-level panel; compare within each common set",
            "observed regional pool includes focal cell; not independent available pool",
            "conditional label null and block intervals are descriptive; no causal test of rail exposure"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("source","niche","expected","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.source),read(a.niche),read(a.expected));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()
