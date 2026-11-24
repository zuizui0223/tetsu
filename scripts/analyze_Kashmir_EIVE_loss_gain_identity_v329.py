#!/usr/bin/env python3
"""v329 exact recorded plant-taxon thermal CWM turnover decomposition in Kashmir.

Paired original 2014/2017 railway Track, Embankment and two Interior plots.
For persistent species P, lost L, gained G, with equal-taxon EIVE CWM:
  delta(2017-2014) = |G|/|P+G|*(mean(G)-mean(P))
                    - |L|/|P+L|*(mean(L)-mean(P)).
No claim of seed dispersal, demographic extinction, causality, physiological
temperature range change, or independent national-climate replication.
"""
import argparse
import json
import math
from collections import defaultdict,Counter
from pathlib import Path
from statistics import mean,median

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import (
    species_names_from_publisher,cohort)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table

ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
SYSTEMS=3
MIN_RECORDED_MATCHED=(1,3)

def m(arr):
    return mean(arr) if arr else None

def vals(indices,names,eive,metric):
    if indices is None:return None
    return {names[i]:float(eive[names[i]][metric]) for i in indices
           if names[i] in eive and eive[names[i]]["systems"]>=SYSTEMS}

def components(a,b):
    """Return additive gain/loss components or None if no persistent anchor."""
    if a is None or b is None or not a or not b:return None
    pa=set(a)&set(b)
    lost=set(a)-set(b)
    gains=set(b)-set(a)
    if not pa:return None
    pm=m([a[x] for x in pa])
    am=m(list(a.values()));bm=m(list(b.values()))
    losing=(-len(lost)/len(a)*(m([a[x] for x in lost])-pm)) if lost else 0.0
    adding=(len(gains)/len(b)*(m([b[x] for x in gains])-pm)) if gains else 0.0
    change=bm-am
    if abs(change-(losing+adding))>1e-10:
        raise AssertionError("Original source CWM loss+gain identity violated")
    return {"loss_component":losing,"gain_component":adding,
            "total_cwm_change":change,"before_n":len(a),"after_n":len(b),
            "n_lost":len(lost),"n_gained":len(gains),"n_persistent":len(pa),
            "mean_persistent":pm,
            "mean_lost_minus_persistent":(
                m([a[x] for x in lost])-pm if lost else None),
            "mean_gained_minus_persistent":(
                m([b[x] for x in gains])-pm if gains else None)}

def site_components(plots,names,eive,site,zone,metric,min_match):
    a=plots.get((site,2014,zone));b=plots.get((site,2017,zone))
    if a is None or b is None:return None,"original_source_missing_taxon_cell"
    a=vals(a,names,eive,metric)
    b=vals(b,names,eive,metric)
    if not a or not b or min(len(a),len(b))<min_match:
        return None,"minimum_EIVE_matched_species"
    if not set(a)&set(b):
        return None,"no_persistent_taxon_to_anchor_additive_partition"
    return components(a,b),"eligible"

def summary(values):
    if not values:return {"n":0,"mean":None,"median":None,"min":None,"max":None}
    return {"n":len(values),"mean":mean(values),"median":median(values),
            "min":min(values),"max":max(values)}

def groupview(records):
    fields=("loss_component","gain_component","total_cwm_change",
            "before_n","after_n","n_lost","n_gained","n_persistent",
            "mean_lost_minus_persistent","mean_gained_minus_persistent")
    out={key:summary([r[key] for r in records if r[key] is not None])
         for key in fields}
    out["count_net_warmer"]=sum(r["total_cwm_change"]>0 for r in records)
    out["count_loss_component_positive"]=sum(r["loss_component"]>0 for r in records)
    out["count_gain_component_positive"]=sum(r["gain_component"]>0 for r in records)
    return out

def zonal(plots,names,eive,metric,min_match):
    sites=cohort(plots)
    details={}
    output={}
    for z in ZONES:
        q=Counter()
        bysite={}
        for site in sites:
            r,reason=site_components(plots,names,eive,site,z,metric,min_match)
            q[reason]+=1
            if r is not None:
                bysite[site]=r
        details[z]=bysite
        output[z]={
            "source_zones":z,"metric":metric,
            "original_source_sites":len(sites),
            "status":dict(q),
            "by_stratum":{
                g:groupview([r for site,r in bysite.items()
                 if g=="all" or ("between" if site.startswith("Between ") else "station")==g])
                for g in ("all","station","between")}}
    return output,details

def paired_zone_comparison(details,reference):
    a=details["Track"]
    b=details[reference]
    shared=sorted(set(a)&set(b))
    effects={}
    for site in shared:
        left,right=a[site],b[site]
        effects[site]={key:left[key]-right[key] for key in (
            "loss_component","gain_component","total_cwm_change")}
        if abs(effects[site]["total_cwm_change"]-
               (effects[site]["loss_component"]+effects[site]["gain_component"]))>1e-10:
            raise AssertionError("Between-zone difference of taxon temporal components not additive")
    out={"reference":reference,"n_same_source_sites":len(shared)}
    for group in ("all","station","between"):
        rows=[v for s,v in effects.items() if group=="all" or
             ("between" if s.startswith("Between ") else "station")==group]
        out[group]={"n":len(rows),
            "track_minus_interior_2017_minus_2014_cwm":summary(
                [r["total_cwm_change"] for r in rows]),
            "track_minus_interior_loss_component":summary(
                [r["loss_component"] for r in rows]),
            "track_minus_interior_gain_component":summary(
                [r["gain_component"] for r in rows]),
            "negative_overall_difference_sites":sum(
                r["total_cwm_change"]<0 for r in rows)}
    return out

def analyze(kashmir_source,eive_source):
    plots,source_qc=parse_xlsx_matrix(Path(kashmir_source))
    names=species_names_from_publisher(kashmir_source)
    eive,eiveqc=eive_table(eive_source)
    out={"version":"v329","question":"additive observed loss vs gain components of recorded Kashmir railway and local interior EIVE",
         "original_kashmir":source_qc,"original_EIVE":eiveqc,
         "minimum_contributing_EIVE_systems":SYSTEMS,
         "metrics":{}}
    for metric in ("position","width"):
        out["metrics"][metric]={}
        for min_n in MIN_RECORDED_MATCHED:
            zones,details=zonal(plots,names,eive,metric,min_n)
            out["metrics"][metric][str(min_n)]={
                "minimum_matched_EIVE_species_each_year_each_zone":min_n,
                "original_zone_components":zones,
                "same_source_site_track_vs_interior_components":[
                    paired_zone_comparison(details,r)
                    for r in ("Interior Plot 1","Interior Plot 2")]}
    out["limitations"]=[
      "Post hoc exploratory diagnosis based on discovered v327-v328 temporal signal, not a newly prespecified confirmatory hypothesis",
      "Only site-years with at least one persistent EIVE species are decomposable with this anchor; full replacement sites are excluded and reported",
      "Loss/gain means recorded taxon absence/presence between 2014 and 2017, not demographic extinction or colonization",
      "No guarantee that source plots or habitats have identical environmental conditions between years",
      "Even n>=3 species in sparse 2017 Track creates noisy CWM; n>=1 panel descriptive and not evidence of persistent environmental filtering",
      "Both references within original railway corridor; cannot estimate counterfactual unaffected regional plant community",
      "EIVE is European-derived consensus trait proxy, not thermal niche physiological width measured in Kashmir",
      "Within single corridor 31 dependent sites: no population-scale field replication p-values or warming-climate causation",
      "No source original species-by-zone rows, GPS or site-level EIVE data published"]
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=analyze(a.kashmir_original,a.eive_original)
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf8")
    compact={"version":"v329","metrics":{}}
    for k,panels in result["metrics"].items():
        compact["metrics"][k]={}
        for cond,v in panels.items():
            compact["metrics"][k][cond]={
               "zone_summary":{z:{
                   "source_status":x["status"],
                   "all":x["by_stratum"]["all"]
               } for z,x in v["original_zone_components"].items()},
               "paired":v["same_source_site_track_vs_interior_components"]}
    print(json.dumps(compact,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
