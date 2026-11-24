#!/usr/bin/env python3
"""v328: Source-authenticated Kashmir temporal thermal attrition and reference.

Replay Dryad 2014/2017 same field sites, with EIVE species preference
indicators. EIVE is Europe-calibrated, not direct physiology in Kashmir.
This study tests descriptive *recorded species turnover*; not individual
mortality, colonization, train seed carriage or causal climate change.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_Kashmir_same_sites_reference_year_v254 import (
    YEARS,ZONES,parse_xlsx_matrix)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import (
    species_names_from_publisher,cohort,site_effect,MAIN_ZONES,ALTERNATE_ZONES,
    profile
)

ALL_ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
GATES=((3,5,.5),(1,3,.3),(3,1,.0))
METRICS=("position","width")
def mean_safe(values):
    return mean(values) if values else None

def quant(values):
    if not values:return {"n":0,"mean":None,"median":None,"min":None,"max":None}
    return {"n":len(values),"mean":mean(values),"median":median(values),
            "min":min(values),"max":max(values)}

def eligible_temporal_site(plots,names,eive,site,zones,metric,systems,min_n,min_c):
    traits={}
    for year in YEARS:
        a,b=site_effect(plots,names,eive,site,year,zones,metric,systems,min_n,min_c)
        if a is None:return None
        traits[year]=b
    track_change=traits[2017]["track"]["cwm"]-traits[2014]["track"]["cwm"]
    interior_change=traits[2017]["reference"]["cwm"]-traits[2014]["reference"]["cwm"]
    if abs((track_change-interior_change)-
           ((traits[2017]["track"]["cwm"]-traits[2017]["reference"]["cwm"])-
            (traits[2014]["track"]["cwm"]-traits[2014]["reference"]["cwm"])))>1e-10:
        raise AssertionError("Within-source temporal difference-of-differences not additive")
    return {"delta_track":track_change,"delta_interior":interior_change,
            "delta_track_minus_interior":track_change-interior_change,
            "track_known_2014":traits[2014]["track"]["n_eive"],
            "track_known_2017":traits[2017]["track"]["n_eive"],
            "interior_known_2014":traits[2014]["reference"]["n_eive"],
            "interior_known_2017":traits[2017]["reference"]["n_eive"]}

def transition_details(plots,names,eive,site,zone,metric,systems):
    a=plots.get((site,2014,zone));b=plots.get((site,2017,zone))
    if a is None or b is None:return None
    x={names[i] for i in a if names[i] in eive and eive[names[i]]["systems"]>=systems}
    y={names[i] for i in b if names[i] in eive and eive[names[i]]["systems"]>=systems}
    if not x or not y:return None
    lost=x-y
    persist=x&y
    gained=y-x
    def val(sets):
        return mean_safe([eive[s][metric] for s in sets])
    return {"source_2014_total":len(a),"source_2017_total":len(b),
            "EIVE_2014_known":len(x),"EIVE_2017_known":len(y),
            "n_lost":len(lost),"n_persist":len(persist),"n_gained":len(gained),
            "mean_lost":val(lost),"mean_persist":val(persist),
            "mean_gained":val(gained),
            "lost_minus_persist":(val(lost)-val(persist)) if lost and persist else None,
            "gained_minus_persist":(val(gained)-val(persist)) if gained and persist else None,
            "delta_cwm_year":val(y)-val(x)}

def describe_temporal(plots,names,eive,systems):
    sites=cohort(plots)
    zone={}
    for z in ALL_ZONES:
        valid=[(site,transition_details(plots,names,eive,site,z,"position",systems))
               for site in sites]
        good=[(s,d) for s,d in valid if d is not None]
        bygroup={}
        for group in ("all","station","between"):
            take=[d for s,d in good if group=="all" or
                  ("between" if s.startswith("Between ") else "station")==group]
            bygroup[group]={
               "n_source_sites_with_both_nonempty_trait_profiles":len(take),
               "source_richness_2014":quant([d["source_2014_total"] for d in take]),
               "source_richness_2017":quant([d["source_2017_total"] for d in take]),
               "source_richness_year_change":quant([
                    d["source_2017_total"]-d["source_2014_total"] for d in take]),
               "EIVE_known_2014":quant([d["EIVE_2014_known"] for d in take]),
               "EIVE_known_2017":quant([d["EIVE_2017_known"] for d in take]),
               "EIVE_position_CWM_2017_minus_2014":quant([
                     d["delta_cwm_year"] for d in take]),
               "lost_taxon_counts":quant([d["n_lost"] for d in take]),
               "persisting_taxon_counts":quant([d["n_persist"] for d in take]),
               "gained_taxon_counts":quant([d["n_gained"] for d in take]),
               "lost_minus_persist_EIVE_position":quant([
                  d["lost_minus_persist"] for d in take
                  if d["lost_minus_persist"] is not None]),
               "gained_minus_persist_EIVE_position":quant([
                  d["gained_minus_persist"] for d in take
                  if d["gained_minus_persist"] is not None])}
        zone[z]=bygroup
    return zone

def paired_temporal(plots,names,eive):
    sites=cohort(plots)
    out=[]
    for zones in (MAIN_ZONES,ALTERNATE_ZONES):
      for systems,minn,cov in GATES:
       for metric in METRICS:
        cohortrows=[]
        for site in sites:
            r=eligible_temporal_site(plots,names,eive,site,zones,metric,systems,minn,cov)
            if r is None:continue
            cohortrows.append((site,r))
        groups={}
        for group in ("all","station","between"):
            sub=[x for site,x in cohortrows if group=="all" or
                 ("between" if site.startswith("Between ") else "station")==group]
            groups[group]={
             "n_same_source_sites":len(sub),
             "Track_EIVE_CWM_temporal_change":quant([x["delta_track"] for x in sub]),
             "Interior_EIVE_CWM_temporal_change":quant([x["delta_interior"] for x in sub]),
             "difference_in_differences":quant([x["delta_track_minus_interior"] for x in sub]),
             "Track_EIVE_matched_species_2014":quant([x["track_known_2014"] for x in sub]),
             "Track_EIVE_matched_species_2017":quant([x["track_known_2017"] for x in sub]),
             "Interior_EIVE_matched_species_2014":quant([x["interior_known_2014"] for x in sub]),
             "Interior_EIVE_matched_species_2017":quant([x["interior_known_2017"] for x in sub]),
             "sites_with_negative_thermal_difference_change":sum(
                    x["delta_track_minus_interior"]<0 for x in sub)}
        out.append({"zone_pair":zones,"metric":metric,"source_systems":systems,
           "min_numeric_taxa_each_zone":minn,"min_trait_coverage_each_zone":cov,
           "matched_reference_is_within_same_railway_corridor":True,
           "paired_field_site_year_changes":groups})
    return out

def analyze(kashmir_xlsx,eive_xlsx):
    plots,qa=parse_xlsx_matrix(Path(kashmir_xlsx))
    names=species_names_from_publisher(kashmir_xlsx)
    traits,eiveqa=eive_table(eive_xlsx)
    out={"version":"v328","source":"Kashmir genuine Dryad paired 2014 and 2017 railway field plots",
         "kashmir_qa":qa,"eive_qa":eiveqa,
         "transition_trait_matched_systems_ge3":describe_temporal(plots,names,traits,3),
         "temporal_pair_coverage_sensitivity":paired_temporal(plots,names,traits),
         "warnings":[
           "These are observed taxon presence losses/gains between sampling years, not individual extinction, dispersal or survival probabilities",
           "All sites lie on one Kashmir railway corridor, observational and spatially correlated",
           "Reference Interior Plot 1/2 lies within 150m of railway, not proved environmentally unaffected",
           "Missing original taxon cell makes source plot incomplete, never zero",
           "EIVE Europe-calibrated temperature preferences are not physiological species temperature tolerances measured in Asia",
           "EIVE n>=3 five-species gate is original 2014 strict primary; a relaxed 2017 temporal comparison is exploratory and not an independent validation",
           "Different interior references change temporal signal: not generalizable until reference controls converge",
           "Same source sites are observed twice: no two independent country comparisons",
           "No direct local measured temperature, fitness, rail seed carriage or experimentally assigned disturbance"]}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=analyze(a.kashmir_original,a.eive_original)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":"v328",
         "transition_trait_matched_systems_ge3":result["transition_trait_matched_systems_ge3"],
         "temporal_pair_coverage_sensitivity":result["temporal_pair_coverage_sensitivity"]
         },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
