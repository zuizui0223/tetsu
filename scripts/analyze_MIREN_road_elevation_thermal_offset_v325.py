#!/usr/bin/env python3
"""v325: Does MIREN ROAD-vs-interior thermal affinity change with elevation?

This analyzes an observational ROAD edge, not railway climate release.
Elevation is an external macro-environmental proxy, NOT a measured
temperature of ballast or a random exposure. Region/road/site/year
matching, source hashes and first-visit constraints are mandatory.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

from analyze_original_MIREN_EIVE_road_vs_interior_v323 import (
    read_miren,select_first_complete,MIREN_MD5,EUROPE_SOURCE_REGIONS)
from analyze_MIREN_EIVE_genus_position_width_v324 import eligible_roads
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Slovenia_thermal_generalist_section_v319 import section_test, section_drop

MIN_ROAD_SITES=4
MIN_ELEV_RANGE_M=50
MIN_ROADS_PER_REGION=2

def float_elev(x):
    try:v=float(x)
    except (ValueError,TypeError):return None
    return v if math.isfinite(v) and -500<v<9000 else None

def altitude_from_source(path):
    if hashlib.md5(Path(path).read_bytes()).hexdigest()!=MIREN_MD5:
        raise ValueError("MIREN original source MD5 mismatch")
    altitude=defaultdict(lambda:defaultdict(set))
    records=0
    with Path(path).open("r",newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f)
        for name in ("Region","Road","Site","Year","Plot","Plot_id","Elevation"):
            if name not in (reader.fieldnames or []):raise ValueError("No source elevation field")
        for row in reader:
            records+=1
            reg=(row["Region"] or "").strip()
            plot=(row["Plot"] or "").strip()
            if reg not in EUROPE_SOURCE_REGIONS or plot not in ("1","3"):continue
            if not (row["Plot_id"] or "").strip().endswith("."+plot):continue
            road=(row["Road"] or "").strip();site=(row["Site"] or "").strip()
            year=(row["Year"] or "").strip()
            if not road or not site or not year.isdigit():continue
            elev=float_elev(row.get("Elevation"))
            if elev is None:continue
            altitude[(reg,road,site,int(year))][plot].add(elev)
    if records!=167804:raise ValueError("Original 167804 MIREN rows drift")
    out={}
    qc={"siteyear_with_both_strict_plot_elevations":0,
        "siteyear_discordance_over_10m":0,
        "siteyear_discordance_over_50m_excluded":0,
        "siteyear_same_plot_internal_discordance_over_5m_excluded":0}
    for key,plotvalues in altitude.items():
        if not all(p in plotvalues for p in ("1","3")):continue
        qc["siteyear_with_both_strict_plot_elevations"]+=1
        if any(max(v)-min(v)>5 for v in plotvalues.values()):
            qc["siteyear_same_plot_internal_discordance_over_5m_excluded"]+=1
            continue
        values=[median(v) for p,v in sorted(plotvalues.items()) if p in ("1","3")]
        spread=abs(values[0]-values[1])
        if spread>10:qc["siteyear_discordance_over_10m"]+=1
        if spread>50:
            qc["siteyear_discordance_over_50m_excluded"]+=1
            continue
        out[key]=sum(values)/len(values)
    return out,qc

def mean(vals):
    return sum(vals)/len(vals) if vals else None

def road_slope(xs,ys):
    if len(xs)!=len(ys) or len(xs)<MIN_ROAD_SITES:
        raise ValueError("Not enough valid road sites")
    mx=mean(xs);my=mean(ys)
    den=sum((x-mx)**2 for x in xs)
    return (100*sum((x-mx)*(y-my) for x,y in zip(xs,ys))/den
            if den>0 else None)

def build_cohort(roads,altitudes):
    out=defaultdict(list)
    qa={"selected_first_survey_road_sites_before_elevation":0,
        "selected_sites_without_same_year_plot_elevation":0,
        "selected_sites_requiring_site_elevation_imputation":0}
    for (reg,road),sites in roads.items():
        for s in sites:
            qa["selected_first_survey_road_sites_before_elevation"]+=1
            key=(reg,road,s["site"],s["year"])
            if key not in altitudes:
                qa["selected_sites_without_same_year_plot_elevation"]+=1
                continue
            out[(reg,road)].append(dict(s,elevation_m=altitudes[key]))
    eligible={}
    for k,sites in out.items():
        if len(sites)<MIN_ROAD_SITES:continue
        elev=[s["elevation_m"] for s in sites]
        if max(elev)-min(elev)<MIN_ELEV_RANGE_M:continue
        eligible[k]=sites
    byreg=defaultdict(int)
    for reg,road in eligible:byreg[reg]+=1
    eligible={k:v for k,v in eligible.items()
              if byreg[k[0]]>=MIN_ROADS_PER_REGION}
    return eligible,qa

def score(codes,eive,metric):
    # Cohort prevalidated by v324: EIVE >=3 and full source plant lists.
    return mean([eive[s][metric] for s in codes])

def model(roads,eive,metric):
    byreg=defaultdict(list)
    roadresults={}
    for (reg,road),sites in sorted(roads.items()):
        elev=[s["elevation_m"] for s in sites]
        delta=[score(s["roadside"],eive,metric)-score(s["inland"],eive,metric)
               for s in sites]
        slope=road_slope(elev,delta)
        if slope is None:continue
        byreg[reg].append(slope)
        roadresults[(reg,road)]=slope
    # Do not assign multiple site/year rows independent region replication.
    used={k:v for k,v in byreg.items()
          if len(v)>=MIN_ROADS_PER_REGION}
    reg_effects={reg:mean(v) for reg,v in used.items()}
    return {"metric":metric,"min_road_elevation_range_m":MIN_ELEV_RANGE_M,
      "eligible_original_road_ids":sum(len(v) for v in used.values()),
      "eligible_original_region_codes":len(used),
      "eligible_site_count":sum(len(sites) for (reg,road),sites in roads.items()
                                 if reg in used),
      "region_equal_mean_delta_EIVE_per_100m_altitude":mean(list(reg_effects.values())),
      "positive_region_slopes":sum(x>0 for x in reg_effects.values()),
      "negative_region_slopes":sum(x<0 for x in reg_effects.values()),
      "region_results":[{"region":reg,"roads":len(used[reg]),
             "slope_per_100m":x,"positive_roads":sum(v>0 for v in used[reg]),
             "negative_roads":sum(v<0 for v in used[reg])}
             for reg,x in sorted(reg_effects.items())],
      "road_slope_signflip_exploratory":section_test(used),
      "leave_one_region_out":section_drop(used)}

def analyze(mirenpath,eivepath):
    plots,sourceqa=read_miren(mirenpath)
    altitude,altqa=altitude_from_source(mirenpath)
    eive,eiveqa=eive_table(eivepath)
    panel=[]
    for minimum_year in (2010,2017):
        candidate=select_first_complete(plots,minimum_year)
        source_roads=eligible_roads(candidate,eive)
        roads,gate=build_cohort(source_roads,altitude)
        for metric in ("position","width"):
            one=model(roads,eive,metric)
            one["minimum_first_original_survey_year"]=minimum_year
            one["source_elevation_qa"]=gate
            panel.append(one)
    return {"version":"v325","source":"MIREN European mountain ROAD Plot1 vs Plot3 paired elevation gradients",
      "MIREN_original_md5":MIREN_MD5,"EIVE_original_md5":eiveqa["EIVE_original_md5"],
      "source_miren_qa":sourceqa,"source_elevation_qa":altqa,
      "panels":panel,
      "warnings":[
       "Not a railway study or direct test of railway climate-offset hypothesis A",
       "Elevation is proxy for broad physical climate, not measured local track or roadside temperature",
       "Roadside versus inland plot differs in canopy/light/source pool as well as local abiotic stress",
       "Slope is descriptive within observed nonrandom mountain roads; no randomized climate treatment",
       "One MIREN source program with four candidate European region codes; regional sign flips are exploratory",
       "Contributing botanical plot-level original protocol may differ between MIREN road regions",
       "Road elevation means are only calculated when BOTH real Plot1 and Plot3 source elevations exist same year",
       "Contiguous sites along a road are spatially correlated: within-road permuted p values should not be taken as causal",
       "Ecological temperature indicator scale is not Celsius; outcome is rail-unrelated species preference"]}

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--miren-csv",required=True)
    a.add_argument("--eive-xlsx",required=True)
    a.add_argument("--output",required=True)
    p=a.parse_args()
    out=analyze(p.miren_csv,p.eive_xlsx)
    path=Path(p.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":"v325","source_elevation_qa":out["source_elevation_qa"],
      "panels":out["panels"]},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
