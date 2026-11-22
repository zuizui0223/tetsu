#!/usr/bin/env python3
"""v326 source-authenticated MIREN Plot1-Plot3 vertical displacement audit.

Plot1 road edge and Plot3 52-102m inland can differ in ALTITUDE on a
mountain transect. Sensitivity restricts true same-site/year plots to
near-equal measured elevation before interpreting EIVE thermal contrasts.
ROAD data only, not a railway versus nonrail experiment.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median

from analyze_original_MIREN_EIVE_road_vs_interior_v323 import (
    MIREN_MD5,EUROPE_SOURCE_REGIONS,read_miren,select_first_complete)
from analyze_MIREN_EIVE_genus_position_width_v324 import eligible_roads
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Slovenia_thermal_generalist_section_v319 import section_test

TOLERANCES=(0,2,5,10,20,50)
MIN_SITES_PER_ROAD=4
MIN_ROADS_PER_REGION=2

def num(x):
    try:n=float(x)
    except (ValueError,TypeError):return None
    return n if math.isfinite(n) and -500<n<9000 else None

def pair_plot_altitudes(path):
    if hashlib.md5(Path(path).read_bytes()).hexdigest()!=MIREN_MD5:
        raise ValueError("Not the exact original MIREN source")
    original=defaultdict(lambda:defaultdict(set))
    n=0
    with Path(path).open("r",encoding="utf-8-sig",newline="") as f:
        rdr=csv.DictReader(f)
        for c in ("Region","Road","Site","Year","Plot","Plot_id","Elevation"):
            if c not in (rdr.fieldnames or []):raise ValueError("Missing source field "+c)
        for row in rdr:
            n+=1
            reg=str(row["Region"] or "").strip()
            plot=str(row["Plot"] or "").strip()
            if reg not in EUROPE_SOURCE_REGIONS or plot not in ("1","3"):continue
            year=str(row["Year"] or "").strip()
            road=str(row["Road"] or "").strip()
            site=str(row["Site"] or "").strip()
            if not (year.isdigit() and road and site):continue
            if not str(row["Plot_id"] or "").strip().endswith("."+plot):continue
            v=num(row["Elevation"])
            if v is None:continue
            original[(reg,road,site,int(year))][plot].add(v)
    if n!=167804:raise ValueError("Original MIREN species rows changed")
    out={}
    qc=Counter()
    for key,side in original.items():
        if "1" not in side or "3" not in side:continue
        qc["siteyears_with_two_actual_plot_elevations"]+=1
        if any(max(v)-min(v)>5 for v in side.values()):
            qc["same_plot_discordant_more_than_5m"]+=1
            continue
        road_alt=median(side["1"])
        inland_alt=median(side["3"])
        delta=road_alt-inland_alt
        if abs(delta)>10:qc["two_plots_separated_vertically_more_than_10m"]+=1
        if abs(delta)>50:qc["two_plots_separated_vertically_more_than_50m"]+=1
        out[key]={"roadside_m":road_alt,"interior_m":inland_alt,"delta_m":delta}
    return out,dict(qc)

def mean(xs):
    return sum(xs)/len(xs) if xs else None

def make_records(source_roads,elev,eive):
    rows=[]
    q=Counter()
    for (reg,road),sites in source_roads.items():
        for site in sites:
            q["original_thermal_evaluable_first_visit_sites"]+=1
            k=(reg,road,site["site"],site["year"])
            if k not in elev:
                q["missing_real_paired_source_elevation"]+=1;continue
            a=site["roadside"];b=site["inland"]
            if not a or not b:raise ValueError("Unexpected unobserved MIREN plot")
            rows.append({"region":reg,"road":road,"site":site["site"],
                  "source_year":site["year"],"plot1_minus_plot3_alt_m":elev[k]["delta_m"],
                  "position_delta":mean([eive[x]["position"] for x in a])-
                                   mean([eive[x]["position"] for x in b]),
                  "width_delta":mean([eive[x]["width"] for x in a])-
                                mean([eive[x]["width"] for x in b])})
    return rows,dict(q)

def summarize(rows,metric,tolerance):
    # Selection is by source *physical elevation* only, never by
    # observed plant EIVE outcome or its sign.
    kept=[r for r in rows if abs(r["plot1_minus_plot3_alt_m"])<=tolerance]
    byroad=defaultdict(list)
    for row in kept:
        byroad[(row["region"],row["road"])].append(row)
    byroad={k:v for k,v in byroad.items() if len(v)>=MIN_SITES_PER_ROAD}
    byreg=defaultdict(list)
    for (reg,road),samples in sorted(byroad.items()):
        byreg[reg].append((road,samples))
    byreg={k:v for k,v in byreg.items() if len(v)>=MIN_ROADS_PER_REGION}
    regionresults={}
    for reg,rr in sorted(byreg.items()):
        values=[mean([s[metric+"_delta"] for s in site]) for _,site in rr]
        regionresults[reg]={"roads":len(rr),
            "sites":sum(len(site) for _,site in rr),
            "road_equal_mean":mean(values),
            "positive_roads":sum(x>0 for x in values),
            "negative_roads":sum(x<0 for x in values)}
    return {"metric":metric,"maximum_absolute_same_site_plot_elevation_difference_m":tolerance,
       "all_pre_gate_matched_sites_within_altitude_band":len(kept),
       "supported_regions":len(regionresults),
       "supported_roads":sum(x["roads"] for x in regionresults.values()),
       "supported_sites":sum(x["sites"] for x in regionresults.values()),
       "region_equal_mean_roadside_minus_interior_EIVE":(
           mean([x["road_equal_mean"] for x in regionresults.values()])
           if regionresults else None),
       "positive_regions":sum(x["road_equal_mean"]>0 for x in regionresults.values()),
       "region_details":regionresults,
       "source_region_sign_flip_exploratory":section_test({
           r:[mean([s[metric+"_delta"] for s in site]) for _,site in rr]
           for r,rr in byreg.items()})}

def analyze(miren_path,eive_path):
    vegetation,srcqc=read_miren(miren_path)
    alt,altqc=pair_plot_altitudes(miren_path)
    eive,eivemeta=eive_table(eive_path)
    first=select_first_complete(vegetation,2010)
    eligible=eligible_roads(first,eive)
    records,q=make_records(eligible,alt,eive)
    deltas=[r["plot1_minus_plot3_alt_m"] for r in records]
    qa={"with_elevation":len(records),
        "plot1_lower_m":sum(x<0 for x in deltas),
        "plot1_higher_m":sum(x>0 for x in deltas),
        "equal_m":sum(x==0 for x in deltas),
        "mean_plot1_minus_plot3_elevation_m":mean(deltas),
        "median_plot1_minus_plot3_elevation_m":median(deltas) if deltas else None,
        "min":min(deltas) if deltas else None,
        "max":max(deltas) if deltas else None}
    panels=[summarize(records,metric,threshold)
           for metric in ("position","width") for threshold in TOLERANCES]
    return {"version":"v326","context":"MIREN ROADS-only physical plot elevation offset falsification",
       "MIREN_MD5":MIREN_MD5,"EIVE_MD5":eivemeta["EIVE_original_md5"],
       "source_site_elevation_integrity":altqc,
       "eligible_site_sample":q,
       "physical_elevation_delta_summary":qa,
       "panels":panels,
       "limits":[
           "The paired plots are both mountain-road survey plots, not railway-versus-natural plots",
           "The signed vertical position difference may be related to landform and plot canopy, not only climate",
           "Elevation-displacement filtering changes the sample and may cause strong loss of roads or regions; no inference when support lost",
           "EIVE temperature centre and width are ecological preference proxies, not sensor temperature or physiological tolerances",
           "Mountain road and railway ballast habitats differ in survey design; never call this a railway-specific mechanism",
           "No sample without a truly recorded plot elevation or plot species record is imputed"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--miren-csv",required=True)
    p.add_argument("--eive-xlsx",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    data=analyze(args.miren_csv,args.eive_xlsx)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
      "version":"v326",
      "source_site_elevation_integrity":data["source_site_elevation_integrity"],
      "physical_elevation_delta_summary":data["physical_elevation_delta_summary"],
      "panels":[{k:p[k] for k in (
           "metric","maximum_absolute_same_site_plot_elevation_difference_m",
           "all_pre_gate_matched_sites_within_altitude_band",
           "supported_regions","supported_roads","supported_sites",
           "region_equal_mean_roadside_minus_interior_EIVE",
           "positive_regions","region_details")} for p in data["panels"]]
      },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
