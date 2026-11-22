#!/usr/bin/env python3
"""v332: station vs between thermal contrast after true source geospatial matching.

Predeclared design: docs/tetsu_v332_Kashmir_station_between_geographic_matching_frozen_design_ja.md
Dryad authentic 31 Kashmir railway sites and EIVE 1.0 temperature preferences,
both checksum verified by original readers/workflow. Straight geodesic distance,
not railway network distance. No original site locations or species rows exported.
"""
import argparse
import json
import math
import re
from collections import defaultdict,Counter
from itertools import product
from pathlib import Path
from statistics import mean,median
from openpyxl import load_workbook
import numpy as np
from scipy.optimize import linear_sum_assignment

from analyze_Kashmir_same_sites_reference_year_v254 import DRYAD_SHA256,parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import (
    species_names_from_publisher,site_effect,cohort)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table

# ALL matching thresholds were fixed *before* this source-matched analysis.
MATCH_GATES={
  "G0_geo_nearest_no_caliper":(math.inf,math.inf,None),
  "G1_geographic_altitude_primary":(15.0,75.0,None),
  "G2_tight_geographic_altitude":(8.0,30.0,None),
  "G3_geo_altitude_richness_coverage":(15.0,75.0,(3,6,0.20)),
}
TRAIT_GATES={
  "strict_3systems_5species_50pct":(3,5,0.50),
  "relaxed_1system_3species_30pct":(1,3,0.30),
  "sparse_3systems_1species":(3,1,0.0)
}
ZONES=("Interior Plot 2","Interior Plot 1")
METRICS=("position","width")
YEARS=(2014,2017)

def parse_dms(value,is_lat):
    txt=str(value or "").strip()
    if not txt: raise ValueError("Original coordinate is empty")
    norm=txt.lower().replace("º","o").replace("°","o").replace("˚","o")
    norm=norm.replace("’","'").replace("′","'").replace("‘","'").replace("＇","'")
    chunks=re.findall(r"\d+(?:\.\d+)?",norm)
    if len(chunks)!=3 or "o" not in norm or "'" not in norm:
        raise ValueError("Unrecognized publisher original coordinate DMS shape")
    deg,minute,second=map(float,chunks)
    if minute<0 or minute>=60 or second<0 or second>=60 or deg<0:
        raise ValueError("Invalid source degree/minute/second")
    val=deg+minute/60+second/3600
    lower,upper=(30,38) if is_lat else (70,82)
    if not lower<=val<=upper:
        raise ValueError("Source geographic location outside broad original Kashmir area")
    return val

def site_metadata(path):
    from hashlib import sha256
    if sha256(Path(path).read_bytes()).hexdigest()!=DRYAD_SHA256:
        raise ValueError("Original Dryad source checksum mismatch")
    book=load_workbook(path,read_only=True,data_only=True)
    metadata={}
    seen=Counter()
    try:
        rows=book.active.values
        h=next(rows)
        if tuple(str(x or "").strip().lower() for x in h[:6])!=(
           "name of the site","latitude","longitude","altitude",
           "year of sampling","location"):
            raise ValueError("Original publisher source geography column headers changed")
        for row in rows:
            name=" ".join(str(row[0] or "").split())
            if not name: raise ValueError("Missing source site ID")
            lat=parse_dms(row[1],True);lon=parse_dms(row[2],False)
            alt=float(row[3])
            if not 500<alt<3000 or not math.isfinite(alt):
                raise ValueError("Original source altitude invalid")
            yr=int(row[4]);zone=" ".join(str(row[5] or "").split())
            if yr not in YEARS or zone not in ("Track","Embankment",*ZONES):
                raise ValueError("Source year or plot zone invalid")
            tup=(lat,lon,alt)
            if name in metadata and any(abs(x-y)>1e-8
                   for x,y in zip(metadata[name]["raw_geo"],tup)):
                raise ValueError("Coordinates or altitude differ among same source site plot/year rows")
            metadata.setdefault(name,{
               "raw_geo":tup,
               "station_type":"between" if name.startswith("Between ") else "station"})
            seen[name]+=1
    finally:
        book.close()
    if len(metadata)!=31 or sum(d["station_type"]=="between" for d in metadata.values())!=15:
        raise ValueError("Original 31-site station-between universe changed")
    if any(n!=8 for n in seen.values()):
        raise ValueError("Original Kashmir site 2yrs x 4source plots missing")
    return metadata,{
      "n_field_sites":len(metadata),"n_original_source_row_records":sum(seen.values()),
      "n_station":16,"n_between":15,
      "lat_range":[min(x["raw_geo"][0] for x in metadata.values()),
                   max(x["raw_geo"][0] for x in metadata.values())],
      "lon_range":[min(x["raw_geo"][1] for x in metadata.values()),
                   max(x["raw_geo"][1] for x in metadata.values())],
      "alt_range":[min(x["raw_geo"][2] for x in metadata.values()),
                   max(x["raw_geo"][2] for x in metadata.values())],
      "coordinates_outside_aggregate_not_republished":True}

def haversine_km(x,y):
    phi1,lam1,phi2,lam2=map(math.radians,(x[0],x[1],y[0],y[1]))
    dp=phi2-phi1;dl=lam2-lam1
    q=math.sin(dp/2)**2+math.cos(phi1)*math.cos(phi2)*math.sin(dl/2)**2
    return 2*6371.0088*math.asin(min(1,math.sqrt(q)))

def site_panel(plots,names,eive,metadata,year,reference,metric,gate):
    systems,min_species,min_cover=TRAIT_GATES[gate]
    selected={}
    q=Counter()
    for site in sorted(metadata):
        effect,info=site_effect(plots,names,eive,site,year,
            ("Track",reference),metric,systems,min_species,min_cover)
        if effect is None:
            q[info["reason"]]+=1
            continue
        q["eligible"]+=1
        a,b=info["track"],info["reference"]
        selected[site]={
          "group":metadata[site]["station_type"],
          "geo":metadata[site]["raw_geo"],
          "temp_total":effect,
          "temp_track":a["cwm"],"temp_interior":b["cwm"],
          "track_richness":a["n_source"],"interior_richness":b["n_source"],
          "track_trait_n":a["n_eive"],"interior_trait_n":b["n_eive"],
          "track_cov":a["coverage"],"interior_cov":b["coverage"]}
    return selected,dict(q)

def source_match(selected,gate):
    maximal_distance,maximal_alt,extra=MATCH_GATES[gate]
    between=sorted((k for k,v in selected.items() if v["group"]=="between"))
    stations=sorted((k for k,v in selected.items() if v["group"]=="station"))
    n=len(between);m=len(stations)
    if not n or not m:
        return [],{"eligible_between_sites":n,"eligible_station_sites":m,
                   "feasible_geographic_edges":0,"matched_count":0}
    # Append one dummy column per between site. Maximizes number of
    # feasible unique stations before minimizing distance/elevation cost.
    matrix=np.full((n,m+n),1e6,dtype=float)
    matrix[:,m:]=1e5
    feasibility=0
    for i,a in enumerate(between):
        A=selected[a]
        for j,b in enumerate(stations):
            B=selected[b]
            dist=haversine_km(A["geo"],B["geo"])
            alt=abs(A["geo"][2]-B["geo"][2])
            if dist>maximal_distance or alt>maximal_alt:
                continue
            if extra is not None:
                track_diff,int_diff,covdiff=extra
                if (abs(A["track_richness"]-B["track_richness"])>track_diff or
                    abs(A["interior_richness"]-B["interior_richness"])>int_diff or
                    abs(A["track_cov"]-B["track_cov"])>covdiff or
                    abs(A["interior_cov"]-B["interior_cov"])>covdiff):
                    continue
            feasibility+=1
            matrix[i,j]=dist+0.03*alt+1e-6*(i*m+j)
    rid,cid=linear_sum_assignment(matrix)
    matched=[]
    for i,j in zip(rid,cid):
        if j>=m or matrix[i,j]>=1e5:
            continue
        a,b=selected[between[i]],selected[stations[j]]
        matched.append({"between":a,"station":b,
            "geodesic_km":haversine_km(a["geo"],b["geo"]),
            "alt_abs_m":abs(a["geo"][2]-b["geo"][2]),
            "alt_signed_between_minus_station_m":a["geo"][2]-b["geo"][2]})
    return matched,{"eligible_between_sites":n,"eligible_station_sites":m,
        "feasible_geographic_edges":feasibility,"matched_count":len(matched),
        "unique_station_reuse_prevented":True}

def mean_safe(v):
    return mean(v) if v else None

def quant(v):
    if not v:return {"n":0,"mean":None,"median":None,"max":None}
    return {"n":len(v),"mean":mean(v),"median":median(v),"max":max(v)}

def exact_signflip(vals):
    if not 3<=len(vals)<=16:return None
    obs=abs(mean(vals))
    n=len(vals);total=1<<n
    c=0
    for mask in range(total):
        stat=abs(sum(x if (mask>>i)&1 else -x
                 for i,x in enumerate(vals))/n)
        if stat>=obs-1e-12:c+=1
    return {"n_matched_pairs":n,
            "two_sided_exact_signflip_fraction":c/total,
            "exact_label_or_pair_exchangeability_NOT_randomized_field_treatment":True}

def summarize_pairs(pairs,selection,match_gate,metric,year,reference,trait_gate,source_q):
    gap=[p["between"]["temp_total"]-p["station"]["temp_total"] for p in pairs]
    track=[p["between"]["temp_track"]-p["station"]["temp_track"] for p in pairs]
    interior=[p["between"]["temp_interior"]-p["station"]["temp_interior"] for p in pairs]
    if any(abs(x-y+z)>1e-10 for x,y,z in zip(gap,track,interior)):
        raise AssertionError("Matched between-station difference must decompose track minus interior")
    out={
      "year":year,"reference":reference,"trait_metric":metric,
      "trait_support_gate":trait_gate,"match_gate":match_gate,
      "original_eligibility":source_q,
      "source_match_qa":selection,
      "matched_between_minus_station":{
        "delta_track_minus_interior":quant(gap),
        "track_CWM_only":quant(track),
        "interior_CWM_only":quant(interior),
        "between_mean_track_minus_interior":mean_safe(
             [p["between"]["temp_total"] for p in pairs]),
        "station_mean_track_minus_interior":mean_safe(
             [p["station"]["temp_total"] for p in pairs]),
        "positive_paired_differences":sum(x>0 for x in gap),
        "negative_paired_differences":sum(x<0 for x in gap),
        "pair_label_sign_flip_diagnostic":exact_signflip(gap)},
      "balance":{
        "absolute_geodesic_distance_km":quant([p["geodesic_km"] for p in pairs]),
        "absolute_station_between_altitude_difference_m":quant([p["alt_abs_m"] for p in pairs]),
        "signed_altitude_between_minus_station_m":quant([p["alt_signed_between_minus_station_m"] for p in pairs]),
        "track_source_richness_between_minus_station":quant([
            p["between"]["track_richness"]-p["station"]["track_richness"] for p in pairs]),
        "interior_source_richness_between_minus_station":quant([
            p["between"]["interior_richness"]-p["station"]["interior_richness"] for p in pairs]),
        "track_trait_coverage_between_minus_station":quant([
            p["between"]["track_cov"]-p["station"]["track_cov"] for p in pairs]),
        "interior_trait_coverage_between_minus_station":quant([
            p["between"]["interior_cov"]-p["station"]["interior_cov"] for p in pairs]),
        "track_trait_count_between_minus_station":quant([
            p["between"]["track_trait_n"]-p["station"]["track_trait_n"] for p in pairs])},
      "limitations":[
       "Geographic distances are straight great-circle distances, not physical rail-route distances",
       "The station/between classification was not randomized and remains confounded with settlement and management",
       "Pair sign flips are not true randomization causal p-values",
       "Strict Kashmir 2014 between support is very limited",
       "EIVE is European ecological indicators transferred to Asian plant occurrence lists",
       "No original site-specific coordinates or source biological names in aggregate"]}
    return out

def analyze(original_path,eive_path):
    plots,plotqa=parse_xlsx_matrix(Path(original_path))
    names=species_names_from_publisher(original_path)
    eive,eiveqa=eive_table(eive_path)
    meta,metaqa=site_metadata(original_path)
    if set(meta)!=set(cohort(plots)):
        raise ValueError("Geographic and species source sites differ")
    out={"version":"v332","scope":"WITHIN ONE KASHMIR RAILWAY CORRIDOR STATION VERSUS BETWEEN; OBSERVATIONAL",
      "source_original_Dryad":"10.5061/dryad.0p2ngf20h",
      "original_sha256":DRYAD_SHA256,
      "independent_EIVE":"10.5281/zenodo.7534792",
      "eive_source_md5":eiveqa["EIVE_original_md5"],
      "original_plot_QA":plotqa,"aggregate_geographic_source_QA":metaqa,
      "matched_panels":[]}
    for year,reference,metric,trait_gate in product(
            YEARS,ZONES,METRICS,TRAIT_GATES):
        valid,src_q=site_panel(plots,names,eive,meta,year,reference,metric,trait_gate)
        for match_gate in MATCH_GATES:
            pairs,matchqa=source_match(valid,match_gate)
            out["matched_panels"].append(summarize_pairs(
                pairs,matchqa,match_gate,metric,year,reference,trait_gate,src_q))
    out["limitations"]=[
     "Only 31 Kashmir source sites, no between-region railway-vs-nonrail climate interaction",
     "Species EIVE-centre and niche width proxies do not measure actual temperature in degrees C",
     "2014 strict eligibility leaves only five between-station sites",
     "Original coordinates are DMS strings converted without any missing-value location imputation",
     "Matching may eliminate most pairs and cannot adjust unmeasured land cover, openness, management or traffic",
     "Richness is potentially a mediator of habitat filtering; richness matching is secondary sensitivity",
     "Nearby interior plots may also experience railway-associated changes",
     "Paired signflip is an exploratory exchangeability diagnostic conditional on matching, not randomized experiment",
     "2017 strict comparison has zero eligible source sites; relaxed sparse-site tests remain nonconfirmatory",
     "Original field site coordinates or source taxon-by-plot table are not redistributed"]
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
    selected=[x for x in result["matched_panels"] if
        x["trait_metric"]=="position" and x["reference"]=="Interior Plot 2" and
        (x["year"]==2014 or x["trait_support_gate"]=="sparse_3systems_1species")]
    print(json.dumps({"version":"v332","metadata":result["aggregate_geographic_source_QA"],
       "primary_and_sensitivity":[{"year":x["year"],"gate":x["trait_support_gate"],
        "match":x["match_gate"],"support":x["source_match_qa"],
        "station_between_CWM_effect":x["matched_between_minus_station"],
        "balance":x["balance"]} for x in selected]},
        ensure_ascii=False),flush=True)

if __name__=="__main__":main()
