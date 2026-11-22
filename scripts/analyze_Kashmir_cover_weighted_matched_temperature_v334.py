#!/usr/bin/env python3
"""v334: verify authentic Dryad species positive source values and weight EIVE CWM.

The published source abundance/cover numeric units are preserved AS IS; no
assumption that source numbers are percentages or biomass. Compare equal
occurrence weight, raw positive source weight, and sqrt source weight while
holding the same v332 source-geo station-between pairs fixed.
"""
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path
from statistics import mean,median
from openpyxl import load_workbook

from analyze_Kashmir_same_sites_reference_year_v254 import DRYAD_SHA256,parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_geographic_station_between_EIVE_v332 import (
    site_metadata,site_panel,source_match,TRAIT_GATES,exact_signflip
)

PANELS=(
  (2014,"relaxed_1system_3species_30pct","G1_geographic_altitude_primary"),
  (2014,"strict_3systems_5species_50pct","G1_geographic_altitude_primary"),
  (2014,"relaxed_1system_3species_30pct","G3_geo_altitude_richness_coverage"),
  (2017,"sparse_3systems_1species","G1_geographic_altitude_primary"),
)
REFERENCE="Interior Plot 2"
WEIGHTINGS=("presence","original_positive_source_value","sqrt_original_source_value")

def read_original_plant_values(path,original_presence):
    if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=DRYAD_SHA256:
        raise ValueError("Original Dryad checksum changed")
    book=load_workbook(path,read_only=True,data_only=True)
    qa=Counter(); data={}; positive_values=[]
    try:
        rows=book.active.values
        header=next(rows)
        if len(header)!=129:raise ValueError("Original source header width changed")
        for row in rows:
            site=" ".join(str(row[0] or "").split())
            year=int(row[4]); zone=" ".join(str(row[5] or "").split())
            key=(site,year,zone)
            if key in data:raise ValueError("Source species-by-site-year-zone repeated")
            recorded={}
            any_missing=False
            if len(row)!=129:raise ValueError("Unexpected original numeric row length")
            for i,val in enumerate(row[6:]):
                qa["source_taxon_cell_fields"]+=1
                if val is None:
                    qa["source_missing_taxon_cells"]+=1
                    any_missing=True
                    continue
                if isinstance(val,bool) or not isinstance(val,(float,int)) or not math.isfinite(val) or val<0:
                    raise ValueError("Original published species numerical record non-numeric or negative")
                if val>0:
                    recorded[i]=float(val)
                    positive_values.append(float(val))
                    qa["source_positive_taxon_cells"]+=1
                    if val>100:qa["source_positive_cells_above100"]=qa.get("source_positive_cells_above100",0)+1
            data[key]=None if any_missing or not recorded else recorded
            source=original_presence.get(key)
            if source is None and data[key] is not None:
                raise ValueError("Original v254 presence identifies plot missing but cover parser did not")
            if source is not None and (data[key] is None or set(data[key])!=source):
                raise ValueError("Original v254 positive presence and cover value identities differ")
    finally:book.close()
    if len(data)!=248 or qa["source_missing_taxon_cells"]!=4:
        raise ValueError("Original Kashmir 248 rows and four missing cells must hold")
    if not positive_values:raise ValueError("Original source lacked positive field plant values")
    return data,{
        "source_plot_rows":len(data),
        **dict(qa),
        "original_positive_source_values_n":len(positive_values),
        "minimum_original_positive_source_value":min(positive_values),
        "median_original_positive_source_value":median(positive_values),
        "maximum_original_positive_source_value":max(positive_values),
        "no_assumption_source_numbers_equal_percent_cover":True}

def weighted_community(record,names,eive,systems,metric,scheme):
    if record is None:return None
    chosen={}
    for species,value in record.items():
        name=names[species]
        if name in eive and eive[name]["systems"]>=systems:
            chosen[name]=value
    if not chosen:return None
    weights={}
    for name,source_value in chosen.items():
        if scheme=="presence":w=1.0
        elif scheme=="original_positive_source_value":w=source_value
        elif scheme=="sqrt_original_source_value":w=math.sqrt(source_value)
        else:raise ValueError("Not a frozen source weighting method")
        if not w>0:raise ValueError("No original positive weight")
        weights[name]=w
    total=sum(weights.values())
    p={s:w/total for s,w in weights.items()}
    return {
        "cwm":sum(p[s]*eive[s][metric] for s in p),
        "taxon_n":len(chosen),
        "effective_trait_species_number":1/sum(v*v for v in p.values()),
        "largest_single_species_weight_fraction":max(p.values()),
        "total_source_positive_weight_in_matched_taxa":sum(chosen.values())}

def group_comparison(values,sites,names,eive,year,traitgate,scheme,metric):
    systems=TRAIT_GATES[traitgate][0]
    paired=[]
    for pair in sites:
        b,st=pair["between_source_site"],pair["station_source_site"]
        group={}
        for key,sitename,zone in (
            ("BT",b,"Track"),("BI",b,REFERENCE),
            ("ST",st,"Track"),("SI",st,REFERENCE)):
            group[key]=weighted_community(values.get((sitename,year,zone)),
                        names,eive,systems,metric,scheme)
            if group[key] is None:
                raise ValueError("Matched source plot lacks cover value for qualifying EIVE taxa")
        paired.append(group)
    deltas=[p["BT"]["cwm"]-p["BI"]["cwm"]-
             p["ST"]["cwm"]+p["SI"]["cwm"] for p in paired]
    track=[p["BT"]["cwm"]-p["ST"]["cwm"] for p in paired]
    interior=[p["BI"]["cwm"]-p["SI"]["cwm"] for p in paired]
    if any(abs(a-b+c)>1e-10 for a,b,c in zip(deltas,track,interior)):
        raise AssertionError("Weight-aware station-between decomposition must be additive")
    out={"weighting_scheme":scheme,"matched_pairs":len(paired),
        "mean_between_minus_station_track_relative_interior":mean(deltas),
        "positive_matched_pairs":sum(x>0 for x in deltas),
        "negative_matched_pairs":sum(x<0 for x in deltas),
        "mean_between_minus_station_track_only":mean(track),
        "mean_between_minus_station_interior_only":mean(interior),
        "conditional_matched_signflip_diagnostic":exact_signflip(deltas),
        "known_trait_effective_species_number":{
           key:{"mean":mean(p[key]["effective_trait_species_number"] for p in paired),
             "min":min(p[key]["effective_trait_species_number"] for p in paired),
             "max":max(p[key]["effective_trait_species_number"] for p in paired)}
           for key in ("BT","BI","ST","SI")},
        "source_largest_species_weight_fraction":{
           key:mean(p[key]["largest_single_species_weight_fraction"] for p in paired)
           for key in ("BT","BI","ST","SI")}}
    return out

def cohort(plots,names,eive,geo,year,traitgate,matchgate):
    from analyze_Kashmir_geographic_station_between_EIVE_v332 import site_panel
    selected,qa=site_panel(plots,names,eive,geo,year,REFERENCE,"position",traitgate)
    pairs,matchqa=source_match(selected,matchgate)
    reverse={id(v):k for k,v in selected.items()}
    matches=[{"between_source_site":reverse[id(p["between"])],
              "station_source_site":reverse[id(p["station"])],
              "baseline_between_minus_station":p["between"]["temp_total"]-p["station"]["temp_total"]}
             for p in pairs]
    return matches,qa,matchqa

def analyze(path,eivepath):
    plots,originalqa=parse_xlsx_matrix(Path(path))
    names=species_names_from_publisher(path)
    eive,eiveqa=eive_table(eivepath)
    geo,geoqa=site_metadata(path)
    cover,sourceqa=read_original_plant_values(path,plots)
    panels=[]
    for year,traitgate,matchgate in PANELS:
        pairs,origqa,matchqa=cohort(plots,names,eive,geo,year,traitgate,matchgate)
        p={"year":year,"trait_gate":traitgate,"match_gate":matchgate,
           "eligibility":origqa,"matched_source":matchqa,"weight_sensitivity":[]}
        if not pairs:
            p["status"]="no qualifying pairs";panels.append(p);continue
        baseline=mean(x["baseline_between_minus_station"] for x in pairs)
        for scheme in WEIGHTINGS:
            result=group_comparison(cover,pairs,names,eive,year,traitgate,scheme,"position")
            if scheme=="presence" and not math.isclose(
                    baseline,result["mean_between_minus_station_track_relative_interior"],
                    rel_tol=0,abs_tol=1e-10):
                raise AssertionError("Original v332 presence-weighted source baseline mismatch")
            p["weight_sensitivity"].append(result)
        panels.append(p)
    return {"version":"v334","source_original_Dryad_SHA256":DRYAD_SHA256,
      "EIVE_original_MD5":eiveqa["EIVE_original_md5"],
      "source_cover_numeric_integrity":sourceqa,
      "aggregate_coordinate_metadata":geoqa,
      "panels":panels,
      "warnings":[
        "Original source numeric cover values are not assumed to be biomass, percent cover or plant counts",
        "These are source abundance-weighted CWM compositional statistics, not measured thermal microclimate",
        "Source column of cover is species recorded value but source units must be confirmed from original protocol before biological interpretation",
        "Matching and trait coverage cohort fixed from v332; no outcome-driven rematching",
        "All site pairs from one Kashmir railway corridor, no independent climatic regions",
        "Original 2017 sparse EIVE trait cohort is not the prespecified strict primary and should not be interpreted as confirmatory",
        "Rail station versus between differs in settlement, human management and botanical source pool",
        "No train carriage, natural-population recruitment, microclimate recording or plant genetic adaptation observed",
        "Only source aggregate QC and weighted contrasts archived, no site identities or species-by-plot source records published"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=analyze(a.kashmir_original,a.eive_original)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":"v334",
        "source_cover_numeric_integrity":out["source_cover_numeric_integrity"],
        "panels":out["panels"]},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
