#!/usr/bin/env python3
"""v327 genuine independent Kashmir RAIL Track-vs-Interior EIVE thermal reanalysis.

31 original field sites x 2014/2017 repeated, 4 source zones and 123 source
species, source-authenticated by frozen Dryad SHA256. Asian source vegetation
traits are matched to European EIVE consensus by EXACT Latin binomials only.
Report missingness BEFORE interpreting and no claims of train carriage or
causal climate-offset. Original sample rows or coordinates are not exported.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean
from openpyxl import load_workbook

from analyze_Kashmir_same_sites_reference_year_v254 import (
    DRYAD_SHA256, YEARS, GROUPS, parse_xlsx_matrix,
)
from analyze_Slovenia_EIVE_temperature_width_v321 import (
    eive_table, normal_name
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial

MAIN_ZONES=("Track","Interior Plot 2")
ALTERNATE_ZONES=("Track","Interior Plot 1")
SENSITIVITIES=((3,5,.5),(1,3,.3))
METRICS=("position","width")
EXPECTED_SPECIES=123

def species_names_from_publisher(path):
    book=load_workbook(path,read_only=True,data_only=True)
    try:
        h=next(book.worksheets[0].values)
    finally:
        book.close()
    if len(h)!=129 or len(h[6:])!=EXPECTED_SPECIES:
        raise ValueError("Unexpected original Kashmiri source 129-column header")
    names=[]
    for value in h[6:]:
        if not isinstance(value,str) or not value.strip():
            raise ValueError("Invalid source biological species header")
        accepted=binomial(value)
        names.append(normal_name(accepted) if accepted else None)
    if len(set(v for v in names if v)) !=len([v for v in names if v]):
        raise ValueError("Duplicate source accepted species name after binomial parsing")
    return names

def cohort(plots):
    allsites=sorted({s for s,_,_ in plots})
    if len(allsites)!=31 or sum(not x.startswith("Between ") for x in allsites)!=16:
        raise ValueError("Original Kashmir spatial-site census drift")
    return allsites

def profile(indices,names,eive,metric,systems):
    if indices is None:
        return None
    # Source parser fails closed on four missing XLSX species cells.
    if not indices:
        return {"n_source":0,"n_eive":0,"coverage":0.0,"cwm":None}
    matched=[eive[names[i]][metric] for i in indices
             if names[i] in eive and eive[names[i]]["systems"]>=systems]
    return {"n_source":len(indices),"n_eive":len(matched),
            "coverage":len(matched)/len(indices),
            "cwm":mean(matched) if matched else None}

def site_effect(plots,names,eive,site,year,zones,metric,systems,min_n,min_coverage):
    a=profile(plots.get((site,year,zones[0])),names,eive,metric,systems)
    b=profile(plots.get((site,year,zones[1])),names,eive,metric,systems)
    if a is None or b is None:
        return None,{"reason":"missing_source_record"}
    if min(a["n_eive"],b["n_eive"])<min_n or min(a["coverage"],b["coverage"])<min_coverage:
        return None,{"reason":"trait_coverage_gate","track":a,"reference":b}
    if a["cwm"] is None or b["cwm"] is None:
        return None,{"reason":"no_valid_EIVE_numeric"}
    return a["cwm"]-b["cwm"],{"reason":"eligible","track":a,"reference":b}

def panel(plots,names,eive,year,zones,metric,systems,min_n,min_coverage):
    sites=cohort(plots)
    eligible=[]
    reasons=Counter()
    coverage=defaultdict(list)
    for site in sites:
        group="between" if site.startswith("Between ") else "station"
        delta,meta=site_effect(plots,names,eive,site,year,zones,metric,systems,min_n,min_coverage)
        reasons[meta["reason"]]+=1
        if meta.get("track"):
            coverage["track"].append(meta["track"]["coverage"])
            coverage["interior"].append(meta["reference"]["coverage"])
        if delta is not None:
            eligible.append((site,group,delta))
    stats={}
    for group in ("all","station","between"):
        samples=[x[2] for x in eligible if group=="all" or x[1]==group]
        stats[group]={"sites":len(samples),"track_higher_sites":sum(v>0 for v in samples),
                      "interior_higher_sites":sum(v<0 for v in samples),
                      "equal_sites":sum(v==0 for v in samples),
                      "mean_track_minus_interior":mean(samples) if samples else None,
                      "min":min(samples) if samples else None,
                      "max":max(samples) if samples else None}
    loo=None
    allvals=[x[2] for x in eligible]
    if len(allvals)>2:
        deleted=[mean(allvals[:i]+allvals[i+1:]) for i in range(len(allvals))]
        loo={"min":min(deleted),"max":max(deleted),
             "all_leave_one_source_site_out_positive":all(x>0 for x in deleted),
             "all_leave_one_source_site_out_negative":all(x<0 for x in deleted)}
    return {
        "source_year":year,"railway_zone":zones[0],"source_interior_reference":zones[1],
        "metric":metric,"minimum_EIVE_source_indicator_systems":systems,
        "minimum_numeric_taxa_per_zone":min_n,"minimum_coverage_fraction_per_zone":min_coverage,
        "source_sites_in_design":len(sites),
        "excluded_source_sites_by_reason":dict(reasons),
        "source_plot_coverage_across_nonmissing_site_pairs":{
           "track_mean":mean(coverage["track"]) if coverage["track"] else None,
           "reference_mean":mean(coverage["interior"]) if coverage["interior"] else None},
        "observed":stats,"leave_one_spatial_site_out":loo,
        "_effects":{site:v for site,group,v in eligible}
    }

def temporal(panel2014,panel2017):
    a=panel2014["_effects"];b=panel2017["_effects"]
    shared=sorted(set(a)&set(b))
    x=[b[s]-a[s] for s in shared]
    group_stats={}
    for group in ("all","station","between"):
        y=[b[s]-a[s] for s in shared
           if group=="all" or ("between" if s.startswith("Between ") else "station")==group]
        group_stats[group]={"n_same_original_sites":len(y),
            "mean_2017_minus_2014_of_track_interior_delta":mean(y) if y else None,
            "more_positive_sites":sum(v>0 for v in y),
            "more_negative_sites":sum(v<0 for v in y)}
    return {"n_common_spatial_sites":len(shared),
        "repeated_year_field_sites_NOT_independent_studies":True,
        "within_site_temporal":group_stats}

def analyze(source_xlsx,eive_xlsx):
    plots,source_qa=parse_xlsx_matrix(Path(source_xlsx))
    names=species_names_from_publisher(source_xlsx)
    eive,eive_qa=eive_table(eive_xlsx)
    source_named={s for s in names if s}
    source_match={s for s in source_named if s in eive}
    source_match3={s for s in source_match if eive[s]["systems"]>=3}
    out={"version":"v327","original_kashmir_source_sha256":DRYAD_SHA256,
         "original_dryad_version":99735,"original_kashmir":source_qa,
         "independent_EIVE_source_qa":eive_qa,
         "original_123_species_headers_with_simple_binomial":len(source_named),
         "original_taxa_matched_exact_EIVE_n>=1":len(source_match),
         "original_taxa_matched_exact_EIVE_n>=3":len(source_match3),
         "matched_source_species_are_European_EIVE_NOT_independent_measured_Asian_thermal_niches":True,
         "panels":[],"repeated_site_changes":[]}
    for zones in (MAIN_ZONES,ALTERNATE_ZONES):
        for systems,minn,coverage in SENSITIVITIES:
            for metric in METRICS:
                byyear={}
                for yr in YEARS:
                    result=panel(plots,names,eive,yr,zones,metric,systems,minn,coverage)
                    byyear[yr]=result
                    out["panels"].append({k:v for k,v in result.items() if k!="_effects"})
                t=temporal(byyear[2014],byyear[2017])
                out["repeated_site_changes"].append({"zones":zones,"metric":metric,
                    "minimum_source_systems":systems,
                    "min_n":minn,"min_coverage":coverage,**t})
    out["limitations"]=[
      "Same railway corridor in Kashmir with narrow original 1591–1741m gradient, not global climate gradient",
      "Interior Plot 1/2 at ca. 50–150m from line and same original 100m2 protocol are local railway-adjacent references, NOT verified pristine/unaffected forest",
      "Source 2014/2017 repeats are same field sites, not two new independent railway studies",
      "Exact EIVE indicator name match for subset of Kashmir flora; source EIVE calibrated for Europe and not validated Asian climatic physiological tolerances",
      "Source missing plot species cells are not zeros; any affected plot failed strict completeness",
      "No direct measured microclimate, seed-carriage, local disturbance controls or causal climate-range release",
      "Ecological EIVE centre and model niche-width indicator do not measure temperature in degrees C",
      "Any site sign/inference is exploratory under spatially dependent source locations"]
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=analyze(a.kashmir_original,a.eive_original)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({k:result[k] for k in ("version",
        "original_taxa_matched_exact_EIVE_n>=1",
        "original_taxa_matched_exact_EIVE_n>=3",
        "panels","repeated_site_changes")},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
