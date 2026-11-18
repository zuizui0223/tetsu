#!/usr/bin/env python3
"""v28.8: Exact missing-life-span partial-identification across all Kashmir sites.

Instead of losing 24/27 sites at the >=70%-known-trait gate, use original
presence sets for every common source site. Unknown taxa have a SINGLE global
binary annual/nonannual classification consistent across all observations.
Analytical coefficient bounds are exact for the source incidence matrix, not
approximate bootstrap confidence intervals.

This DOES NOT estimate causal effects or species genetics, and does not make
the two field-survey years independent studies.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_Kashmir_annual_microsites_v285 import life_span_crosswalk, original_taxon_headers
from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix

YEARS=(2014,2017)
ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
COMPARISONS=(
    ("Track","Interior Plot 1"),
    ("Track","Interior Plot 2"),
    ("Track","Embankment"),
    ("Embankment","Interior Plot 1"),
    ("Embankment","Interior Plot 2"),
)
WEIGHTS={
    "2014": {2014:1.0},
    "2017": {2017:1.0},
    "two_year_mean": {2014:0.5,2017:0.5},
    "2017_minus_2014": {2014:-1.0,2017:1.0},
}


def common_source_sites(plots: dict)->list[str]:
    return sorted(s for s in {k[0] for k in plots}
                  if all(plots.get((s,y,z)) is not None for y in YEARS for z in ZONES))


def group_of(site: str)->str:
    return "between" if site.startswith("Between ") else "station"


def zone_fractions(plot: set[int], classified: dict[int,str],
                   total_taxa: int)->tuple[float,dict[int,float]]:
    if not plot:raise ValueError("Zero-richness original plot; no annual fraction defined")
    if any(not (0<=idx<total_taxa) for idx in plot):
        raise ValueError("Source plot contains taxon index outside genuine headers")
    k=len(plot)
    baseline=sum(classified.get(i)=="annual" for i in plot)/k
    unknown={i:1/k for i in plot if i not in classified}
    return baseline,unknown


def accumulate(dest:dict[int,float], source:dict[int,float], scale:float)->None:
    for idx,coef in source.items():
        dest[idx]=dest.get(idx,0.0)+scale*coef


def one_site_linear_contrast(plots:dict, classified:dict[int,str],
                             total_taxa:int, site:str,
                             zone:str, reference:str,
                             year_weights:dict[int,float])->tuple[float,dict[int,float]]:
    base=0.0
    terms={}
    for year,weight in year_weights.items():
        a0,a=zone_fractions(plots[site,year,zone],classified,total_taxa)
        b0,b=zone_fractions(plots[site,year,reference],classified,total_taxa)
        base+=weight*(a0-b0)
        accumulate(terms,a,weight)
        accumulate(terms,b,-weight)
    return base,terms


def partial_identification(site_values:list[tuple[float,dict[int,float]]])->dict:
    """Exact min/max of an affine function on species-consistent [0,1] labels.

    Constraints: each source-unknown species has one LIFE-HISTORY value across
    all sites, zones and years (not independent per observation). The
    optimization is separable: for each taxon j choose x_j=0 if c_j>=0,
    else x_j=1 for the lower bound; reverse for upper.
    """
    if not site_values:
        return {"source_sites":0,"status":"NO_SOURCE_SUPPORT",
                "sharp_lower":None,"sharp_upper":None,"sign_identified":None}
    n=len(site_values)
    fixed=sum(x[0] for x in site_values)/n
    global_coeff=collections.defaultdict(float)
    for _,this in site_values:
        accumulate(global_coeff,this,1.0/n)
    nonzero={i:c for i,c in global_coeff.items() if abs(c)>1e-12}
    lower=fixed+sum(min(0,c) for c in global_coeff.values())
    upper=fixed+sum(max(0,c) for c in global_coeff.values())
    all_annual=fixed+sum(global_coeff.values())
    if lower>1e-12:sign="POSITIVE_FOR_ALL_UNKNOWN_TAXON_LABELS"
    elif upper< -1e-12:sign="NEGATIVE_FOR_ALL_UNKNOWN_TAXON_LABELS"
    else:sign="UNIDENTIFIED_SIGN"
    return {"source_sites":n,
            "status":"PARTIALLY_IDENTIFIED",
            "sharp_lower":lower,"sharp_upper":upper,
            "all_unlabelled_nonannual_scenario":fixed,
            "all_unlabelled_annual_scenario":all_annual,
            "sign_identified":sign,
            "unlabelled_taxa_with_nonzero_global_coefficient":len(nonzero),
            "unknown_species_coefficient_l1":sum(abs(c) for c in nonzero.values()),
            "bound_model":"one binary annual status per original unknown taxon across all years/sites/zones",
            "not_a_field_confidence_interval":True}


def analyze(plots:dict,classified:dict[int,str],
            total_taxa:int=123)->dict:
    if any(not (0<=k<total_taxa) or v not in ("annual","biennial","perennial")
           for k,v in classified.items()):
        raise ValueError("Invalid life-span crosswalk or original source species index")
    source_sites=common_source_sites(plots)
    if not source_sites: raise ValueError("No common source sites across all 4 zones and both years")
    result={"original_taxon_columns":total_taxa,
            "original_taxa_classified":len(classified),
            "original_taxa_unclassified":total_taxa-len(classified),
            "strict_annual_classified_count":sum(x=="annual" for x in classified.values()),
            "known_perennial_classified_count":sum(x=="perennial" for x in classified.values()),
            "known_biennial_classified_count":sum(x=="biennial" for x in classified.values()),
            "same_all_zone_two_year_original_sites":len(source_sites),
            "station_source_sites":sum(group_of(x)=="station" for x in source_sites),
            "between_source_sites":sum(group_of(x)=="between" for x in source_sites),
            "source_years_in_one_original_study":list(YEARS),
            "study_replicates":1,
            "interpretation":"strict annual fraction of ALL source-recorded species; unknown original taxa globally consistent",
            "comparisons":{}}
    for zone,ref in COMPARISONS:
        label=zone+"__minus__"+ref
        result["comparisons"][label]={}
        for group in ("all","station","between"):
            sites=[s for s in source_sites if group=="all" or group_of(s)==group]
            result["comparisons"][label][group]={}
            for target,year_w in WEIGHTS.items():
                values=[one_site_linear_contrast(plots,classified,total_taxa,site,zone,ref,year_w)
                        for site in sites]
                result["comparisons"][label][group][target]=partial_identification(values)
    return result



def unknown_taxon_priority(plots:dict,classified:dict[int,str],
                           original_names:list[str],zone:str="Track",
                           ref:str="Interior Plot 1",
                           year_weights:dict[int,float] | None=None)->dict:
    """Rank original unresolved species by exact contribution to uncertainty width.

    This is a curation queue, not evidence that species are annual/perennial.
    It reports each species name, never its original plot-by-site occurrences.
    """
    weights=year_weights or WEIGHTS["two_year_mean"]
    sites=common_source_sites(plots)
    values=[one_site_linear_contrast(plots,classified,len(original_names),
                                     site,zone,ref,weights) for site in sites]
    aggregate=collections.defaultdict(float)
    for _,coeff in values:
        accumulate(aggregate,coeff,1/len(values))
    missing=set(range(len(original_names)))-set(classified)
    ranked=sorted(
        ({"original_taxon":original_names[i],
          "source_header_index":i,
          "coefficient_if_annual":aggregate.get(i,0.0),
          "uncertainty_width_pp":100*abs(aggregate.get(i,0.0))}
         for i in missing),
        key=lambda row:(-row["uncertainty_width_pp"],row["original_taxon"]))
    full_width=sum(row["uncertainty_width_pp"] for row in ranked)
    for rank,item in enumerate(ranked,1):
        item["priority_rank"]=rank
    return {"unknown_original_taxa":len(missing),"full_sharp_range_width_pp":full_width,
            "top5_width_contribution_pp":sum(r["uncertainty_width_pp"] for r in ranked[:5]),
            "top5_fraction_of_width":(sum(r["uncertainty_width_pp"] for r in ranked[:5])/full_width
                                      if full_width else None),
            "ranked_unmatched_species":ranked,
            "method":"pre-existing source original name; life-span annotation still unverified"}

def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rail",type=Path,required=True,help="SHA-authenticated 2021 Dryad Matrix_plot.xlsx")
    ap.add_argument("--traits",type=Path,required=True,help="independently authenticated global trait crosswalk")
    ap.add_argument("--output",type=Path,required=True,help="Aggregate summary only; source species rows excluded")
    args=ap.parse_args()
    plots,raw_meta=parse_xlsx_matrix(args.rail)
    names=original_taxon_headers(args.rail)
    classified,trait_meta=life_span_crosswalk(args.traits,names)
    out={"version":"v28.8","source_original_rail_sha256":raw_meta["original_file_sha256"],
         "trait_crosswalk_sha256":trait_meta["crosswalk_sha256"],
         "partial_identification_not_observed_trait_complete_effect":True,
         "author_species_by_plot_records_not_redistributed":True,
         "original_study_effect":analyze(plots,classified,len(names)),
         "unclassified_original_taxa_curation_priority":unknown_taxon_priority(plots,classified,names)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    main=out["original_study_effect"]["comparisons"]["Track__minus__Interior Plot 1"]
    print(json.dumps({"version":out["version"],
       "source_original_sites":out["original_study_effect"]["same_all_zone_two_year_original_sites"],
       "source_classified":out["original_study_effect"]["original_taxa_classified"],
       "strict_annual_taxa":out["original_study_effect"]["strict_annual_classified_count"],
       "Track_minus_Interior1":main,
       "Track_minus_Interior2":out["original_study_effect"]["comparisons"]["Track__minus__Interior Plot 2"]["all"],
       "Track_minus_Embankment":out["original_study_effect"]["comparisons"]["Track__minus__Embankment"]["all"],
       "unknown_trait_priority":out["unclassified_original_taxa_curation_priority"],
       "output_path":str(args.output)},ensure_ascii=False))


if __name__=="__main__":
    main()
