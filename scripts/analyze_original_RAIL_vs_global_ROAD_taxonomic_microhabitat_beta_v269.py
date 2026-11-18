#!/usr/bin/env python3
"""v26.9 source-authentic RAIL and global mountain ROAD matched-plot taxonomic beta.

Original sources:
- Kashmir railway 2014,2017: Rashid et al (2021) Dryad Matrix_plot.xlsx,
  SHA256=529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474.
  Source 31 sites, 4 original 100 m2 zones, common 27 complete sites.
- Global MIREN mountain ROADS: source 167804-row CSV v3, original MD5
  b640195cc44f252ccb6ddb28eb6fff13. Plot1 roadside vs Plot3 interior,
  both original 100 m2 vegetation quadrats. First site-year with both plots
  only, each stable original site used once. Main subset: earliest observed
  qualifying first survey >=2012, guarding 2007/08 alien-focused protocol.

DO NOT pool the studies. Distance to interior, substrate, vegetation, years
and protocol differ. The beta decomposition describes recorded identities,
not the movement of populations through time or train-seed dispersal.

For each original site compare same-site corridor and local reference:
  a=number shared, b=corridor-only, c=reference-only
  beta_jac=(b+c)/(a+b+c)
  beta_jtu=2min(b,c)/(a+2min(b,c))
  beta_jne=beta_jac-beta_jtu
  beta_turnover_fraction=beta_jtu/beta_jac when denominator >0
These are pairwise Baselga Jaccard replacement vs nestedness-resultant
algebraic components, NOT evidence of adaptation or actual taxon succession.

Road site mean -> equal road -> equal source Region code, min 4 sites
per Road & 2 Roads per Region. Main observed all-species lists, and
79 technical within-site richness-equalized repeated draws (not field CI).
The original Excel/CSV, taxon-by-plot rows, site/road IDs and GPS are NOT
written to output, only regional aggregate codes and summaries.
"""
from __future__ import annotations

import argparse
import json
import random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix,matching_sites
from audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 import source_read

METRICS=("corridor_richness","reference_richness","shared_taxa",
         "corridor_only_taxa","reference_only_taxa","beta_jac",
         "beta_jtu_replacement","beta_jne_nestedness","replacement_fraction",
         "corridor_unshared_fraction","corridor_strict_subset")
RARE_DRAWS=79
SEED=269118
MIN_SITES=4
MIN_ROADS=2
KASHMIR_ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
KASHMIR_REFERENCES=("Embankment","Interior Plot 1","Interior Plot 2")


def beta(a:set,b:set)->dict:
    if not a or not b:
        raise ValueError("Both original plots must be recorded and nonempty")
    shared=len(a&b)
    corridor_only=len(a-b)
    reference_only=len(b-a)
    union=shared+corridor_only+reference_only
    if union!=len(a|b) or union==0:
        raise AssertionError("Original source species-incidence union inconsistent")
    jac=(corridor_only+reference_only)/union
    small=min(corridor_only,reference_only)
    replacement=2*small/(shared+2*small) if shared+2*small else 0.
    nested=jac-replacement
    if not -1e-12<=nested<=1 or not 0<=replacement<=1 or not 0<=jac<=1:
        raise AssertionError("Original source beta partition outside unit interval")
    if abs(jac-(nested+replacement))>1e-12:
        raise AssertionError("Replacement and nestedness do not add to beta_jac")
    return {
        "corridor_richness":float(len(a)),
        "reference_richness":float(len(b)),
        "shared_taxa":float(shared),
        "corridor_only_taxa":float(corridor_only),
        "reference_only_taxa":float(reference_only),
        "beta_jac":jac,
        "beta_jtu_replacement":replacement,
        "beta_jne_nestedness":max(0.,nested),
        "replacement_fraction":(replacement/jac) if jac>0 else None,
        "corridor_unshared_fraction":corridor_only/len(a),
        "corridor_strict_subset":float(corridor_only==0),
    }


def average(records:list[dict])->dict:
    if not records:raise ValueError("Missing genuine field observations")
    out={}
    for k in METRICS:
        vals=[r[k] for r in records if r.get(k) is not None]
        out[k]=mean(vals) if vals else None
    return out


def rarefied_site(a:set,b:set,site_seed:int,draws:int=RARE_DRAWS)->dict:
    if draws<5:raise ValueError("Need at least five technical rarefaction draws")
    n=min(len(a),len(b))
    if n<1:raise ValueError("Empty original source plot")
    rng=random.Random(site_seed)
    aa=sorted(a);bb=sorted(b)
    vals=[beta(set(rng.sample(aa,n)),set(rng.sample(bb,n)))
          for _ in range(draws)]
    r=average(vals)
    if abs(r["corridor_richness"]-n)>1e-10 or abs(r["reference_richness"]-n)>1e-10:
        raise AssertionError("Source equal richness technical randomization failed")
    return r


def road_first_surveys(source:dict)->tuple[list[dict],dict]:
    plots=source["plants"];proto=source["protocol"]
    same=defaultdict(dict)
    for (region,road,site,year,zone),taxa in plots.items():
        if zone in ("1","3"):
            same[(region,road,site,year)][zone]=taxa
    complete=defaultdict(list)
    n_complete=0
    for (reg,road,site,year),zones in same.items():
        if not all(z in zones for z in ("1","3")):continue
        n_complete+=1
        if proto[(reg,road,site,year,"1")]!=proto[(reg,road,site,year,"3")]:
            raise ValueError("Original roadside and interior source differ in protocol in year")
        complete[(reg,road,site)].append(year)
    first=[]
    for (reg,road,site),years in sorted(complete.items()):
        year=min(years)
        a=set(plots[reg,road,site,year,"1"])
        b=set(plots[reg,road,site,year,"3"])
        if not a or not b:raise ValueError("Unexpected empty original source paired ROAD plant plot")
        first.append({"region":reg,"road":road,"site":site,"year":year,
                      "corridor":a,"reference":b})
    if n_complete!=1876 or len(first)!=1005:
        raise ValueError(f"Frozen MIREN v9.9 genuinely paired 1876 site-years/1005 unique original sites drift: {n_complete}/{len(first)}")
    yearcounts=Counter(x["year"] for x in first)
    return first,{
        "all_source_complete_strict_plot1_plot3_site_years":n_complete,
        "all_source_stable_sites_with_first_complete_year":len(first),
        "original_first_survey_year_counts":dict(sorted(yearcounts.items())),
        "first_source_survey_2012plus_sites":sum(x["year"]>=2012 for x in first),
        "first_source_survey_2017plus_sites":sum(x["year"]>=2017 for x in first),
        "no_source_repeats_as_independent_sites":True
    }


def select_miren_sites(first:list[dict],period:str)->list[dict]:
    if period=="first_anyyear": return first
    if period=="first_2012plus":return [r for r in first if r["year"]>=2012]
    if period=="first_2017plus":return [r for r in first if r["year"]>=2017]
    raise ValueError("Unspecified original MIREN first-year cohort")


def road_macro(rows:list[dict],year_condition:str,draws:int=RARE_DRAWS)->dict:
    rows=select_miren_sites(rows,year_condition)
    byroad=defaultdict(list)
    for item_index,x in enumerate(rows):
        raw=beta(x["corridor"],x["reference"])
        rare=rarefied_site(x["corridor"],x["reference"],SEED+101*item_index,draws)
        byroad[(x["region"],x["road"])].append({"raw":raw,"equalized":rare})
    valid_roads={k:records for k,records in byroad.items() if len(records)>=MIN_SITES}
    regions=defaultdict(list)
    for (region,road),records in valid_roads.items():
        regions[region].append({
            "n_source_sites":len(records),
            "raw":average([v["raw"] for v in records]),
            "equalized":average([v["equalized"] for v in records])
        })
    qualified={region:rs for region,rs in sorted(regions.items()) if len(rs)>=MIN_ROADS}
    if year_condition=="first_2012plus":
        if len(qualified)<8 or sum(len(x) for x in qualified.values())<15:
            raise ValueError("Original modern MIREN source macro region/road support missing")
    byregion={}
    for region,rs in qualified.items():
        byregion[region]={
            "n_source_roads":len(rs),
            "n_source_sites":sum(x["n_source_sites"] for x in rs),
            "raw":average([x["raw"] for x in rs]),
            "equalized":average([x["equalized"] for x in rs])
        }
    if not byregion:
        return {"identified":False,"source_regions":0,
                "original_candidate_qualified_sites":len(rows)}
    overall={}
    for mode in ("raw","equalized"):
        records=[v[mode] for v in byregion.values()]
        estimated=average(records)
        estimated["n_region_codes"] = len(records)
        estimated["replacement_fraction_of_aggregate_beta"]=(estimated["beta_jtu_replacement"]/estimated["beta_jac"] if estimated["beta_jac"] else None)
        estimated["regions_more_replacement_than_nestedness"] = sum(x["beta_jtu_replacement"]>x["beta_jne_nestedness"] for x in records)
        estimated["regions_corridor_taxa_only_gt_reference_taxa_only"] = sum(x["corridor_only_taxa"]>x["reference_only_taxa"] for x in records)
        estimated["regions_corridor_richness_gt_reference_richness"] = sum(x["corridor_richness"]>x["reference_richness"] for x in records)
        overall[mode]=estimated
    return {
        "identified":True,"cohort":year_condition,
        "original_candidate_first_site_count":len(rows),
        "original_roads_meeting_four_site_support_before_region_filter":len(valid_roads),
        "original_regions_meeting_two_road_support":len(byregion),
        "original_roads_in_supported_regions":sum(v["n_source_roads"] for v in byregion.values()),
        "original_sites_in_supported_regions":sum(v["n_source_sites"] for v in byregion.values()),
        "observed_mean_equal_road_equal_region":overall,
        "original_source_region_code_aggregates":byregion,
        "primary_units":"site -> equal original road -> equal original Region code",
        "original_source_region_codes_not_separate_randomly_selected_countries":True,
        "not_railway_or_train_effect":True
    }


def kashmir_panels(path:Path,draws:int=RARE_DRAWS)->dict:
    plots,source=parse_xlsx_matrix(path)
    results={}
    for group in ("station","between"):
        sites=matching_sites(plots,group,KASHMIR_ZONES)
        if len(sites)!=(15 if group=="station" else 12):
            raise ValueError("Original source complete 4-zone both-year Kashmir spatial cohort drift")
        results[group]={}
        for year in (2014,2017):
            results[group][str(year)]={}
            for z,ref in enumerate(KASHMIR_REFERENCES):
                source_vals=[];equalized=[]
                for i,site in enumerate(sites):
                    a=plots[site,year,"Track"];b=plots[site,year,ref]
                    if a is None or b is None:raise ValueError("Kashmir source site and year zone missing")
                    source_vals.append(beta(a,b))
                    equalized.append(rarefied_site(a,b,SEED+100003*i+1001*year+z,draws))
                results[group][str(year)][ref]={
                    "source_spatial_sites":len(sites),
                    "raw":average(source_vals),
                    "equalized":average(equalized),
                    "original_two_zone_100m2_site_pair":True
                }
    return {
        "genuine_source_QA":source,
        "same_original_complete_spatial_sites_station":15,
        "same_original_complete_spatial_sites_between":12,
        "results":results,
        "train_traffic_measured":False,
        "Kashmir_RAIL_never_pooled_with_MIREN_ROADS":True
    }


def run(miren:Path,kashmir:Path,draws:int=RARE_DRAWS)->dict:
    if draws<5:raise ValueError("Need source rarefaction draws >=5")
    roads,md5=source_read(miren)
    first,firstqa=road_first_surveys(roads)
    source= {
        p:road_macro(first,p,draws)
        for p in ("first_2012plus","first_2017plus","first_anyyear")
    }
    if not source["first_2012plus"]["identified"]:
        raise ValueError("Modern full flora MIREN 2012+ original source failed")
    rail=kashmir_panels(kashmir,draws)
    return {
        "version":"v26.9",
        "observed_real_original_sources":2,
        "original_MIREN_ROAD_Zenodo_DOI":"10.5281/zenodo.15749536",
        "original_Kashmir_RAIL_Dryad_DOI":"10.5061/dryad.0p2ngf20h",
        "source_MIREN_QA":md5,
        "MIREN_original_first_site_year_QA":firstqa,
        "MIREN_ROAD_source_matched_microhabitat_beta":source,
        "Kashmir_RAIL_same_source_microhabitat_beta":rail,
        "site_area_both_source_studies_m2":100,
        "species_binomials_NOT_directly_combined_across_databases":True,
        "baselga_pairwise_jaccard_beta_jac_equals_beta_jtu_plus_beta_jne":True,
        "technical_equal_richness_draws_per_original_site":draws,
        "technical_rarefaction_NOT_real_95pct_abstract_populations":True,
        "two_original_studies_not_multiple_independent_study_replicates":True,
        "ROAD_and_RAIL_interior_distance_source_protocol_vegetation_temporal_confounded":True,
        "no_field_causal_train_seed_transport_mowing_gust_heritable_adaptation":True,
        "no_original_plant_by_site_name_matrix_or_GPS_republished":True
    }


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--miren-csv",type=Path,required=True)
    p.add_argument("--kashmir-xlsx",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=RARE_DRAWS)
    a=p.parse_args()
    result=run(a.miren_csv,a.kashmir_xlsx,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for cohort,score in result["MIREN_ROAD_source_matched_microhabitat_beta"].items():
        print(json.dumps({"version":"v26.9","study":"MIREN_ROAD","cohort":cohort,
          "source_regions":score.get("original_regions_meeting_two_road_support"),
          "source_roads":score.get("original_roads_in_supported_regions"),
          "source_sites":score.get("original_sites_in_supported_regions"),
          "raw":score.get("observed_mean_equal_road_equal_region",{}).get("raw"),
          "equalized":score.get("observed_mean_equal_road_equal_region",{}).get("equalized")},ensure_ascii=False))
    for group in ("station","between"):
        for year in (2014,2017):
            raw=result["Kashmir_RAIL_same_source_microhabitat_beta"]["results"][group][str(year)]["Interior Plot 2"]
            print(json.dumps({"version":"v26.9","study":"Kashmir_RAIL",
              "group":group,"year":year,"reference":"Interior Plot 2",
              "sites":raw["source_spatial_sites"],"raw":raw["raw"],"equalized":raw["equalized"]},ensure_ascii=False))


if __name__=="__main__":main()
