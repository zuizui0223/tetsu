#!/usr/bin/env python3
"""v10.1: Baselga (2012) Jaccard species-turnover vs nestedness within
genuine globally standardized MIREN mountain ROAD vegetation plots.

Site × survey-year matched Alien species sets in original fixed 100m²
Plot1 (near road) vs Plot3 (52–102m interior). This is NOT an observed
RAILWAY vegetation dataset or a proof of vehicle seed carriage.

For shared a, road-only b, interior-only c:
  beta_jac=(b+c)/(a+b+c)
  beta_jtu=2min(b,c)/(a+2min(b,c))
  beta_jne=beta_jac-beta_jtu
Decompose ONLY when both plots have >=1 recorded group taxa; if
a plot has zero Alien spp, show a separate source support category
and never invent a biological turnover component or source absence.
The first observed same-site/year pair per unique road site is primary.
Region-level sign comparisons require >=10 / >=20 BOTH-positive pairs.
Native plants in 2017+ protocols are a non-road-specific compositional
comparison; not a randomized biological negative control.

Source: original MIREN Zenodo 15749536 v3, 167804 original occurrence
rows; original MD5 b640195cc44f252ccb6ddb28eb6fff13. Reconcile
independently validated v9.9 1876 same-year pair /1005 unique road
sites and v10.0 484 repeat road-site results before analysis.

No mowing schedules, wind/slipstream, actual freight seeds, rail
plants, plasticity or evolved adaptation are measured.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from statistics import median,mean
import pandas as pd
import numpy as np

from analyze_MIREN_standardized_road_edge_vegetation_v99 import (
    DOI,SOURCE_MD5,SOURCE_ROWS,SITE_YEAR,SITE,
    native_status_safe)
from analyze_MIREN_roadside_alien_resurvey_turnover_v100 import (
    per_site_turnover)

PLOTS=("1","3")
SOURCE_GROUPS=("Alien","Native")
REGION_BOTH_POSITIVE_MIN=(10,20)
REGION_CODE_NOT_INDEPENDENT_COUNTRIES=True


def jaccard_components(edge,interior):
    """Both sides must have a recorded taxon to define a proper
    replacement/nestedness allocation, otherwise report invalid scope.
    """
    a=len(edge&interior)
    b=len(edge-interior)
    c=len(interior-edge)
    out={"shared_species":a,"roadside_only_species":b,
         "interior_only_species":c,
         "roadside_recorded_richness":len(edge),
         "interior_recorded_richness":len(interior),
         "roadside_minus_interior_richness":len(edge)-len(interior),
         "has_both_recorded_species_pools":bool(edge) and bool(interior)}
    if not edge or not interior:
        out.update({"jaccard_total":None,"jaccard_turnover":None,
                    "jaccard_nestedness_resultant":None,
                    "turnover_higher_than_nestedness":None})
        return out
    denom=a+b+c
    replace=a+2*min(b,c)
    if denom<=0 or replace<=0:
        raise ValueError("Both-positive species pools but undefined Baselga decomposition")
    total=(b+c)/denom
    turnover=2*min(b,c)/replace
    nestedness=total-turnover
    if not all(-1e-10<=v<=1+1e-10 for v in (total,turnover,nestedness)):
        raise ValueError("Jaccard beta decomposition nonphysical components")
    if not math.isclose(total,turnover+nestedness,abs_tol=1e-10):
        raise ValueError("Jaccard total is not conserved between turnover and nestedness")
    out.update({
       "jaccard_total":total,
       "jaccard_turnover":turnover,
       "jaccard_nestedness_resultant":max(0.0,nestedness),
       "turnover_higher_than_nestedness":turnover>nestedness})
    return out


def original_spatial_incidence(original):
    """Build actual field-collected SITE_YEAR Plot1 &3 taxon incidence.
    A recorded Native plant in one plot means Alien count zero can be
    observed under same plot; missing entire plot never means zero.
    """
    clean,qc=native_status_safe(original)
    plot_taxa=clean.groupby(list(SITE_YEAR)+["Plot","Status"])[
        "Accepted.Name.MIREN"].agg(lambda v:set(v)).to_dict()
    observed=clean.groupby(list(SITE_YEAR))["Plot"].agg(set)
    paired=sorted([key for key,v in observed.items() if set(PLOTS).issubset(v)])
    rows=[]
    for siteyear in paired:
        row=dict(zip(SITE_YEAR,siteyear))
        for group in SOURCE_GROUPS:
            edge=plot_taxa.get((*siteyear,"1",group),set())
            inner=plot_taxa.get((*siteyear,"3",group),set())
            result=jaccard_components(edge,inner)
            row.update({group+"_"+k:v for k,v in result.items()})
        rows.append(row)
    table=pd.DataFrame(rows).sort_values(list(SITE_YEAR)).reset_index(drop=True)
    if table.empty:
        raise ValueError("No real same-road/site/year MIREN Plot 1 & Plot 3 source pair")
    first=table.sort_values("Year").drop_duplicates(list(SITE),keep="first")
    if first.duplicated(list(SITE)).any():
        raise ValueError("Original mountain road sites are not de-duplicated")
    return clean,qc,table,first


def group_summary(table,group):
    if group not in SOURCE_GROUPS:
        raise ValueError("Only independently source-classified Alien or Native taxon group")
    pref=group+"_"
    a=table[pref+"roadside_recorded_richness"]
    b=table[pref+"interior_recorded_richness"]
    cohort=table[table[pref+"has_both_recorded_species_pools"]].copy()
    t=cohort[pref+"jaccard_turnover"].astype(float)
    n=cohort[pref+"jaccard_nestedness_resultant"].astype(float)
    total=cohort[pref+"jaccard_total"].astype(float)
    if not np.allclose(total,t+n,atol=1e-10):
        raise ValueError("Source total beta dissimilarity sum mismatch")
    decomposition={
        "n_same_original_site_year_plot1_vs_plot3_comparisons":len(table),
        "n_both_original_plot_group_taxa_nonzero":len(cohort),
        "n_roadside_positive_interior_zero_group_taxa":int(((a>0)&(b==0)).sum()),
        "n_roadside_zero_interior_positive_group_taxa":int(((a==0)&(b>0)).sum()),
        "n_both_plots_group_taxa_zero_but_plots_observed":int(((a==0)&(b==0)).sum()),
        "n_source_siteyears_without_defined_Jaccard_component":
            len(table)-len(cohort),
        "mean_roadside_only_recorded_group_taxa_per_pair":round(float(
            table[pref+"roadside_only_species"].mean()),6),
        "mean_interior_only_recorded_group_taxa_per_pair":round(float(
            table[pref+"interior_only_species"].mean()),6),
        "median_roadside_alien_richness_if_group_Alien_else_native":
            round(float(a.median()),6),
        "median_interior_alien_richness_if_group_Alien_else_native":
            round(float(b.median()),6),
        "group_both_positive_mean_Jaccard_total":round(float(total.mean()),6) if len(total) else None,
        "group_both_positive_mean_Jaccard_turnover":round(float(t.mean()),6) if len(t) else None,
        "group_both_positive_mean_Jaccard_nestedness":round(float(n.mean()),6) if len(n) else None,
        "group_both_positive_turnover_share_of_TOTAL_dissimilarity_sum":
            round(float(t.sum()/total.sum()),6) if total.sum()>0 else None,
        "n_both_positive_pairs_turnover_greater_nestedness":int((t>n).sum()),
        "n_both_positive_pairs_nestedness_greater_turnover":int((t<n).sum()),
        "no_zero_group_taxon_plot_imputed_as_biological_absence":True,
        "no_p_values_from_nonindependent_siteyear_comparisons":True
    }
    return decomposition


def region_summary(first,group,min_support):
    # Only compare genuine both-positive group taxon pools with at least
    # min_support separate ORIGINAL road-site sources per region.
    pref=group+"_"
    eligible=first[first[pref+"has_both_recorded_species_pools"]]
    regions=[]
    for code,x in eligible.groupby("Region"):
        if len(x)<min_support:continue
        turnover=x[pref+"jaccard_turnover"].astype(float)
        nested=x[pref+"jaccard_nestedness_resultant"].astype(float)
        total=x[pref+"jaccard_total"].astype(float)
        regions.append({
            "Region_code":code,
            "n_distinct_original_road_sites_both_taxon_pools_positive":len(x),
            "mean_Jaccard_turnover":round(float(turnover.mean()),6),
            "mean_Jaccard_nestedness":round(float(nested.mean()),6),
            "mean_Jaccard_total":round(float(total.mean()),6),
            "mean_turnover_minus_nestedness":round(float((turnover-nested).mean()),6),
            "proportion_total_dissimilarity_attributed_to_turnover":
                round(float(turnover.sum()/total.sum()),6) if total.sum()>0 else None,
            "n_sites_turnover_greater_nestedness":int((turnover>nested).sum()),
            "no_roadside_mowing_or_gust_measurement":True
        })
    vals=[r["mean_turnover_minus_nestedness"] for r in regions]
    return {
        "n_region_codes_with_min_support":len(regions),
        "min_both_positive_unique_first_site_pairs_per_region":min_support,
        "n_region_codes_turnover_dominant":sum(x>0 for x in vals),
        "n_region_codes_nestedness_dominant":sum(x<0 for x in vals),
        "n_region_codes_equal":sum(x==0 for x in vals),
        "median_region_equal_mean_turnover_minus_nestedness":round(
            float(median(vals)),6) if vals else None,
        "n_total_original_road_sites_from_these_region_codes":sum(
            r["n_distinct_original_road_sites_both_taxon_pools_positive"]
            for r in regions),
        "region_code_details":regions,
        "not_countries_or_randomized_experimental_treatments":True
    }


def run(original,prior99,prior100):
    if (prior99.get("complete") is not True or
        prior100.get("complete") is not True or
        prior99.get("geography",{}).get("n_matched_original_plot1_plot3_same_site_year_pairs")!=1876 or
        prior99.get("geography",{}).get("n_distinct_site_IDs_matching_both_plots")!=1005 or
        prior100.get("results",{}).get("n_original_matched_first_to_last_road_sites")!=484):
        raise ValueError("Original independently audited MIREN 1876 site-years, 1005 sites or 484 revisits drift")
    clean,qc,all_pairs,first=original_spatial_incidence(original)
    if len(all_pairs)!=1876 or len(first)!=1005:
        raise ValueError("Original MIREN same road site-year geography changed")
    regional={}
    for group in SOURCE_GROUPS:
        regional[group]={}
        for m in REGION_BOTH_POSITIVE_MIN:
            regional[group][f"min_{m}_both_positive_unique_first_sites"]=region_summary(
                first,group,m)
    # Protocol comparison for Native restricted to recent >=2017 survey
    # when both original Plot 1/Plot 3 used same year.
    after2017=first[first.Year>=2017]
    post2017=all_pairs[all_pairs.Year>=2017]
    out={
        "complete":True,
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_original_Zenodo_DOI":DOI,
        "actual_independent_field_survey":"global mountain ROAD not RAILWAY roadside vegetation",
        "source_original_MD5_reconciled_to_v9_8":SOURCE_MD5,
        "source_quality":qc,
        "original_matched_plot1_plot3_siteyears":len(all_pairs),
        "original_independent_unique_ROADsites_first_survey":len(first),
        "source_Alien_species_spatial_beta":{},
        "source_Native_species_spatial_beta_comparator":{},
        "region_first_original_roadsite_turnover_vs_nestedness":regional,
        "2017_or_later_protocol_compatibility":{},
        "NO_railway_seed_delivery_measured_mowing_wind_or_adaptation":True
    }
    for key,label in ((first,"first_unique_site"),(all_pairs,"all_repeated_siteyears")):
        out["source_Alien_species_spatial_beta"][label]=group_summary(key,"Alien")
        out["source_Native_species_spatial_beta_comparator"][label]=group_summary(key,"Native")
    for key,label in ((after2017,"first_unique_site_surveyed_2017plus"),(post2017,"all_siteyears_2017plus")):
        out["2017_or_later_protocol_compatibility"][label]={
          g:group_summary(key,g) for g in SOURCE_GROUPS}
    out["2017_or_later_protocol_compatibility"][
        "Alien_first_2017plus_min10_region_turnover_signs"]=region_summary(
           after2017,"Alien",10)
    out["2017_or_later_protocol_compatibility"][
        "Native_first_2017plus_min10_region_turnover_signs"]=region_summary(
           after2017,"Native",10)
    out["scientific_limitations"]=[
       "Data are original MIREN mountain-ROAD standardized 100m² vegetation plot species lists, NOT the ecological flora of an active RAIL line.",
       "MIREN source alien status is region-coded, not globally fixed for a taxon; two region-taxon nativity conflicts are conservatively unknown.",
       "Jaccard beta components only estimated when BOTH plots have >=1 source-recorded taxon in a given group; 0-group plots must be kept as separate environmental source states.",
       "Jaccard replacement vs nestedness are descriptive patterns of observed plant community composition; species turnover is NOT direct dispersal flux, mowing selection or evolved adaptation.",
       "Many independent recorded region codes are not statistically independent countries; 1005 nested road sites are not 1005 globally randomized experimental units.",
       "Native early source survey protocols (<2012) can under-record native plants, so 2017plus Native subset is the safer comparator; it still differs from an experimental control.",
       "No active railway freight activity, known viability of propagules, road-vehicle seeds, grass mowing schedule, gust velocity, heat exposure, seed mass, regrowth capacity, plasticity or genetics measured.",
       "Do not claim this reanalysis of MIREN consortium survey data as original field collection; cite data DOI and published methods."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for name in ("original-csv","verified99","verified100","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    raw=Path(a.original_csv)
    if hashlib.md5(raw.read_bytes()).hexdigest()!=SOURCE_MD5:
        raise ValueError("Original Zenodo MIREN v3 source MD5 incorrect")
    flora=pd.read_csv(raw,low_memory=False)
    if len(flora)!=SOURCE_ROWS:
        raise ValueError("MIREN original source data 167804-row identity mismatch")
    load=lambda x:json.loads(Path(getattr(a,x)).read_text(encoding="utf-8"))
    out=run(flora,load("verified99"),load("verified100"))
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "original_matched_siteyear_source":out["original_matched_plot1_plot3_siteyears"],
        "unique_first_site_source":out["original_independent_unique_ROADsites_first_survey"],
        "Alien":out["source_Alien_species_spatial_beta"]["first_unique_site"],
        "region_min10":out["region_first_original_roadsite_turnover_vs_nestedness"][
            "Alien"]["min_10_both_positive_unique_first_sites"],
        "region_min20":out["region_first_original_roadsite_turnover_vs_nestedness"][
            "Alien"]["min_20_both_positive_unique_first_sites"],
        "native_first_2017plus":out["2017_or_later_protocol_compatibility"][
            "first_unique_site_surveyed_2017plus"]["Native"],
        "no_railway_mowing_gust_adaptive_evolution_inferred":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
