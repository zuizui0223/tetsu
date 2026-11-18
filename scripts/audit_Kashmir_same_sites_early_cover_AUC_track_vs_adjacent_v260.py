#!/usr/bin/env python3
"""v26.0: within SAME original Kashmir spatial site, how much more strongly
does 2014 numeric abundance rank predict 2017 re-recording in railway Track
than Embankment or Interior 100m2 plots?

Paired comparisons only on original spatial sites with BOTH later-retained
and later-unrecorded 2014 taxa in EACH compared plot, fixed 2014/2017 site
cohort; numeric source cover is ordinal, no claimed percent/biomass unit.

The within-plot AUC in each environment is not a causal measure of track
stress or natural selection; lower-cover species may be harder to detect,
plots differ in richness/physical conditions, and plants are not marked
individuals. No source taxon/cover/site IDs are exported.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from statistics import mean,median

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix,matching_sites
from audit_Kashmir_2014_species_cover_persistence_2017_v259 import (
    early_numeric_cover, site_early_persistence
)

GROUPS=("station","between")
ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
COMPARATORS=("Embankment","Interior Plot 1","Interior Plot 2")
SEED=260118


def qtile(v:list[float],p:float) -> float:
    a=sorted(v);t=(len(a)-1)*p;i=int(t)
    return a[i] if i==len(a)-1 else a[i]*(1-t+i)+a[i+1]*(t-i)


def paired_summary(records:dict[str,dict[str,dict]],draws:int=1999,seed:int=SEED)->dict:
    if draws<99:raise ValueError("At least 99 field-unit bootstrap draws required")
    sites=sorted(records)
    if len(sites)<10:raise ValueError("Insufficient independent source spatial sites")
    outcomes={}
    for reference in COMPARATORS:
        eligible=[s for s in sites
                  if records[s]["Track"]["eligible_auc"] and records[s][reference]["eligible_auc"]]
        base={"source_reference_zone":reference,
              "original_common_site_cohort":len(sites),
              "paired_original_site_count":len(eligible),
              "minimum_paired_sites_for_estimation":5}
        if len(eligible)<5:
            base.update({"identified":False,"insufficient_identical_site_support":True})
            outcomes[reference]=base
            continue
        deltas={
            s:records[s]["Track"]["source_2014_cover_AUC"] -
                records[s][reference]["source_2014_cover_AUC"] for s in eligible
        }
        diffs=list(deltas.values())
        rng=random.Random(seed+sum(map(ord,reference)))
        boots=[mean(rng.choices(diffs,k=len(diffs))) for _ in range(draws)]
        leave=[mean(d for key,d in deltas.items() if key!=excluded) for excluded in eligible]
        base.update({
            "identified":True,
            "source_2014_cover_AUC_track_on_same_original_sites":
                mean(records[s]["Track"]["source_2014_cover_AUC"] for s in eligible),
            "source_2014_cover_AUC_reference_on_same_original_sites":
                mean(records[s][reference]["source_2014_cover_AUC"] for s in eligible),
            "source_2014_cover_AUC_track_minus_reference_paired_mean":mean(diffs),
            "source_2014_cover_AUC_track_minus_reference_paired_median":median(diffs),
            "n_sites_track_auc_gt_reference_auc":sum(d>0 for d in diffs),
            "n_sites_track_auc_eq_reference_auc":sum(d==0 for d in diffs),
            "n_sites_track_auc_lt_reference_auc":sum(d<0 for d in diffs),
            "site_cluster_resampling_technical": {
                "draws":draws,"q025_not_randomized_field_CI":qtile(boots,.025),
                "q975_not_randomized_field_CI":qtile(boots,.975),
                "fraction_positive":sum(x>0 for x in boots)/len(boots)
            },
            "leave_one_original_site_out":{
                "min":min(leave),"max":max(leave)},
            "per_plot_2014_source_taxa_start_count_in_matched_sites":{
                "track":sum(records[s]["Track"]["n_original_2014_taxa"] for s in eligible),
                "reference":sum(records[s][reference]["n_original_2014_taxa"] for s in eligible)
            },
            "per_plot_2014_retained_site_taxa_count_in_matched_sites":{
                "track":sum(records[s]["Track"]["n_persistent_recorded_2017"] for s in eligible),
                "reference":sum(records[s][reference]["n_persistent_recorded_2017"] for s in eligible)
            },
            "no_randomized_zone_assignment_or_shared_plant_risk_set":True
        })
        outcomes[reference]=base
    return outcomes


def run(path:Path,draws:int=1999)->dict:
    presence,source=parse_xlsx_matrix(path)
    numeric,cover_qc=early_numeric_cover(path,presence)
    groups={}
    for i,group in enumerate(GROUPS):
        sites=matching_sites(presence,group,ZONES)
        if len(sites)!=(15 if group=="station" else 12):
            raise ValueError("Original 27 repeated same-zone sites unsupported")
        rec={}
        for site in sites:
            rec[site]={}
            for zone in ZONES:
                early=numeric[(site,2014,zone)]
                late=presence[(site,2017,zone)]
                if early is None or late is None:
                    raise ValueError("Not genuinely complete source zone and year")
                rec[site][zone]=site_early_persistence(early,late)
        groups[group]={
            "original_independent_spatial_sites":len(sites),
            "strict_same_site_2014_2017_all_four_original_100m2_zones":True,
            "matched_zone_contrasts":paired_summary(rec,draws,SEED+i*611)
        }
    return {
        "version":"v26.0",
        "source_original_Dryad_matrix_SHA256":source["original_file_sha256"],
        "source_site_support":source,
        "source_cover_read_QA":cover_qc,
        "results":groups,
        "metric":"Within-plot original 2014 cover rank AUC predicting 2017 original species re-recording",
        "no_railway_traffic_or_mowing_or_wind_causes_tested":True,
        "not_a_proven_survival_or_recruitment_measure":True,
        "pairwise_comparison_keeps_same_original_site_not_identical_taxa":True,
        "bootstrap_intervals_are_source_site_resampling_descriptive_not_unbiased_field_CI":True,
        "no_raw_cover_or_species_site_records_output":True
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--bootstrap",type=int,default=1999)
    args=parser.parse_args()
    out=run(args.source,args.bootstrap)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    for group,section in out["results"].items():
        for ref,item in section["matched_zone_contrasts"].items():
            print(json.dumps({"version":"v26.0","group":group,"reference":ref,
                "source_paired_sites":item["paired_original_site_count"],
                "track_auc":item.get("source_2014_cover_AUC_track_on_same_original_sites"),
                "reference_auc":item.get("source_2014_cover_AUC_reference_on_same_original_sites"),
                "paired_auc_difference":item.get("source_2014_cover_AUC_track_minus_reference_paired_mean"),
                "technical_resampling":item.get("site_cluster_resampling_technical")
            },ensure_ascii=False))


if __name__=="__main__":
    main()
