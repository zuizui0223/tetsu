#!/usr/bin/env python3
"""v10.0: same-plot repeated MIREN mountain-road Alien taxa lists.

Source is the ORIGINAL MD5-fixed 2007–2023 MIREN Zenodo v3
167804-record independent field survey, not GBIF general presence.
Original Plot1 roadside and Plot3 52–102m interior, same Region+
Road+Site+Year, no replacement plot IDs. For each original road-site
with >=2 observed pair years, compare earliest/latest recorded
Alien species sets at each plot:

recorded GAIN = LATE speciesKey name NOT in EARLY positive plot list
recorded LOSS = EARLY taxon name NOT in LATE positive plot list
relative GAIN = GAIN / UNION(early,late) of each plot (undefined if
UNION=0, NEVER impute zero).

This is a FIELD-RECORDED apparent temporal gain, NOT new naturalization
or actual viable seeds from vehicles; first/last plot surveys can
miss present plant species. Within-region summary of 11 regions
prevents claiming hundreds of independent regional replicates.
No real measured road mowing/wind, and this is not railway field data.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from statistics import median,mean
import numpy as np
import pandas as pd

from analyze_MIREN_standardized_road_edge_vegetation_v99 import (
    SITE_YEAR,SITE,SOURCE_MD5,SOURCE_ROWS,DOI,
    native_status_safe,plot_species_richness,matched_siteyears)

STRICT_PLOTS=("1","3")
MAX_PLOT_REGIONS_WITH_REPEATED_SOURCE=100
MIN_REGION_REPEAT_SITES=10


def per_site_turnover(plant):
    """ALL original Plantae source plots can be selected as a sampled
    site/year, with Alien count ZERO when Native records exist; missing
    whole plot/year is not taken as true absence or zero.
    """
    src=plant[plant["Plot"].isin(STRICT_PLOTS)].copy()
    lists={}
    for key,subset in src.groupby(list(SITE_YEAR)+["Plot"]):
        if key in lists:raise ValueError("Original plot/year taxon set duplicated after aggregate")
        lists[key]=set(subset.loc[subset["Status"]=="Alien","Accepted.Name.MIREN"])
    bysite=defaultdict(lambda:defaultdict(set))
    for (reg,road,site,year,plot) in lists:
        bysite[(reg,road,site)][year].add(plot)
    results=[]
    for (region,road,site),times in sorted(bysite.items()):
        complete=sorted(y for y,plots in times.items() if set(STRICT_PLOTS)<=plots)
        if len(complete)<2:continue
        early,late=complete[0],complete[-1]
        item={"Region":region,"Road":road,"Site":site,
              "initial_year":early,"latest_year":late,
              "years_between_first_and_last":late-early}
        for plot in STRICT_PLOTS:
            old=lists[(region,road,site,early,plot)]
            new=lists[(region,road,site,late,plot)]
            gained=len(new-old)
            lost=len(old-new)
            retained=len(new&old)
            denom=len(new|old)
            item.update({
                "gain"+plot:gained,
                "loss"+plot:lost,
                "retained"+plot:retained,
                "union"+plot:denom,
                "early"+plot:len(old),
                "late"+plot:len(new),
                "gain_fraction_union"+plot:gained/denom if denom>0 else np.nan,
                "loss_fraction_union"+plot:lost/denom if denom>0 else np.nan,
            })
        item["gain_roadside_minus_interior"]=item["gain1"]-item["gain3"]
        item["loss_roadside_minus_interior"]=item["loss1"]-item["loss3"]
        item["difference_in_gain_fraction_roadside_minus_interior"]=(
           item["gain_fraction_union1"]-item["gain_fraction_union3"])
        item["difference_in_earliest_vs_latest_roadside_alien_advantage"]=(
           (item["late1"]-item["late3"])-(item["early1"]-item["early3"]))
        results.append(item)
    return pd.DataFrame(results)


def summarize_change(table):
    if table.empty:raise ValueError("No genuine repeated MIREN original 1/3 site-year source pairs")
    if table.duplicated(list(SITE)).any():
        raise ValueError("Same original repeated road-site used more than once")
    if table["years_between_first_and_last"].min()<=0:
        raise ValueError("No time separation between actual first/last MIREN field surveys")
    required=("gain_roadside_minus_interior",
              "loss_roadside_minus_interior",
              "difference_in_gain_fraction_roadside_minus_interior",
              "difference_in_earliest_vs_latest_roadside_alien_advantage")
    regions=[]
    for reg,grp in table.groupby("Region"):
        if len(grp)<MIN_REGION_REPEAT_SITES:continue
        details={
          "Region_code":reg,
          "n_original_repeated_road_sites":len(grp)}
        for k in required:
            y=grp[k].dropna()
            details["mean_"+k]=round(float(y.mean()),6) if len(y) else None
            details["n_eligible_"+k]=len(y)
            details["median_"+k]=round(float(y.median()),6) if len(y) else None
        regions.append(details)
    if len(regions)>MAX_PLOT_REGIONS_WITH_REPEATED_SOURCE:
        raise ValueError("More than predeclared 100 individual region codes unexpectedly used")
    def signs(metric):
        names=["mean_"+metric]
        vals=[r[names[0]] for r in regions if r[names[0]] is not None]
        return {
            "n_region_codes_with_>=10_repeated_sites_and_outcome":len(vals),
            "n_region_codes_mean_positive":sum(v>0 for v in vals),
            "n_region_codes_mean_negative":sum(v<0 for v in vals),
            "n_region_codes_mean_exact_zero":sum(v==0 for v in vals),
            "median_regional_average":round(float(np.median(vals)),6) if vals else None,
        }
    gain=table["gain_roadside_minus_interior"]
    loss=table["loss_roadside_minus_interior"]
    proportion=table["difference_in_gain_fraction_roadside_minus_interior"].dropna()
    return {
        "n_original_matched_first_to_last_road_sites":len(table),
        "n_original_matched_repeated_source_region_codes":table["Region"].nunique(),
        "n_original_matched_repeated_road_source_ids":len(table[["Region","Road"]].drop_duplicates()),
        "first_vs_last_year_gap_distribution":{
            str(k):int(v) for k,v in sorted(Counter(table["years_between_first_and_last"]).items())},
        "roadside_mean_gain_taxon_categories":round(float(table["gain1"].mean()),6),
        "interior_mean_gain_taxon_categories":round(float(table["gain3"].mean()),6),
        "roadside_mean_loss_taxon_categories":round(float(table["loss1"].mean()),6),
        "interior_mean_loss_taxon_categories":round(float(table["loss3"].mean()),6),
        "pooled_site_mean_gain_roadside_minus_interior":round(float(gain.mean()),6),
        "pooled_site_median_gain_roadside_minus_interior":round(float(gain.median()),6),
        "n_original_sites_more_gain_roadside":int((gain>0).sum()),
        "n_original_sites_equal_gain":int((gain==0).sum()),
        "n_original_sites_more_gain_interior":int((gain<0).sum()),
        "pooled_site_mean_loss_roadside_minus_interior":round(float(loss.mean()),6),
        "n_original_sites_with_nonzero_union_in_BOTH_plots_for_fraction":len(proportion),
        "pooled_site_mean_relative_gain_fraction_roadside_minus_interior":
            round(float(proportion.mean()),6) if len(proportion) else None,
        "pooled_site_median_relative_gain_fraction_roadside_minus_interior":
            round(float(proportion.median()),6) if len(proportion) else None,
        "region_code_equal_signed_mean_absolute_gain_source":
            signs("gain_roadside_minus_interior"),
        "region_code_equal_signed_mean_relative_gain_source":
            signs("difference_in_gain_fraction_roadside_minus_interior"),
        "region_code_equal_signed_mean_roadside_alien_advantage_change_source":
            signs("difference_in_earliest_vs_latest_roadside_alien_advantage"),
        "region_level_actual_source_outcome_table":regions,
        "no_statepair_pseudoreplication_or_local_seed_colonization_rate":True,
    }


def main():
    parser=argparse.ArgumentParser()
    for name in ("original-csv","verified-v99","output"):
        parser.add_argument("--"+name,required=True)
    a=parser.parse_args()
    original=Path(a.original_csv)
    if hashlib.md5(original.read_bytes()).hexdigest()!=SOURCE_MD5:
        raise ValueError("Original 2025 Zenodo MIREN field source md5 mismatch")
    data=pd.read_csv(original,low_memory=False)
    if len(data)!=SOURCE_ROWS:
        raise ValueError("MIREN original row population source drift")
    prior=json.loads(Path(a.verified_v99).read_text(encoding="utf-8"))
    if (prior.get("complete") is not True or
        prior.get("source_original_MIREN_Zenodo_DOI")!=DOI or
        prior.get("repeat_same_site_earliest_vs_latest_alien_roadside_advantage",{}).get(
            "n_original_road_sites_with_two_or_more_matched_survey_years")!=484):
        raise ValueError("Independently verified v9.9 ROAD site-year source provenance missing")
    clean,qc=native_status_safe(data)
    source_matched=matched_siteyears(plot_species_richness(clean))[1]
    if len(source_matched)!=prior["geography"][
        "n_matched_original_plot1_plot3_same_site_year_pairs"]:
        raise ValueError("Original strict Plot 1 vs 3 same site/year source plot matches changed")
    table=per_site_turnover(clean)
    outcome=summarize_change(table)
    if outcome["n_original_matched_first_to_last_road_sites"]!=484:
        raise ValueError("V9.9 original repeated MIREN road-site support changed")
    # Mandatory independently coded invariance from plot data:
    # net roadside change = (gain1-loss1)-(gain3-loss3).
    if not np.array_equal(
         table["difference_in_earliest_vs_latest_roadside_alien_advantage"].to_numpy(),
         (table["gain1"]-table["loss1"]-table["gain3"]+table["loss3"]).to_numpy()):
        raise ValueError("Observed MIREN species-set turnover does not conserve net richness change")
    report={
        "complete":True,
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_2025_Zenodo_MIREN_v3_md5":SOURCE_MD5,
        "original_MIREN_dataset_DOI":DOI,
        "actual_standardized_ROADside_original_plot_year_not_railway":True,
        "same_167804_original_source_rows_and_1876_paired_road_site_years_reconciled":True,
        "same_484_repeated_original_site_IDs_reconciled":True,
        "recorded_Alien_species_gains_NOT_proven_new_plant_colonizations":True,
        "results":outcome,
        "NO_actual_train_seed_dispersion_or_mowing_gust_heritable_evolution":True,
        "limits":[
            "A new taxon in an actually resurveyed MIREN original vegetation plot is a stronger local field-record change than a first GBIF state record, but it is NOT verified biological immigration; source imperfect plant detection remains.",
            "The result counts changed accepted names in strict original Plot1/3 equal-area MIREN road survey (2m×50m), NOT railway trackside surveys.",
            "Absolute gains scale naturally with the size of early/late species pools; separate GAIN / UNION controls that denominator but is still not colonization probability.",
            "If neither period in a plot contained an Alien taxon, gain/union is UNDEFINED, not zero; relative rates use only sites with nonempty Alien species pools in BOTH plots.",
            "Only 11 geographical region codes have >=10 repeated original road sites. These are not 484 independent mountain regions and not randomized traffic/stress assignments.",
            "Original 2007–08 MIREN surveys sometimes focused on nonnative plants, and plots destroyed or relocated are labeled with A suffix and excluded. The vegetation observation frame changes with protocol and time.",
            "No direct cut/mowing frequency, railway air gust measurement, biomass recovery, fitness, genetic heritability or cargo propagule load were collected in the original source.",
        ]
    }
    out=Path(a.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "MIREN_exact_roadside_revisit_result":{
            k:v for k,v in outcome.items() if k!="region_level_actual_source_outcome_table"},
        "no_railway_seed_or_stress_adaptation_causality":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
