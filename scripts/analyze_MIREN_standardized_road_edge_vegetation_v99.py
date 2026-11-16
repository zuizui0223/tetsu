#!/usr/bin/env python3
"""v9.9 independently surveyed roadside-vs-interior plant flora macro.

Original actual MIREN v3 public SOURCE CSV (2007–23), Zenodo
DOI 10.5281/zenodo.15749536; published as Haider et al. (2026)
Global Ecology and Biogeography DOI 10.1111/geb.70284.

The independent survey design at a fixed mountain ROAD site has
three 2m×50m plots: Plot 1 roadside, Plot 2 2–52m, Plot 3 52–102m.
Use same Region+Road+Site+Year pairing before estimating any
alien-plant richness difference. Never compare GBIF iNaturalist
positive-only sites as biological absences.

Pre-original data integrity:
- Source MD5 must match original Zenodo b640195cc44f252ccb6ddb28eb6fff13.
- The original file has 167804 CSV occurrence rows; may contain
  duplicate species per identical plot-year. Deduplicate taxa.
- Exactly 2 region-species show contradictory Native/Alien status:
  conservatively classify them UNKNOWN for all plots rather than
  choosing favorable nativity.
- Strict original Plot '1','2','3' main result excludes replacement
  identifiers '1A','2A','3A'. They are reported as source omissions.
- Missing plot-year from original positive-only source is NOT a
  zero-species surveyed plot. Only complete same site-year paired
  plots are compared; missing region/plot groups are not zeros.
- First survey year for each road-site is the primary independent
  local site sensitivity, so five-year resurvey repeat isn't
  pseudoreplicated. Region codes are not randomized regions.

This is a ROAD disturbance/establishment-filter observational BASELINE
for future train-versus-road comparisons. NOT observed train-mediated
seed transfer, measured cutting/gust, locally evolved adaptation.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

DOI="10.5281/zenodo.15749536"
PAPER_DOI="10.1111/geb.70284"
SOURCE_MD5="b640195cc44f252ccb6ddb28eb6fff13"
SOURCE_ROWS=167804
REQUIRED=("Plot_id","Region","Road","Site","Plot","Year",
          "Elevation","Latitude","Longitude","Accepted.Name.MIREN",
          "Abundance","Cover","Status","Protocol")
SITE_YEAR=("Region","Road","Site","Year")
SITE=("Region","Road","Site")
STRICT_PLOTS=("1","2","3")
GROUPS=("Alien","Native","Unknown")
MIN_REGION_PAIRED_SITEYEARS=20
MIN_REGION_FIRST_SURVEY_SITES=20
MIN_TEMPORAL_REPEAT_SITES_PER_REGION=10
REGION_NOT_INDEPENDENT_TREATMENTS=True


def source_integrity(source_csv,source_audit):
    raw=Path(source_csv).read_bytes()
    md5=hashlib.md5(raw).hexdigest()
    if md5!=SOURCE_MD5:
        raise ValueError("MIREN source original MD5 does not match frozen Zenodo v3")
    audit=json.loads(Path(source_audit).read_text(encoding="utf-8"))
    declared=audit["source_original_binary_metadata"]
    if (audit.get("complete") is not True or
        declared.get("source_original_zenodo_md5")!=md5 or
        declared.get("source_original_Zenodo_DOI")!=DOI or
        audit["actual_CSV_source_schema"]["n_original_csv_rows"]!=SOURCE_ROWS):
        raise ValueError("Independent v9.8 original Zenodo original CSV source audit not reconciled")
    source=pd.read_csv(source_csv,low_memory=False)
    if len(source)!=SOURCE_ROWS or set(REQUIRED)-set(source.columns):
        raise ValueError("2025 MIREN original source rows/columns changed")
    return source,md5


def native_status_safe(source):
    """No double-counting the same species recorded twice at a plot-year;
    region-level contradictory nativity never decided by outcome.
    """
    src=source[source["Plot"].isin(STRICT_PLOTS) &
               source["Accepted.Name.MIREN"].notna()].copy()
    src["Status"]=src["Status"].fillna("Unknown").replace({
        "":"Unknown","NA":"Unknown"})
    if not set(src["Status"]).issubset(set(GROUPS)):
        raise ValueError("Original regional species nativity has new unknown labels")
    conflict=(src.groupby(["Region","Accepted.Name.MIREN"])["Status"]
                .nunique(dropna=True))
    conflicting=set(conflict[conflict>1].index)
    mask=[(region,name) in conflicting for region,name in zip(
        src["Region"],src["Accepted.Name.MIREN"])]
    src.loc[mask,"Status"]="Unknown"
    before=len(src)
    src=src.drop_duplicates(list(SITE_YEAR)+["Plot","Accepted.Name.MIREN"])
    if src.duplicated(list(SITE_YEAR)+["Plot","Accepted.Name.MIREN"]).any():
        raise ValueError("Same taxon counted multiple times within exact plot-year")
    return src,{
        "raw_original_source_rows":len(source),
        "n_original_STRICT_plot_1_2_3_rows_before_dedup":before,
        "n_removed_repeated_species_same_plot_year":before-len(src),
        "n_region_species_with_ambiguous_native_alien_status_excluded_from_both":
            len(conflicting),
        "n_status_unknown_records_after_dedup":int((src["Status"]=="Unknown").sum()),
        "n_original_nonstandard_A_replacement_plot_rows_excluded":
            int((~source["Plot"].isin(STRICT_PLOTS)).sum()),
        "n_original_region_codes_with_any_plot_rows":int(source["Region"].nunique()),
        "strict_plot_1_2_3_no_replacement_assumptions":True,
        "unobserved_plots_not_filled_with_zero_species":True,
    }


def plot_species_richness(src):
    groups=(src.groupby(list(SITE_YEAR)+["Plot","Status"])
              .size().unstack("Status",fill_value=0).reset_index())
    for kind in GROUPS:
        if kind not in groups:groups[kind]=0
    groups["total"]=groups[list(GROUPS)].sum(axis=1)
    if (groups[list(GROUPS)]<0).any().any():
        raise ValueError("Invalid independent repeated-plot richness source")
    return groups


def matched_siteyears(groups):
    wide=groups.pivot(index=list(SITE_YEAR),columns="Plot",
                      values=["Alien","Native","Unknown","total"])
    wide.columns=[str(a)+"_"+str(b) for a,b in wide.columns]
    wide=wide.reset_index()
    pairs=wide.dropna(subset=["Alien_1","Alien_3"]).copy()
    for name in ("Alien","Native","total"):
        pairs[name+"_roadside_minus_interior"]=pairs[name+"_1"]-pairs[name+"_3"]
    for plot in ("1","3"):
        d=pairs["Alien_"+plot]+pairs["Native_"+plot]
        pairs["alien_fraction_"+plot]=pairs["Alien_"+plot]/d.replace(0,np.nan)
    pairs["alien_fraction_roadside_minus_interior"]=(
        pairs["alien_fraction_1"]-pairs["alien_fraction_3"])
    triplets=wide.dropna(subset=["Alien_1","Alien_2","Alien_3"]).copy()
    return wide,pairs,triplets


def region_summary(data):
    if len(data)==0:return []
    entries=[]
    for code,chunk in data.groupby("Region"):
        alien=chunk["Alien_roadside_minus_interior"]
        native=chunk["Native_roadside_minus_interior"]
        frac=chunk["alien_fraction_roadside_minus_interior"]
        entries.append({
            "Region_code":code,
            "n_matched_site_years_or_unique_sites":len(chunk),
            "n_distinct_road_codes":int(chunk["Road"].nunique()),
            "n_distinct_road_site_IDs":len(chunk[list(SITE)].drop_duplicates()),
            "mean_roadside_minus_interior_alien_richness":round(float(alien.mean()),6),
            "median_roadside_minus_interior_alien_richness":round(float(alien.median()),6),
            "n_siteyears_roadside_alien_richer":int((alien>0).sum()),
            "n_siteyears_alien_equal":int((alien==0).sum()),
            "n_siteyears_interior_alien_richer":int((alien<0).sum()),
            "mean_roadside_minus_interior_native_richness":round(float(native.mean()),6),
            "mean_roadside_minus_interior_alien_species_fraction":
                round(float(frac.mean()),6) if frac.notna().any() else None,
            "region_mountain_road_source_not_independent_rail_treatment":True,
        })
    return sorted(entries,key=lambda r:r["Region_code"])


def region_signs(regions,min_matched):
    qualified=[z for z in regions if
               z["n_matched_site_years_or_unique_sites"]>=min_matched]
    vals=[z["mean_roadside_minus_interior_alien_richness"] for z in qualified]
    return {
        "min_original_matched_source_pairs_per_region":min_matched,
        "n_qualified_region_codes_not_independent_countries":len(qualified),
        "n_codes_roadside_alien_richness_higher":sum(v>0 for v in vals),
        "n_codes_equal":sum(v==0 for v in vals),
        "n_codes_interior_alien_richness_higher":sum(v<0 for v in vals),
        "median_region_equal_mean_alien_richness_delta":(
            round(float(np.median(vals)),6) if vals else None),
        "mean_region_equal_mean_alien_richness_delta":(
            round(float(np.mean(vals)),6) if vals else None),
        "NO_pseudoreplicated_siteyear_p_values_or_rail_causes":True,
    }


def repeated_site_temporal(pairs):
    ordered=pairs.sort_values("Year").copy()
    n=ordered.groupby(list(SITE)).size()
    eligible=n[n>=2].index
    if len(eligible)==0:return {"n_repeat_sites":0}
    repeated=(ordered.set_index(list(SITE)).loc[eligible].reset_index())
    first=repeated.groupby(list(SITE),as_index=False).first()
    last=repeated.groupby(list(SITE),as_index=False).last()
    paired=first.merge(last,on=list(SITE),suffixes=("_earliest","_latest"),
                       validate="one_to_one")
    paired["delta_roadside_alien_advantage"]=(
        paired["Alien_roadside_minus_interior_latest"]
       -paired["Alien_roadside_minus_interior_earliest"])
    result=[]
    for region,sub in paired.groupby("Region"):
        if len(sub)<MIN_TEMPORAL_REPEAT_SITES_PER_REGION:continue
        diff=sub["delta_roadside_alien_advantage"]
        result.append({
            "Region_code":region,
            "n_repeated_original_road_sites":len(sub),
            "mean_change_alien_roadside_advantage":round(float(diff.mean()),6),
            "n_sites_increased":int((diff>0).sum()),
            "n_sites_equal":int((diff==0).sum()),
            "n_sites_decreased":int((diff<0).sum()),
        })
    return {
        "n_original_road_sites_with_two_or_more_matched_survey_years":len(paired),
        "n_region_codes_with_at_least_10_repeated_sites":len(result),
        "n_region_codes_increasing_mean_roadside_alien_advantage":
            sum(r["mean_change_alien_roadside_advantage"]>0 for r in result),
        "n_region_codes_decreasing_mean_roadside_alien_advantage":
            sum(r["mean_change_alien_roadside_advantage"]<0 for r in result),
        "median_original_road_site_change_latest_minus_first":
            round(float(paired["delta_roadside_alien_advantage"].median()),6),
        "n_site_increasing":int((paired["delta_roadside_alien_advantage"]>0).sum()),
        "n_site_decreasing":int((paired["delta_roadside_alien_advantage"]<0).sum()),
        "regions":result,
        "earlier_protocols_native_survey_incomplete_temporal_alien_focus_only":True,
        "NO_road_management_or_railway_development_causal_trend":True,
    }


def compute(original):
    dedup,qc=native_status_safe(original)
    plot_table=plot_species_richness(dedup)
    wide,pairs,triplets=matched_siteyears(plot_table)
    primary_first_site=(pairs.sort_values("Year")
                         .drop_duplicates(list(SITE),keep="first"))
    last_site=(pairs.sort_values("Year")
               .drop_duplicates(list(SITE),keep="last"))
    region_all=region_summary(pairs)
    region_first=region_summary(primary_first_site)
    region_last=region_summary(last_site)
    region_2017=region_summary(pairs[pairs["Year"]>=2017])
    region_2022=region_summary(pairs[pairs["Year"]>=2022])
    if len(pairs)==0 or len(primary_first_site)==0:
        raise ValueError("No independent ROAD plot-1 / plot-3 matched vegetation source")
    median_triplet={f"plot_{j}":round(float(triplets[f"Alien_{j}"].median()),6)
                    for j in STRICT_PLOTS}
    report={
        "complete":True,
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_original_MIREN_Zenodo_DOI":DOI,
        "source_2026_GEB_data_paper_DOI":PAPER_DOI,
        "data_type":"actual standardized mountain ROAD flora original fixed 2m×50m vegetation plots, NOT railway vegetation",
        "source":qc,
        "geography":{
            "n_original_region_codes":int(original["Region"].nunique()),
            "n_region_codes_with_original_road_plot1_and_interior_plot3_same_site_year":
                len(region_all),
            "region_codes_without_matched_Plot1_Plot3":sorted(set(original["Region"])-set(
                r["Region_code"] for r in region_all)),
            "n_matched_original_plot1_plot3_same_site_year_pairs":len(pairs),
            "n_distinct_site_IDs_matching_both_plots":len(primary_first_site),
            "n_distinct_region_plus_road_code_sources":len(
                pairs[["Region","Road"]].drop_duplicates()),
            "n_observed_original_STRICT_plot_year_source_rows":len(plot_table),
            "n_original_STRICT_complete_three_plot_site_year_sets":len(triplets),
        },
        "primary_first_survey_per_original_site":{
            "n_unique_original_region_road_site_units":len(primary_first_site),
            "region_equal_alien_richness_roadside_minus_interior_source":
                region_signs(region_first,MIN_REGION_FIRST_SURVEY_SITES),
            "region_codes_all_matched_first_surveys":region_signs(region_first,1),
            "regions":region_first,
            "NO_repeated_5year_site_visits_counted_as_independent_sites":True,
        },
        "all_matched_siteyear_source_sensitivity":{
            "n_matched_same_site_and_year_Plot1_vs_Plot3":len(pairs),
            "mean_alien_richness_plot1_minus_plot3":round(float(
                pairs["Alien_roadside_minus_interior"].mean()),6),
            "median_alien_richness_plot1_minus_plot3":round(float(
                pairs["Alien_roadside_minus_interior"].median()),6),
            "n_siteyears_roadside_alien_more":int((pairs[
                "Alien_roadside_minus_interior"]>0).sum()),
            "n_siteyears_same_alien_richness":int((pairs[
                "Alien_roadside_minus_interior"]==0).sum()),
            "n_siteyears_interior_alien_more":int((pairs[
                "Alien_roadside_minus_interior"]<0).sum()),
            "mean_native_richness_plot1_minus_plot3":round(float(
                pairs["Native_roadside_minus_interior"].mean()),6),
            "regions_matched_source":region_signs(region_all,MIN_REGION_PAIRED_SITEYEARS),
            "regions_all_source":region_signs(region_all,1),
        },
        "complete_plot1_plot2_plot3_same_siteyear_source":{
            "n_complete_three_plot_siteyears":len(triplets),
            "median_alien_species_richness_per_100m2_plot":median_triplet,
            "mean_alien_richness_plot1_minus_plot2":round(float(
                (triplets["Alien_1"]-triplets["Alien_2"]).mean()),6),
            "mean_alien_richness_plot2_minus_plot3":round(float(
                (triplets["Alien_2"]-triplets["Alien_3"]).mean()),6),
            "mean_alien_richness_plot1_minus_plot3":round(float(
                (triplets["Alien_1"]-triplets["Alien_3"]).mean()),6),
        },
        "mandatory_period_and_first_last_sensitivity":{
            "first_original_source_per_site":region_signs(region_first,10),
            "last_original_source_per_site":region_signs(region_last,10),
            "original_year_at_least_2017":region_signs(region_2017,10),
            "original_year_at_least_2022":region_signs(region_2022,10),
            "region_2017_first_source_codes":region_2017,
            "region_2022_last_source_codes":region_2022,
        },
        "repeat_same_site_earliest_vs_latest_alien_roadside_advantage":
            repeated_site_temporal(pairs),
        "region_all_years":region_all,
        "NO_train_seed_dispersal_or_mowing_wind_heritable_adaptation_conclusions":True,
        "limitations":[
            "This is authentic standardized ROADside vegetation, not railway rail-yard or trackside surveys. It establishes no comparative RAIL-vs-ROAD causal effect.",
            "MIREN survey sites and mountain roads were not randomly assigned traffic, mowing, slipstream gust, gravel, soil or climate; paired plot differences are observational associations, not isolated disturbance mechanisms.",
            "Plot1 roadside, Plot2 2–52m and Plot3 52–102m are 100 m² equal-area patches. Their local vegetation, light and habitat differ naturally.",
            "Plot-year with zero reported plant species may not appear in long format; absent plot is NOT imputed as zero alien/native species.",
            "Original Zenodo v3 raw CSV has duplicate plot-year taxon lines, which must be deduplicated before biological species richness calculations; ambiguous region nativity keys go UNKNOWN.",
            "The 26 raw region CODES are not automatically 26 statistically independent countries; two have no strict original plot1/plot3 same-site-year pair. 5-year repeated sites must not be pseudoindependent.",
            "Early 2007/08 MIREN protocols often inventoried nonnative flora more completely than native flora; native-share analyses from all years are not a verified native-richness experimental control.",
            "The source Zenodo dataset and 2026 data paper are published original consortium survey data, NOT novel data collected by this GitHub project. Any manuscript must cite the consortium and comply with its data sharing/authorship guidelines.",
            "This dataset contains NO measured direct mowing frequency, field railway air gusts, seed-carriage on trains, plasticity or evolutionary genetic adaptation.",
        ],
    }
    return report


def main():
    p=argparse.ArgumentParser()
    for x in ("source-csv","source-audit-v98","output"):p.add_argument("--"+x,required=True)
    a=p.parse_args()
    raw,md5=source_integrity(a.source_csv,a.source_audit_v98)
    res=compute(raw)
    res["source"]["independently_reverified_original_Zenodo_MD5"]=md5
    path=Path(a.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,"source_original_rows":res["source"]["raw_original_source_rows"],
        "geography":res["geography"],
        "first_survey_primary":res["primary_first_survey_per_original_site"][
            "region_equal_alien_richness_roadside_minus_interior_source"],
        "all_same_siteyear":res["all_matched_siteyear_source_sensitivity"],
        "plot1_2_3":res["complete_plot1_plot2_plot3_same_siteyear_source"],
        "last_vs_first_original_repeated_sites":{
            k:v for k,v in res["repeat_same_site_earliest_vs_latest_alien_roadside_advantage"].items()
            if k!="regions"},
        "strict_source_no_actual_railway_transport_or_adaptation_effect":True,
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
