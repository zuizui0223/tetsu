#!/usr/bin/env python3
"""v13.0: direct physically measured TRAIN-SPILLED GRAIN at 19 real
Canadian Pacific railway sites, published numeric source summaries.

Source ORIGINAL peer-reviewed field experiment: Gangadharan A,
Pollock S, Gilhooly P, Friesen A, Dorsey B, St. Clair CC (2017).
Grain spilled from moving trains create a substantial wildlife
attractant in protected areas. Animal Conservation 20:391–400.
DOI 10.1111/acv.12336.

AUTHORS directly collected and COUNTED six types of kernels spilled
by moving FREIGHT TRAINS at 19 permanent actual field sites along
134 km of Banff+Yoho Canadian Pacific railway 2013–2015, 1979
author samples, 2-week visits and one 2014–15 winter. Here tetsu
only reanalyzes the exact PUBLISHED source-level seasonal summary
means/SD & author's explicitly MODEL-EXTRAPOLATED totals; no original
site/visit raw CSV is publicly downloaded or invented.

BASKET original warm-season Apr-Oct 2013/2014/2015 author summary:
mean±SD gram/m²/day 1.12±1.23, .87±.71, 2.75±5.52, respectively;
author 134km track area projections 46.0,35.6,112.8 tons.
WINTER Nov2014–Mar2015 SCREEN 1.52±2.37 g/m²/day,
author projected 44.6 tons.
Author annual c.110 tons is PROJECTION 3-year mean
Apr–Oct 64.8t + one observed winter 44.6t ≈109.4t;
NOT tonnes personally collected by scientists.
All warm-season pooled samples 1.64±3.60 g/m²/day;
NOT same as simple 3-year 1.58g/m²/day average
because of original sampling/weight distribution.
Original grain composition 46% wheat, 25% CANOLA Brassica spp,
18% peas, 10% lentils, 1% canary,1% flax (sum 101% due
author rounding, must NOT renormalize to exactly 100).
Devices: 652 screen + 614 basket + 614 box + 99 lid = 1979;
year replicate original season rate is not identical to
per site/visit or per seed observations.
Original author models explained 9–31% variation.
True field seed viability/germination and subsequent seedling
establishment NOT measured; do not call 110 tonnes viable plants.

No claims of direct *wind* diaspore transport by passenger
train or grass cutting/local evolutionary adaptation. This study
is cargo grain drop, distinct mechanism from aerodynamic seeds.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter
from pathlib import Path
from statistics import mean,stdev

PAPER_DOI="10.1111/acv.12336"
ORIGINAL_PUBLISHER_FULLTEXT="https://onlinelibrary.wiley.com/doi/10.1111/acv.12336/full"
SPECIES_BY_GRAIN_MASS_SOURCE_PERCENT={
    "wheat_and_barley":46,
    "Brassica_canola":25,
    "peas_including_soybean":18,
    "lentils":10,
    "canary_seed":1,
    "flax":1,
}
ORIGINAL_COLLECTOR_SAMPLE_N={
    "screen":652,"basket":614,"box":614,"lid":99,
}
SOURCE_CSV={
    ("April-October",2013,2013,"basket"):(1.12,1.23,46.0),
    ("April-October",2014,2014,"basket"):(0.87,0.71,35.6),
    ("April-October",2015,2015,"basket"):(2.75,5.52,112.8),
    ("November-March",2014,2015,"screen"):(1.52,2.37,44.6),
}


def read_published_author_rail_kernels_summary(path):
    expected=["original_survey_season","source_calendar_start_year",
              "source_calendar_end_year","source_grain_collector",
              "original_mean_all_grain_deposition_g_m2_day",
              "original_sd_deposition_g_m2_day",
              "author_extrapolated_grain_mass_tonnes_for_134km"]
    rows={}
    with open(path,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=expected:
            raise ValueError("Published original train kernel source summary schema mismatch")
        for original in reader:
            k=(original["original_survey_season"],
               int(original["source_calendar_start_year"]),
               int(original["source_calendar_end_year"]),
               original["source_grain_collector"])
            if k in rows or k not in SOURCE_CSV:
                raise ValueError("Unexpected/duplicate original author train study year/season/method")
            values=tuple(float(original[col]) for col in expected[4:])
            if not all(math.isfinite(z) and z>=0 for z in values):
                raise ValueError("Original source train deposited grain numeric value invalid")
            if any(not math.isclose(z,actual,abs_tol=1e-11)
                   for z,actual in zip(values,SOURCE_CSV[k])):
                raise ValueError("Published original author Canadian train deposition summary changed "+str(k))
            rows[k]={
                "season":k[0],"start_year":k[1],"end_year":k[2],
                "original_collector":k[3],
                "author_reported_mean_g_per_m2_day":values[0],
                "author_reported_sd_g_per_m2_day":values[1],
                "author_extrapolated_tonnes_134km_NOT_physically_recovered":values[2],
            }
    if set(rows)!=set(SOURCE_CSV):
        raise ValueError("Missing published 2013–2015 warm and Nov2014–Mar2015 winter sample summaries")
    return rows


def actual_original_annual_source(rows):
    warm=[x for k,x in rows.items() if k[0]=="April-October"]
    winter=next(x for k,x in rows.items() if k[0]=="November-March")
    if len(warm)!=3:raise ValueError("Three source warm season repeats missing")
    annual_paper_extrap_warm=mean(x["author_extrapolated_tonnes_134km_NOT_physically_recovered"]
                                 for x in warm)
    annual_paper_extrap_total=annual_paper_extrap_warm+winter[
        "author_extrapolated_tonnes_134km_NOT_physically_recovered"]
    if not math.isclose(annual_paper_extrap_warm,64.8,abs_tol=1e-11):
        raise ValueError("Published author warm-season 3-year extrap mean 64.8t mismatch")
    if not math.isclose(annual_paper_extrap_total,109.4,abs_tol=1e-11):
        raise ValueError("Published author approximate annual 110t not reconciled")
    return {
      "n_original_real_railway_deposition_site_units":19,
      "n_original_track_distance_km":134,
      "original_sample_period_warm_seasons":3,
      "original_exclusive_observed_winter_seasons":1,
      "reported_original_field_sample_count":1979,
      "original_sampling_collector_method_counts":ORIGINAL_COLLECTOR_SAMPLE_N,
      "original_grain_type_source_reported_mass_shares_percent":
          SPECIES_BY_GRAIN_MASS_SOURCE_PERCENT,
      "author_rounded_grain_type_share_total_percent":
          sum(SPECIES_BY_GRAIN_MASS_SOURCE_PERCENT.values()),
      "original_all_warm_season_pooled_source_mean_g_m2_day":1.64,
      "original_all_warm_season_pooled_source_sd_g_m2_day":3.60,
      "original_winter_reported_mean_g_m2_day":1.52,
      "original_winter_reported_sd_g_m2_day":2.37,
      "source_unweighted_mean_of_3_reported_yearwise_warm_rates_NOT_POOLED":
          mean(x["author_reported_mean_g_per_m2_day"] for x in warm),
      "author_original_warm_yearly_deposition_mean_extrapolated_tonnes_134km":
          annual_paper_extrap_warm,
      "author_original_warm_yearly_deposition_SD_extrapolated_tonnes_134km":
          stdev(x["author_extrapolated_tonnes_134km_NOT_physically_recovered"]
                 for x in warm),
      "author_one_observed_winter_model_projection_tonnes_134km":
          winter["author_extrapolated_tonnes_134km_NOT_physically_recovered"],
      "mean_warm_projection_plus_one_observed_winter_projection_tonnes_134km":
          annual_paper_extrap_total,
      "author_rounded_annual_TOTAL_model_projection_tonnes_134km":110,
      "source_2015_warm_rate_divided_by_2014_warm_rate":
          SOURCE_CSV[("April-October",2015,2015,"basket")][0]/
          SOURCE_CSV[("April-October",2014,2014,"basket")][0],
      "source_2015_warm_masses_divided_by_2014_model_extrapolated_masses":
          112.8/35.6,
      "source_author_models_min_max_explained_deposition_variance_percent":[9,31],
      "source_18_sample_method_question_four_original_devices_2013_then_three":
          "Original four equipment designs evaluated 2013, only screen active winter; basket captured most during snow-free period. Do NOT compare winter screen rates with summer basket rates as randomized seasonal effect.",
      "directly_measured_train_spilled_cargo_grain_kernels":True,
      "directly_measured_grain_seed_viability_or_germination":False,
      "directly_measured_individual_spilled_grain_establishing_in_ballast":False,
      "source_observed_1979_samples_are_REPEATED_nested_within_original_19_sites":True,
    }


def run(path):
    original=read_published_author_rail_kernels_summary(path)
    yearly=actual_original_annual_source(original)
    return {
        "complete":True,"version":"v13.0",
        "source_peer_reviewed_original_article":
          "Gangadharan, Pollock, Gilhooly, Friesen, Dorsey & St. Clair 2017 Animal Conservation 20:391–400",
        "original_peer_reviewed_article_DOI":PAPER_DOI,
        "original_publisher_source_article_URL":ORIGINAL_PUBLISHER_FULLTEXT,
        "source_data_reuse_status":
            "manually transcribed true author's published seasonal summary statistics, NOT source raw 19site×visit seed counts or 1979 row matrix; original raw site-level data NOT acquired",
        "author_original_season_by_year_and_method_results":[
            original[k] for k in SOURCE_CSV],
        "original_author_real_train_GR AIN_deposition_source".replace(" ",""):
            yearly,
        "critical_limitations":[
            "Actual author field collectors counted spilled grain/cargo kernels physically recovered at 19 RAILWAY sites; viability of individual collected kernels/ germination / establishment has NOT been measured.",
            "Annual source approximately 110 tonnes is AUTHORS' extrapolation using study 134km railway geometry, not 110 tons physically collected; warm-year and winter rates derive from different equipment and observation windows.",
            "Original 1979 samples are repeated within 19 source sites at about 2-week sampling intervals; no independence of sample rows, date/site or annual effects can be assumed.",
            "The six grain type mass fractions 46,25,18,10,1,1 add to 101% due original source rounding, not a normalized six-part true mixture or species abundance.",
            "Canola (Brassica spp.) represents 25% author-reported deposited grain MASS, NOT 25% viable seeds germinating or population invasions along rail.",
            "The original 2017 goal concerns grizzly/wildlife food attractants, not a source plant colonization experiment; railway grade/cargo and management predictors remain confounded.",
            "No direct train slipstream effect on wind-dispersed noncrop flora, no genetic adaptation, flower phenology or rail mowing field data.",
            "Author-source original 19site×visit raw source data and appendix row-level observations NOT acquired, no spurious site-level mixed model or significance test fitted.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--published-2017-season-summary",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=run(a.published_2017_season_summary)
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v13.0",
        "seasonal_author_original_published_source":out[
            "author_original_season_by_year_and_method_results"],
        "direct_real_train_grain_source":out[
            "original_author_real_train_GRAIN_deposition_source"],
        "no_raw_19source_site_data_or_spilled_viability_observations_claimed":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()
