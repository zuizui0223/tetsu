#!/usr/bin/env python3
"""v15.9 independent source-study HOLDOUT influence audit from already VERIFIED
GitHub Actions original-source aggregate artifacts. NO new plant sampling.

Slovenia Šabić et al (2026) 89 100m original transect field sites/15 route
sections, original GPS, 2024 Geofabrik OSM rail+road source.
Czech Kutlvašr et al (2024) 35 original whole railway stations/5 stretches,
2026 SR70 CENTROID PROXIES not original plant GPS, 2023 Geofabrik OSM.
Both field-studies' original 2022-23 and 2020-22 respectively.
Nonindependent site pairs and fold errors, shared source national regions.

Outputs:
- exact reproduction of each source experiment's heldout equal-section
  average MSE and road-to-rail increment;
- one-original-source-section-HELDOUT-FOLD deletion sensitivity, same
  prefit error ledger, NOT new retraining;
- source-test-pair-weighted alternative aggregation, with duplicated
  partner sites and pairs across folds explicitly flagged;
- NEVER pool percentages/coefficients as common train seed dispersal meta.
The 2 observations remain from separate differently-grained field studies.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from statistics import mean

SLO_DOI="10.1007/s44473-026-00155-z"
CZ_DOI="10.23855/preslia.2024.247"
SLO_PBF_SHA="b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39"
CZ_PBF_SHA="26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"
SLO_V151_BENCHMARK=-0.5370955189438025
CZ_KNOWN_V156_BENCHMARK=-2.633035077679033
CZ_ALL_V156_BENCHMARK=-3.426664292963867

def check_close(x,y,label,atol=1e-10):
    if not(math.isfinite(x) and math.isfinite(y) and abs(x-y)<=atol):
        raise ValueError(label+" original verified author-source summary no longer reproduces")

def fold_summary(folds,road_field,rail_field,n_field,sec_field,
                 expected_n,expected_pct,name):
    if len(folds)!=expected_n or len({x[sec_field] for x in folds})!=expected_n:
        raise ValueError(name+" original route fold count/unique labels changed")
    if any(not x["identified"] for x in folds if "identified" in x):
        raise ValueError(name+" actual original section not fully estimable")
    rows=[]
    for x in folds:
        n=int(x[n_field])
        a=float(x[road_field]);b=float(x[rail_field])
        if n<1 or not all(math.isfinite(y) and y>=0 for y in (a,b)):
            raise ValueError(name+" invalid original source heldout pair number/MSE")
        rows.append({"section":str(x[sec_field]),"n_test_site_pairs":n,
                     "road_MSE":a,"rail_MSE":b,"rail_gain":a-b})
    road=mean(x["road_MSE"] for x in rows)
    rail=mean(x["rail_MSE"] for x in rows)
    pct=100*(road-rail)/road
    check_close(pct,expected_pct,name+" original overall rail gain")
    leave_fold=[]
    for omitted in rows:
        rest=[v for v in rows if v["section"]!=omitted["section"]]
        p=100*(mean(v["road_MSE"] for v in rest)-
               mean(v["rail_MSE"] for v in rest))/mean(v["road_MSE"] for v in rest)
        leave_fold.append({
            "excluded_original_heldout_source_section":omitted["section"],
            "rail_mse_gain_pct_from_other_original_fold_errors":p,
            "total_remaining_original_folds":len(rest)})
    weight=sum(x["n_test_site_pairs"] for x in rows)
    a=sum(x["road_MSE"]*x["n_test_site_pairs"] for x in rows)/weight
    b=sum(x["rail_MSE"]*x["n_test_site_pairs"] for x in rows)/weight
    return {
        "original_source_study":name,
        "original_source_route_section_folds":len(rows),
        "source_section_equal_weight_road_mse":road,
        "source_section_equal_weight_rail_mse":rail,
        "source_section_equal_weight_rail_increment_MSE_gain_pct":pct,
        "n_original_source_folds_with_rail_improvement":sum(v["rail_gain"]>0 for v in rows),
        "original_test_pair_weighted_road_mse":a,
        "original_test_pair_weighted_rail_mse":b,
        "source_test_pair_weighted_rail_increment_MSE_gain_pct":100*(a-b)/a,
        "source_original_test_pair_weights_sum_DUPLICATED_ACROSS_FOLDS":weight,
        "source_leave_one_PREDICTION_FOLD_OUT_influence":{
            "n_source_fold_deletions":len(leave_fold),
            "n_original_fold_deletions_with_positive_rail_increment":sum(
                t["rail_mse_gain_pct_from_other_original_fold_errors"]>0 for t in leave_fold),
            "minimum_rail_gain_pct":min(t["rail_mse_gain_pct_from_other_original_fold_errors"] for t in leave_fold),
            "maximum_rail_gain_pct":max(t["rail_mse_gain_pct_from_other_original_fold_errors"] for t in leave_fold),
            "individual_original_stretch_folds_removed":leave_fold,
            "NOT_original_retraining_after_removing_source_study_section":True,
        },
        "all_original_section_fold_error_rows":rows,
        "folds_and_pairs_NOT_independent_field_studies":True,
    }

def run(slo,cze):
    if slo.get("version")!="v15.1_real_topology_road_vs_rail_plant_identity":
        raise ValueError("Slovenia original v15.1 source version is not verified")
    if cze.get("version")!="v15.6_Czech2023_real_rail_road_conditional_station_flora_CV":
        raise ValueError("Czech original v15.6 source version is not verified")
    slgeo=slo["original_2026_field_and_2024_Geofabrik_source"]
    czgeo=cze["original_physical_network_provenance"]
    if slgeo["original_vegetation_DOI"]!=SLO_DOI or czgeo["original_botanical_source_DOI"]!=CZ_DOI:
        raise ValueError("Original two published cohort DOIs changed")
    if (slgeo["original_OSM_pbf_SHA256"]!=SLO_PBF_SHA or
        czgeo["actual_download_source_SHA256"]!=CZ_PBF_SHA):
        raise ValueError("Original source date-specific Geofabrik rail/road PBF SHA changed")
    slmain=slo["main_section_heldout_rail_vs_road_predictive_result"]
    sl=fold_summary(slmain["details_per_original_section"],
        "mse_plus_road","mse_plus_road_and_rail","n_holdout","source_section",
        15,SLO_V151_BENCHMARK,"Slovenia_original_89_100m_transects")
    checks={
       "known":("both_Pladias_known",CZ_KNOWN_V156_BENCHMARK),
       "all":("all_original_S1_current",CZ_ALL_V156_BENCHMARK)
    }
    cz={}
    for name,(key,expected) in checks.items():
        original=cze["source_frozen_geographic_guild_control_inference"][key]
        main=original["original_stretch_heldout_CV"]
        c=fold_summary(main["source_stretch_holdout"],
          "MSE_M1_plus_actual_2023_road_route",
          "MSE_M2_plus_actual_2023_road_and_rail_ratio",
          "n_test_pairs","heldout_original_source_stretch",5,expected,
          "Czech_original_35_whole_stations_"+name)
        if original["source_original_railroad_pair_support"]!=594:
            raise ValueError("Czech authentic 594 pair support drifted")
        cz[name]=c
    return {
       "version":"v15.9_authentic_two_original_studies_distinct_source_validation_influence",
       "original_source_DOIs":[SLO_DOI,CZ_DOI],
       "true_original_date_fixed_geofabrik_SHA":[SLO_PBF_SHA,CZ_PBF_SHA],
       "source_Slovenia_v151":sl,
       "source_Czech_v156":cz,
       "two_source_cohorts_NOT_common_meta_effect":True,
       "only_within_source_holdout_error_influence_NOT_permutation_inference":True,
       "research_limitations":[
          "Both original studies have non-independent station/site pairs and overlapping partner sites across heldout folds; no confidence intervals pretending each heldout pair is independent.",
          "The five Czech source stretches are nested within ONE country; 15 Slovenian source sections also belong to ONE country.",
          "Czech whole-station survey and 2026 centroid proxy vs Slovenian true original 100m field GPS/2024 network differ in spatial grain, taxa/guild ontology, sample design and metadata.",
          "Source prediction-section fold deletion removes only an evaluation fold from the aggregate; it does not re-estimate each omitted route's 4 or 14-section training models.",
          "Any change from negative aggregate to positive after one fold deletion signals fragility, NOT the existence of railway-mediated seed dispersal.",
          "A geometric railroad route is not measured active freight movement, seed carriage, cutting or evolved wind/flowering traits.",
          "Study-specific predictive MSE improvement percentages may be discussed side by side as methods diagnostics, never pooled into a common biological effect."
       ],
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-slovenia-v151",required=True)
    p.add_argument("--source-czech-v156",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=run(
      json.loads(Path(a.source_slovenia_v151).read_text(encoding="utf-8")),
      json.loads(Path(a.source_czech_v156).read_text(encoding="utf-8")))
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    brief={}
    for name,d in (("Slovenia",result["source_Slovenia_v151"]),
                   ("Czech_known",result["source_Czech_v156"]["known"]),
                   ("Czech_all",result["source_Czech_v156"]["all"])):
        brief[name]={
            "source_original_route_folds":d["original_source_route_section_folds"],
            "rail_gain_equal_section_pct":d["source_section_equal_weight_rail_increment_MSE_gain_pct"],
            "rail_gain_test_pair_weighted_pct":d["source_test_pair_weighted_rail_increment_MSE_gain_pct"],
            "n_positive_route_folds":d["n_original_source_folds_with_rail_improvement"],
            "original_fold_deletion_sensitivity":d["source_leave_one_PREDICTION_FOLD_OUT_influence"]
        }
    print(json.dumps({"complete":True,"version":result["version"],"true_source_influence":brief,
                     "not_crosscountry_pooled_meta_or_causal_seed_dispersal":True},ensure_ascii=False))

if __name__=="__main__":main()
