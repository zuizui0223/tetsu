#!/usr/bin/env python3
"""v25.1 independent methodological evidence gate for Nov 17 genuine tetsu
country-held-out railway flora results. NO NEW BIOLOGICAL OBSERVATIONS.

Audit the November-17 source-generated (not manually invented) repository
JSON under research/validation_20261117, source Git commit:
a8138a50c7db3abad7ea2f717dc25f374a34e87c.
Test whether a foreign-country station-frequency model is not only
associated with same-country railway taxa on a SELECTED shared species
support, but (1) beats the domestic disturbed-vegetation control at
species rank and (2) absolute-frequency calibration and (3) retains
the ranking advantage with unrecorded species and (4) leave-one-training
source out. The eight-gate designation is an EXPLORATORY RECONCILIATION,
NOT preregistered significance or evidence of causal railway fidelity.

Source evidence caveats:
 - station-wide railway lists (CZ, IT, HR, BA, SK) vs small surveyed
   roadside plots, tree pits and urban disturbed relevés: unequal
   sampling area and effort; controls not comparable or randomly drawn.
 - Slovakia 30 stations are rounded ORIGINAL 2025 station frequencies;
   independent 2016 East Slovakia 19 stations are a separate study.
 - HR and BA are ONE original paper; neither is independent of the
   same paper's other-country railway source.
 - Railway 2020/22, source trackside Slovenia 89 original 100m plots,
   Kashmir 2014 and 2017 same 16 stations, AU/CN only 3 corridors each.
 - genus bootstrap is botanical species resampling, not a country-
   level 95% confidence interval.
 - global sPlotOpen site distance to GeoNames station does not verify
   actual station footprint or trackbed; 1km and 2km cohorts are NOT
   nested, and different CWM definitions (TRY abundance-weighted
   plant height vs LEDA species-equal canopy height) must not pool.

Outputs only compact numerical stress audit with hard source Git blob
SHA gates. Does not redistribute raw source species-by-plot data.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from statistics import mean

SOURCE_COMMIT="a8138a50c7db3abad7ea2f717dc25f374a34e87c"
EXPECTED_BLOBS={
    "country_validation.json":"df8b84ba40c5eb6ee70536b188e2a4a5cfabeb9a",
    "global_station_proximity.json":"7001abd56d182e37187ae391d1ceba2d2c472ada",
    "slovakia_results.json":"3857854785298cee2fb23cdd0a707d7210d44e1d",
    "kashmir_results.json":"75dd78dd2b66a67543c248ff686c3a9dc016ea73",
    "replication_audit.json":"88864debadf9f8c3d41dc937c3fb79a426f70856",
}
MIN_COMMON_TAXA=50
MIN_TARGET_RECORDED_FRACTION=0.5

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+b"\x00"+raw).hexdigest()

def audited_source(root):
    out={}
    checks={}
    for filename,expected in EXPECTED_BLOBS.items():
        f=Path(root)/filename
        raw=f.read_bytes()
        actual=git_blob(raw)
        if actual!=expected:
            raise ValueError(f"Source original Nov-17 evidence JSON drifted: {filename}: {actual} != {expected}")
        out[filename]=json.loads(raw)
        checks[filename]={"git_blob_sha1":actual,"bytes":len(raw)}
    return out,checks

def evaluate_comparison(tag,row):
    common=row["common"];broad=row["extended"]
    leave=row["delta_excluding_training_study_fixed_panel"]
    if not leave or not isinstance(leave,dict):
        raise ValueError("Country source lacks leave-one-training-original-study prediction sensitivity")
    if common["n"]<common["recorded"] or broad["n"]<broad["recorded"]:
        raise ValueError("Recorded species count exceeds source evaluated panel")
    if broad["n"]<common["n"]:
        raise ValueError("Broader holdout source panel cannot be smaller than the common supported panel")
    if common["n"]<=0:raise ValueError("Zero source species support")
    ci=common["genus_interval"]
    if not isinstance(ci,list) or len(ci)!=2 or ci[0]>ci[1]:
        raise ValueError("Source original genus bootstrap quantile interval missing")
    if abs(common["delta"]-(common["foreign"]-common["domestic"]))>1e-12:
        raise ValueError("Original restricted model rank delta inconsistent")
    if abs(broad["delta"]-(broad["foreign"]-broad["domestic"]))>1e-12:
        raise ValueError("Original broad model rank delta inconsistent")
    criteria={
        "at_least_50_common_source_species":common["n"]>=MIN_COMMON_TAXA,
        "at_least_half_source_species_recorded_in_holdout":(
             common["recorded"]/common["n"]>=MIN_TARGET_RECORDED_FRACTION),
        "foreign_rank_beats_domestic_in_selected_common_taxa":common["delta"]>0,
        "genus_resampled_2point5pct_rank_delta_positive":ci[0]>0,
        "foreign_rank_beats_domestic_in_extended_support":broad["delta"]>0,
        "foreign_absolute_frequency_MAE_better_in_common_taxa":(
            common["mae_foreign"]<common["mae_domestic"]),
        "foreign_absolute_frequency_MAE_better_in_extended_support":(
            broad["mae_foreign"]<broad["mae_domestic"]),
        "rank_advantage_after_each_single_training_paper_removed":all(v>0 for v in leave.values()),
    }
    return {
        "target_reference_comparison":tag,
        "original_test_panel_shared_taxa":common["n"],
        "original_target_taxa_recorded_in_selected_panel":common["recorded"],
        "original_recorded_support_fraction":common["recorded"]/common["n"],
        "foreign_minus_domestic_rank_rho_common":common["delta"],
        "foreign_minus_domestic_rank_rho_extended":broad["delta"],
        "genus_bootstrap_source_lower_quantile":ci[0],
        "MAE_gain_domestic_minus_foreign_common":common["mae_domestic"]-common["mae_foreign"],
        "MAE_gain_domestic_minus_foreign_extended":broad["mae_domestic"]-broad["mae_foreign"],
        "worst_foreign_rank_advantage_leave_one_training_study_out":min(leave.values()),
        "n_source_original_training_paper_removal_scenarios":len(leave),
        "eight_exploratory_descriptive_reliability_conditions":criteria,
        "all_eight_conditions_met":all(criteria.values()),
        "not_a_preregistered_hypothesis_or_independent_biological_p_value":True,
    }

def station_proximity(station):
    if station.get("source_plots")!=95104 or station.get("source_station_records")!=66856:
        raise ValueError("Global independent sPlotOpen / GeoNames source record support changed")
    out={}
    for band in ("0.25","1.0","2.0","1.0_strict"):
        record=station["bands_km"][band]
        named=record["supported"]
        cm=[r["cwm"] for r in named]
        leda=[r["leda"]["equal"] for r in named if r["leda"].get("n",0)>=4]
        if not cm:
            out[band]={"n_countries_with_at_least_four_true_station_nearby_ref_pairs":0,
                       "no_cwm_estimate_not_zero_effect":True,
                       "n_original_station_georegistry_matches":record["matched_stations"]}
        else:
            out[band]={
               "n_contributing_original_country_datasets":len(named),
               "n_original_matched_station_plot_pairs":sum(r["stations"] for r in named),
               "TRY_cover_weighted_plantheight_seedmass_near_station_minus_other_mean_equal_country":
                   mean(cm),
               "n_country_datasets_positive_CWM_relative_distance":sum(v>0 for v in cm),
               "n_country_datasets_negative_CWM_relative_distance":sum(v<0 for v in cm),
               "country_CWM_sign_and_value":{x["country"]:x["cwm"] for x in named},
               "LEDA_species_equal_canopyheight_seedmass_n_contributing_datasets":len(leda),
               "LEDA_species_equal_canopyheight_seedmass_mean_on_subset_not_same_plots":
                   mean(leda) if leda else None,
               "nearby_station_plot_not_demonstrated_in_actual_trackbed":True}
    if abs(out["1.0"]["TRY_cover_weighted_plantheight_seedmass_near_station_minus_other_mean_equal_country"]-.2221)>.0002:
        raise ValueError("Nov-17 global station nearest original 1km CWM source point estimate drift")
    return out

def analyze(source):
    c=source["country_validation.json"]
    if len(c["comparisons"])!=10 or c["sites"].get("IN")!=16 or c["sites"].get("CZ")!=39:
        raise ValueError("Genuine original country-held-out count or source station cohort drift")
    evidence={k:evaluate_comparison(k,row) for k,row in sorted(c["comparisons"].items())}
    strong=[key for key,v in evidence.items() if v["all_eight_conditions_met"]]
    if "IT_road" not in evidence or not strong:
        raise ValueError("Expected genuine published country comparison results were not evaluated")
    typed={}
    for key,val in evidence.items():
        category=key.rsplit("_",1)[-1]
        typed.setdefault(category,[]).append(key)
    cz={}
    for b in ("road","urban"):
        row=c["comparisons"]["CZ_"+b].get("same_species_both_domestic_controls")
        if row is None or row["n"]!=52:raise ValueError("Exact 52-species matched Czech controls missing")
        cz[b]={"n_common_source_species":row["n"],
               "foreign_predicted_rho":row["foreign"],
               "domestic_control_predicted_rho":row["domestic"],
               "rank_advantage_foreign_minus_domestic":row["delta"],
               "genus_bootstrap_rank_advantage_interval":row["genus_interval"]}
    sk=source["slovakia_results.json"]
    if sk["source_taxa"]!=318 or sk["coverage"]["SK"]["source_total_author_station_n"]!=19:
        raise ValueError("Genuine independent Slovakia 19 actual station trait field source changed")
    m=sk["results"]["functional_CWM_centroid_distance_all_stationpairs"]
    if abs(m["observed_minus_null"]-(-.0023216720504834754))>1e-10:
        raise ValueError("New independent 19 Slovakia railway station functional beta null finding drifted")
    k=source["kashmir_results.json"]
    if k["source_taxa"]!=123 or k["both_traits"]!=76 or k["source_spatial_sites"]!=31:
        raise ValueError("New authentic Dryad original railway 31 sites with 76 traits changed")
    summary={
       "country_reference_comparison_count":len(evidence),
       "n_original_country_groups_including_two_domestic_controls_for_Czech":
           len({x.rsplit("_",1)[0] for x in evidence}),
       "n_common_panel_foreign_rank_more_predictive":sum(v["foreign_minus_domestic_rank_rho_common"]>0 for v in evidence.values()),
       "n_original_country_comparisons_genus_resample_lower_quantile_positive":sum(
           v["genus_bootstrap_source_lower_quantile"]>0 for v in evidence.values()),
       "n_common_and_extended_frequency_absolute_error_better":sum(
           v["MAE_gain_domestic_minus_foreign_common"]>0 and
           v["MAE_gain_domestic_minus_foreign_extended"]>0 for v in evidence.values()),
       "n_meeting_all_eight_EXPLORATORY_reliability_gates":len(strong),
       "country_reference_comparisons_meeting_all_eight":strong,
       "critical_exploratory_not_pre_registered":True,
    }
    return {"version":"v25.1",
       "status":"AUDIT_OF_LATEST_NOV17_2026_SOURCE_JSON_COUNTRY_BACKGROUND_DEPENDENCY",
       "original_source_Git_commit_before_this_analysis":SOURCE_COMMIT,
       "source_summary":summary,
       "country_heldout_original_reference_sensitivity":evidence,
       "Czech_same_52_species_comparison_of_road_vs_urban_domestic_references":cz,
       "two_domestic_reference_types_are_not_same_biological_control_habitat":True,
       "independent_2016_eastern_Slovak_19stations_LEDA_internal_null_source":{
         "n_genuine_author_stations":19,"source_true_318_station_species":318,
         "two_trait_incidence_coverage":sk["coverage"]["SK"]["source_incidence_joint_trait_coverage_fraction"],
         "CWM_distance_observed":m["original_observed"],
         "CWM_distance_expected_source_internal_null":m["internal_source_geography_fixed_margin_null_mean"],
         "observed_minus_source_internal_null":m["observed_minus_null"],
         "technical_null_interval_NOT_biological_CI":[
           m["technical_null_q025_NOT_field_CI"],m["technical_null_q975_NOT_field_CI"]],
         "not_a_fifth_railway_study_showing_same_positive_direction":True},
       "original_India_Kashmir_railway_Dryad_data_source_QA":{
          "yearly_repeat_is_same_railway_sites_not_new_independent_study":True,
          "actual_31_T_sites":k["source_spatial_sites"],
          "source_123_real_taxa":k["source_taxa"],
          "original_2D_LEDA_joined_taxa":k["both_traits"],
          "original_source_missing_cells_not_zero_imputed":k["missing_cells"]},
       "global_sPlotOpen_and_GeoNames_station_nearby_NOT_true_railway_plot":
            station_proximity(source["global_station_proximity.json"]),
       "scientific_limits":[
          "All 8 gates set after seeing results as a descriptive negative robustness audit, not confirmatory p-values",
          "Original study and country availability nonrandom; 10 country-control comparisons share papers, habitats and training floras",
          "Railway whole-station flora vs road mountain-site plot and urban disturbed/tree-pit small plots differ in area/years/effort",
          "Czech road versus city control on exactly same 52 taxa changes foreign-model advantage from strong to nearly zero",
          "Source genus resampling is species resampling, not a between-study population confidence interval",
          "Absolute MAE is affected by sample grain and non-detection, rather than a mechanistic propagule-flow probability",
          "Original near-station 250m sPlotOpen matching has zero country dataset with >=4 sites; 1km plots are not verified railway track or station land",
          "Original two nearby station buffers 1km and 2km select different country/dataset/plot cohorts; difference not physical distance dose response",
          "Station-wide and 100m ballast location flora, Slovakia rounded frequencies, India repeated years, AU/CN three corridors must not be pooled as identical samples",
          "Real station-level train origin-destination seed rain, survival, mowing/wind phenotypes and adaptive genetics remain unobserved"]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sources",default="research/validation_20261117")
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    data,qa=audited_source(a.sources)
    result=analyze(data)
    result["original_source_Git_blob_provenance"]=qa
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    summary={"version":result["version"],"status":result["status"],
        "source_summary":result["source_summary"],
        "Czech_same_52_species_comparison_of_road_vs_urban_domestic_references":
          result["Czech_same_52_species_comparison_of_road_vs_urban_domestic_references"],
        "independent_2016_eastern_Slovak_19stations_LEDA_internal_null_source":
          result["independent_2016_eastern_Slovak_19stations_LEDA_internal_null_source"],
        "global_sPlotOpen_and_GeoNames_station_nearby_NOT_true_railway_plot":
          result["global_sPlotOpen_and_GeoNames_station_nearby_NOT_true_railway_plot"]}
    print("V251_NOV17_COUNTRY_REFERENCE_EVIDENCE_AUDIT="+json.dumps(summary,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
