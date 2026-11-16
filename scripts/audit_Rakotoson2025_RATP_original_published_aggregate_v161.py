#!/usr/bin/env python3
"""v16.1 primary published AGGREGATE study audit for RATP railway
vegetation vs management alerts in Paris. No authored plot-matrix data.

Original: Rakotoson, Valantin-Morison, Meynard & Mediene (2025),
Weed Research 65(3) e70025, DOI 10.1111/wre.70025.
Published original full text:
https://onlinelibrary.wiley.com/doi/full/10.1111/wre.70025

The file data/literature/Rakotoson2025_original_published_aggregate_not_raw_v161.json
is a TRANSPARENt HUMAN-TRANSCRIBED MANIFEST of the ORIGINAL AUTHOR-PUBLISHED
summary, not a separately SHA-pinned publisher manuscript or plant plot data.
Only SOURCE-REPORTED means, counts, research design and limits are used.

Do not conflate 2021 plots with independent management events 2020–22.
Alerts (reports/requests) are NOT treatments; proximity <200m not true
individual plot mowing intensity, herbicide cessation in 2019 not a
longitudinal experimental baseline, and biodiversity indices alone do not
establish evolved adaptation, aerodynamic effects, seed carriage or fitness.

Reported source means lack ORIGINAL SD, plot-level species occurrence,
paired plots, spatial coordinates and local covariates. Compute only
author aggregate ratio/difference, NEVER Hedges g, confidence intervals,
meta-analysis weight or causal effect.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path

DOI="10.1111/wre.70025"
SOURCE_URL="https://onlinelibrary.wiley.com/doi/full/10.1111/wre.70025"
SOURCE_COUNTRY="France"
EXPECTED_N_PLOTS=394
EXPECTED_N_TRACK=128
EXPECTED_N_EMBANKMENT=266
EXPECTED_N_TAXA=255
EXPECTED_WEED_ALERTS=170
PLOT_AREA_M2=2
MANAGEMENT_RADIUS_M=200

# These metrics are SOURCE PUBLISHED summaries, not a site×taxon matrix.
SOURCE_MEAN_KEY_MAP={
 "mean_species_richness":(7.,2.4),
 "mean_original_abundance_sum_scale":(52.6,5.0),
 "mean_shannon":(1.55,.88),
 "mean_pielou":(.57,.61)
}

def assert_manifest(manifest):
    if not isinstance(manifest,dict) or manifest.get("version")!="v16.1_FR_RATP_author_original_published_AGGREGATES_only":
        raise ValueError("No valid authors' original published RATP summary manifest")
    if (manifest.get("original_author_DOI")!=DOI or
        manifest.get("original_author_publisher_full_article_url")!=SOURCE_URL):
        raise ValueError("Primary original publisher RATP DOI or source URL drift")
    source=manifest["field_author_samples_2021"]
    if ((source["embankment_n_surveys"],source["track_n_surveys"],source["total_n_surveys"],
         source["distinct_species_reported"],source["fixed_area_m2_per_vegetation_survey"])
        != (EXPECTED_N_EMBANKMENT,EXPECTED_N_TRACK,EXPECTED_N_PLOTS,EXPECTED_N_TAXA,PLOT_AREA_M2)):
        raise ValueError("Original Paris published rail vegetation source count/plot area drifted")
    if source["embankment_n_surveys"]+source["track_n_surveys"]!=source["total_n_surveys"]:
        raise ValueError("394 original vegetation sampling plots not conserved")
    if not all(source[k] is True for k in (
        "independent_pairing_IDs_NOT_published",
        "original_plot_by_taxon_presence_matrix_NOT_acquired",
        "original_plot_coordinates_NOT_acquired")):
        raise ValueError("Fake claim of original vegetation plot-level data acquired")
    management=manifest["railway_management_alerts"]
    if (management["source_weeding_alert_events"]!=EXPECTED_WEED_ALERTS or
        management["paper_weed_alert_proximity_threshold_m"]!=MANAGEMENT_RADIUS_M or
        management["source_time_window"]!="2020-2022"):
        raise ValueError("Original reported RATP management alert source changed")
    if (management["source_nearest_station_distance_half_of_alerts_below_m"]!=235 or
        not math.isclose(management["source_original_200m_radius_chosen_to_represent_fraction_of_alerts_approx"],.8)):
        raise ValueError("Original management risk/source proximity design changed")
    if not (management["year_2019_herbicide_policy_change_NO_before_after_plot_cohort"]
            and management["source_2021_plant_inventories_AND_2020_22_alerts_NOT_identical_grain"]):
        raise ValueError("Source 2019 glyphosate policy or independent event grain erased")
    if not(manifest["170_alerts_NOT_independent_mowing_treatments"]
           and manifest["394_surveys_are_NOT_394_distinct_railway_network_corridors"]
           and manifest["no_original_SD_or_plot_heterogeneity_for_Hedges_g"]):
        raise ValueError("Original source limits must be preserved to avoid false meta effect")
    if manifest["author_data_availability_statement"]!="Data will be made available on request.":
        raise ValueError("False claim that original publisher already released 394×255 species plot matrix")
    group=manifest["original_author_published_zone_average"]
    if set(group)!=set(SOURCE_MEAN_KEY_MAP):
        raise ValueError("Author-published aggregate metric set drifted")
    for key,(emb,track) in SOURCE_MEAN_KEY_MAP.items():
        metric=group[key]
        if not all(math.isclose(metric[k],v,rel_tol=0,abs_tol=1e-12)
                   for k,v in (("embankment",emb),("track",track))):
            raise ValueError("Published original RATP "+key+" changed")
    if group["mean_original_abundance_sum_scale"]["NOT_plant_biomass"] is not True:
        raise ValueError("Original railway abundance-index sum cannot become plant biomass")
    overlap=manifest["authors_original_descriptive_species_overlap"]
    if (not math.isclose(overlap["reported_fraction_with_at_least_one_shared_species"],.8) or
        not math.isclose(overlap["mean_copresence_shared_species_percent"],12.8) or
        not overlap["paired_plot_by_species_matrix_NOT_obtained"] or
        not overlap["species_overlap_NOT_direct_seed_flow"]):
        raise ValueError("Original source reporting co-presence is not proof of seed carriage")
    return True


def derive_author_published_aggregates(source):
    assert_manifest(source)
    orig=source["original_author_published_zone_average"]
    vals={}
    for name in sorted(SOURCE_MEAN_KEY_MAP):
        d=orig[name]
        a,b=d["embankment"],d["track"]
        vals[name]={
          "ORIGINAL_AUTHOR_PUBLISHED_EMBANKMENT_MEAN":a,
          "ORIGINAL_AUTHOR_PUBLISHED_TRACK_MEAN":b,
          "derived_simple_unstandardized_embankment_minus_track_mean":a-b,
          "derived_descriptive_embankment_over_track_mean_RATIO":a/b if b else None,
          "source_authors_reported_group_comparison":d["publisher_p_model_comparison"],
          "NOT_meta_analysis_effect_size_or_railway_mowing_treatment_effect":True,
        }
    return vals


def run(source,raw_bytes):
    aggregate=derive_author_published_aggregates(source)
    study=source["field_author_samples_2021"]
    result={
      "complete":True,"version":"v16.1_RATP2025_PRIMARY_AUTHOR_PUBLISHED_AGGREGATE_ONLY",
      "source_study_DOI":DOI,"primary_publisher_url":SOURCE_URL,
      "source_manifest_original_author_transcription_sha256":hashlib.sha256(raw_bytes).hexdigest(),
      "source_manifest_is_human_transcribed_NOT_downloaded_original_raw_PLOT_DATA":True,
      "original_source_plot_counts":{k:study[k] for k in (
          "total_n_surveys","track_n_surveys","embankment_n_surveys",
          "fixed_area_m2_per_vegetation_survey","distinct_species_reported")},
      "source_original_2020_22_reported_alert_events_NOT_2019_weed_mowing_experiments":EXPECTED_WEED_ALERTS,
      "published_200m_alert_PROXIMITY_classification_NOT_actual_mowing_exposure":True,
      "source_2019_herbicide_policy_NO_prepost_plant_comparison":True,
      "author_published_summary_metric_contrasts_NOT_study_meta_hedges_g":aggregate,
      "author_source_measured_80pct_plant_communities_sharing_species_NOT_seed_flux":source[
          "authors_original_descriptive_species_overlap"],
      "source_qualitative_functional_association_not_genetic_evolution":source[
          "authors_original_observational_functional_conclusion"],
      "independent_original_census_sites_observed_by_author_NOT_tetsu":True,
      "same_original_39station_Czech_or_89site_Slovenia_field_study":False,
      "true_394_plot_x_255_taxon_author_matrix_currently_obtained":False,
      "can_estimate_Hedges_g_or_valid_meta_variance_from_rounded_means":False,
      "can_estimate_mowing_causal_effect_with_170_alerts_not_treatment":False,
      "research_next_source_requirements":[
        "Original 394 plot×taxon identities/abundances, 2m2 plot coordinates, source line and paired observation IDs (currently on-request only)",
        "Original 170 dated geolocated weed-alert requests distinct from actual dated mowing, pulling, grubbing and herbicide treatments at plot scale",
        "Plot-level CWM annual/ruderal, environment, track/embankment, actual management history and observer detection effort",
        "Independent source road-embankment, nonrail meadow and region-climate comparator at same plot area before connecting the C disturbance axis",
        "Repeated same source plots before/after actual mechanical management to separate environment filtering from plasticity; common-garden/Fst for genetic adaptation",
      ],
      "interpretation":"Published source supports rail-track vs adjacent embankment habitat contrast and management-ALERT association; NOT causal effects of mowing/wind, train-borne seed arrival, genetic adaptation, or meta analysis."
    }
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    raw=Path(a.manifest).read_bytes()
    source=json.loads(raw)
    result=run(source,raw)
    dst=Path(a.output)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k in (
        "version","source_study_DOI","source_manifest_original_author_transcription_sha256",
        "original_source_plot_counts",
        "source_original_2020_22_reported_alert_events_NOT_2019_weed_mowing_experiments",
        "author_published_summary_metric_contrasts_NOT_study_meta_hedges_g",
        "true_394_plot_x_255_taxon_author_matrix_currently_obtained",
        "can_estimate_Hedges_g_or_valid_meta_variance_from_rounded_means",
        "can_estimate_mowing_causal_effect_with_170_alerts_not_treatment"
    )},ensure_ascii=False))

if __name__=="__main__":main()
