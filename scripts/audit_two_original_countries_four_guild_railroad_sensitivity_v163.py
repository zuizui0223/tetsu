#!/usr/bin/env python3
"""v16.3 compare original v16.2 Czech vs Slovenian real-source functional
guild rail/road predictive robustness WITHOUT rehosting authors' raw studies.

Inputs are exactly the two separately successful AUTHENTIC ORIGINAL-SOURCE
v16.2 GitHub Actions JSON artifacts. Require source DOI & PBF SHA, same
sites/pairs and full 4-guild partition invariants from each publisher
source before comparing conditional prediction signatures.

Each class has TWO SAME original site-pair and guild-turnover outcomes with
whole-source vs guild-specific site richness adjustments. Main non-causal
stability diagnostic: rail-added heldout prediction improvement is POSITIVE
under BOTH adjustment choices. Leave-one-scored source heldout SECTION fold
omission is an influence diagnostic, not refitting or independent inference.

Do NOT pool the site Jaccard values, MSE improvements, 15 SI sections and
five Czech stretches into a global rail transport effect. Source SI strict
annual vs CZ Pladias therophyte are different original ontology concepts;
physical 2023/2024 OSM route distances are not train movement/actual seed flow.
This is an exploratory/falsification response to earlier results, not a
preregistered genuine multi-country railway seed-dispersal model.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from statistics import mean

V162="v16.2_TWO_ORIGINAL_COUNTRIES_DISTINCT_2x2_FUNCTIONAL_GUILD_RAILROAD"
STUDY_GATES={
 "SI":{"doi":"10.1007/s44473-026-00155-z",
       "pbf":"b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39",
       "sites":89,"segments":15,"pairs":3896,
       "guild_taxa":261,"guild_incidences":2729},
 "CZ":{"doi":"10.23855/preslia.2024.247",
       "pbf":"26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747",
       "sites":35,"segments":5,"pairs":594,
       "guild_taxa":533,"guild_incidences":4489}
}
GROUP_NAMES=("annual_wind","annual_nonwind","nonannual_wind","nonannual_nonwind")
RICHNESS_MODES=("guild_specific_log1p","whole_original_known_flora")
TOL=1.e-10

def evaluate_cv(cv,expected_segments):
    if not cv.get("identified"):
        return {"identified":False,"reason":cv.get("reason","No source subgroup common support"),
                "source_support":cv.get("source_support")}
    fold=cv["original_source_section_heldout_guild_fold_score"]
    n=cv["n_original_source_sections_eligible"]
    if n!=expected_segments or len(fold)!=n or len({r["original_author_stretch"] for r in fold})!=n:
        raise ValueError("Original source section-fold count or unique labels drifted")
    road=cv["source_original_section_holdout_MSE_M1_road"]
    rail=cv["source_original_section_holdout_MSE_M2_road_and_rail"]
    gain=cv["source_original_section_equal_weight_rail_after_road_MSE_gain_pct"]
    if not (math.isfinite(road) and math.isfinite(rail) and road>0 and rail>=0 and math.isfinite(gain)):
        raise ValueError("Original source crossvalidation real MSE not finite")
    if not math.isclose(100*(road-rail)/road,gain,abs_tol=TOL,rel_tol=1e-10):
        raise ValueError("Original first source original-fold MSE aggregate percentage inconsistent")
    values=[float(f["rail_increment_after_road_source_fold_MSE_reduction"]) for f in fold]
    if not math.isclose(mean(values),road-rail,abs_tol=TOL,rel_tol=1e-9):
        raise ValueError("Source observed original fold specific rail effect mean does not reproduce main MSE")
    if sum(x>0 for x in values)!=cv["n_original_source_sections_rail_addition_improved"]:
        raise ValueError("Source segment label and original rail gain signs inconsistent")
    # A fold deletion is an evaluation-only mean, NOT a new model training.
    omits=[{
        "omitted_original_section":f["original_author_stretch"],
        "n_remaining_source_folds":n-1,
        "rail_delta_MSE_average_after_eliminating_one_scored_fold":
            (sum(values)-values[i])/(n-1),
        "rail_addition_improves_remaining_source_scored_folds":
            (sum(values)-values[i])>0,
    } for i,f in enumerate(fold)]
    n_positive=sum(x["rail_addition_improves_remaining_source_scored_folds"] for x in omits)
    return {
       "identified":True,
       "source_original_heldout_folds":n,
       "original_road_MSE":road,
       "original_road_plus_rail_MSE":rail,
       "rail_increment_MSE_improvement_pct":gain,
       "rail_increment_positive":bool(gain>0),
       "source_sections_with_rail_addition_improvement":sum(x>0 for x in values),
       "n_scored_fold_deletions_retaining_positive_mean_rail_gain":n_positive,
       "n_scored_fold_deletions_reversing_original_gain_direction":sum(
           x["rail_addition_improves_remaining_source_scored_folds"]!=(gain>0)
           for x in omits),
       "original_source_fold_score_leaveout_not_refit":omits,
       "dependent_fold_errors_NO_independent_pvalues":True,
    }


def evaluate_study(study,d):
    fixed=STUDY_GATES[study]
    if d.get("version")!=V162 or d.get("source_original_country")!=study or not d.get("complete"):
        raise ValueError("Original v16.2 source study/version not recovered")
    ref=d["source_reference_assemblage_road_vs_rail_CV_reproduced"]
    if (ref["author_original_DOI"]!=fixed["doi"]
        or ref["original_field_source_2023_or_2024_OSM_SHA"]!=fixed["pbf"]
        or ref["n_original_study_sites"]!=fixed["sites"]
        or ref["n_original_line_sections"]!=fixed["segments"]
        or ref["n_original_comparable_rail_and_road_site_pairs"]!=fixed["pairs"]
        or not ref["no_source_assemblage_reference_change"]):
        raise ValueError("Original source study DOI, historic OSM SHA, site or rail/road pair frame changed")
    groups=d["original_source_fully_known_4guild_partition"]
    if set(groups)!=set(GROUP_NAMES):
        raise ValueError("Original source 2x2 four guilds lost or added")
    if (sum(x["source_distinct_observed_taxa"] for x in groups.values())!=fixed["guild_taxa"]
        or sum(x["source_guild_site_taxon_incidences"] for x in groups.values())!=fixed["guild_incidences"]):
        raise ValueError("Author original four-guild taxonomy/observed plant source incidences changed")
    results={}
    for name in GROUP_NAMES:
        grp=groups[name]
        variants=grp["richness_source_sensitivity"]
        if set(variants)!=set(RICHNESS_MODES):
            raise ValueError("Different author source group richness control methods than planned")
        saved={}
        for mode in RICHNESS_MODES:
            x=variants[mode]
            q=x["source_comparison_eligibility"]
            if q["source_full_original_comparable_site_pairs"]!=fixed["pairs"]:
                raise ValueError("Original source physical road-and-rail pair frame changed")
            if q["source_pairs_excluded_because_both_sites_have_zero_original_guild_species"]+q["source_guild_pair_comparisons_with_at_least_one_original_species"]!=fixed["pairs"]:
                raise ValueError("Original group double-zero similarity incorrectly imputed")
            a=evaluate_cv(x["original_source_heldout_rail_vs_road_prediction"],fixed["segments"])
            saved[mode]={
               "n_evaluable_source_pairs":q["source_guild_pair_comparisons_with_at_least_one_original_species"],
               "n_double_zero_group_pairs_excluded":q["source_pairs_excluded_because_both_sites_have_zero_original_guild_species"],
               **a,
            }
        if saved[RICHNESS_MODES[0]]["n_evaluable_source_pairs"]!=saved[RICHNESS_MODES[1]]["n_evaluable_source_pairs"]:
            raise ValueError("Original study paired 2-richness controls use different source sample frames")
        positive_both=all(saved[k]["identified"] and saved[k]["rail_increment_positive"]
                          for k in RICHNESS_MODES)
        positive_both_leaveout=positive_both and all(
            saved[k]["n_scored_fold_deletions_retaining_positive_mean_rail_gain"]
            ==fixed["segments"] for k in RICHNESS_MODES)
        results[name]={
            "source_group_taxa":grp["source_distinct_observed_taxa"],
            "source_site_taxon_incidence":grp["source_guild_site_taxon_incidences"],
            "source_sites_with_group_observed":grp["source_sites_with_at_least_one_group_taxon"],
            "n_original_source_sites_total":fixed["sites"],
            "two_original_richness_control_modes":saved,
            "positive_under_BOTH_richness_controls":positive_both,
            "positive_under_BOTH_and_every_source_fold_deletion":positive_both_leaveout,
            "not_original_transport_effect":True,
        }
    return results


def run(si,cz):
    studies={"SI":evaluate_study("SI",si),"CZ":evaluate_study("CZ",cz)}
    summaries={}
    for group in GROUP_NAMES:
        s=studies["SI"][group]
        c=studies["CZ"][group]
        summaries[group]={
          "Slovenia_original_reported_annual_IS_strict_lifespan_Czech_IS_therophyte_proxy":True,
          "n_truly_independent_original_field_studies":2,
          "n_original_countries":2,
          "n_country_sources_positive_under_both_richness_controls":sum(
              studies[k][group]["positive_under_BOTH_richness_controls"] for k in ("SI","CZ")),
          "both_original_countries_have_stable_positive_rail_increment":
              s["positive_under_BOTH_richness_controls"] and c["positive_under_BOTH_richness_controls"],
          "both_original_countries_have_stable_negative_rail_increment_for_both_richness_modes":
              all(studies[country][group]["two_original_richness_control_modes"][mode]["identified"]
                   and not studies[country][group]["two_original_richness_control_modes"][mode]["rail_increment_positive"]
                   for country in ("SI","CZ") for mode in RICHNESS_MODES),
          "same_group_not_equivalent_source_trait_definition_across_countries":True,
        }
    all_evaluated=[(ct,g,m,studies[ct][g]["two_original_richness_control_modes"][m])
                   for ct in ("SI","CZ") for g in GROUP_NAMES for m in RICHNESS_MODES]
    all_primary=[x for x in all_evaluated if x[2]==RICHNESS_MODES[0]]
    both_count=sum(studies[c][g]["positive_under_BOTH_richness_controls"]
                   for c in ("SI","CZ") for g in GROUP_NAMES)
    return {
      "complete":True,
      "version":"v16.3_two_source_real_4guild_road_rail_MSE_control_sensitivity",
      "original_independent_studies":["10.1007/s44473-026-00155-z","10.23855/preslia.2024.247"],
      "real_source_2x2_country_guild_details":studies,
      "two_original_country_qualitative_pattern_by_functional_guild":summaries,
      "full_original_source_country_guild_richness_evaluated":len(all_evaluated),
      "n_original_study_guild_pairs":2*len(GROUP_NAMES),
      "n_country_study_guild_combinations_positive_for_primary_specific_richness":sum(
          x[3]["identified"] and x[3]["rail_increment_positive"] for x in all_primary),
      "n_country_study_guild_combinations_positive_under_BOTH_richness_controls":both_count,
      "n_country_study_guild_combinations_positive_under_BOTH_and_all_fold_deletions":sum(
          studies[c][g]["positive_under_BOTH_and_every_source_fold_deletion"]
          for c in ("SI","CZ") for g in GROUP_NAMES),
      "n_group_types_positive_in_both_independent_original_studies":sum(
          x["both_original_countries_have_stable_positive_rail_increment"] for x in summaries.values()),
      "n_group_types_negative_across_both_studies_and_both_richness_controls":sum(
          x["both_original_countries_have_stable_negative_rail_increment_for_both_richness_modes"]
          for x in summaries.values()),
      "important_limits":[
        "This is posthoc 4 guilds × 2 original botanical studies, not 8 or 16 independent studies or randomized train experiments.",
        "SI author strict annual lifespan is NOT the same field trait as Czech Raunkiaer therophyte proxy; class names are approximate and no pooling of coefficients is statistically justified.",
        "Original wind-diaspore strategy labels are not measured train slipstream, on-train propogules, survival, fitness or evolved mowing resistance.",
        "SI 2022-23 original field transect GPS differs from Czech 2020-22 whole-station vegetation with independent 2026 registry CENTROIDS.",
        "Both OSM 2023/2024 are historical map shortest paths; train traffic OD, species source pools, actual mowing episodes and detection effort remain unknown.",
        "Source fold-deletion analysis removes an already-scored validation fold, NOT model retraining; shared original plant sites make folds/pairs non-independent.",
        "No genuine per-species train seed catch->site recruitment->offspring or genetic adaptation can be estimated from these original site-flora matrices.",
        "Positive MSE improvement even if source-consistent would be observational predictive association, not actual train carried-seed causality.",
      ],
      "NO_pooled_global_railway_dispersal_causal_effect":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-SI-v162",required=True)
    p.add_argument("--original-CZ-v162",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    x=run(json.loads(Path(a.original_SI_v162).read_text(encoding="utf-8")),
          json.loads(Path(a.original_CZ_v162).read_text(encoding="utf-8")))
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
       "version":x["version"],
       "primary_subgroup_cv_n_positive":x["n_country_study_guild_combinations_positive_for_primary_specific_richness"],
       "robust_two_richness_positive_n":x["n_country_study_guild_combinations_positive_under_BOTH_richness_controls"],
       "robust_two_richness_plus_original_fold_deletion_positive_n":x["n_country_study_guild_combinations_positive_under_BOTH_and_all_fold_deletions"],
       "cross_original_study_both_positive_guild_types":x["n_group_types_positive_in_both_independent_original_studies"],
       "cross_original_study_negative_both_controls_guild_types":x["n_group_types_negative_across_both_studies_and_both_richness_controls"],
       "actual_source_country_guild_details":{
          st:{
            g:{
              "taxa":d["source_group_taxa"],
              "source_incidences":d["source_site_taxon_incidence"],
              "site_present":d["source_sites_with_group_observed"],
              "robust_positive":d["positive_under_BOTH_richness_controls"],
              "modes":{
                mode:{
                  "CV_positive_rail_gain_pct":y.get("rail_increment_MSE_improvement_pct"),
                  "fold_score_dropone_positive_count":y.get("n_scored_fold_deletions_retaining_positive_mean_rail_gain"),
                  "n_double_zero_source_pair_excluded":y["n_double_zero_group_pairs_excluded"],
                } for mode,y in d["two_original_richness_control_modes"].items()
              },
            } for g,d in val.items()
          } for st,val in x["real_source_2x2_country_guild_details"].items()
       },
       "no_source_meta_transport_effect":True
    },ensure_ascii=False))


if __name__=="__main__":main()
