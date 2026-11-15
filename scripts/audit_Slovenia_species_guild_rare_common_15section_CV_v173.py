#!/usr/bin/env python3
"""v17.3 SOURCE-AUTHENTIC Slovenia original 89 sites, 15 railway sections:
does the new original-source v17.2 cross-route species-identity predictive
gain depend on just the most frequent plants or the shrinkage parameter?

Original author publication: Šabić, Grašič & Rovanšek (2026),
DOI 10.1007/s44473-026-00155-z. Original publisher 401,279 byte XLSX
SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f.
Original source 89 paired BS/BT field sites, 15 original route sections,
261 joint-known source life-history+wind strategy taxa, 2729 observed
site×taxon presence opportunities (BS_only/BT_only/both).

This is a RE-TRAINED leave-one-ORIGINAL-section-out sensitivity:
  (i) same source data but remove top 1, top 5, top 10 most frequent
      original taxa, ranked by site incidence only (not by heldout outcome);
  (ii) only original taxa present in >=3 or >=5 author line sections;
  (iii) species-level shrinkage prior 3 and 30 vs frozen 10.
The species/source-site frequency filters are selected on the complete
original source data, including heldout sections: these are DESCRIPTIVE
selection sensitivities, not a second confirmatory holdout or external test.
M0 global vs M1 4-guild, M2 within-species ID-shrunk to M1. Original
source biological unit and 15-fold blocked structure never change.

Two separate original response definitions retain all source observations
(BS-only/BT-only/BOTH) and only exclusively located site×species events.
Phenotypic within-species response is NOT measured; species identity and
location-frequency persistence might reflect environmental microhabitat,
source pool, detection and nonuniform BS vs BT survey zone areas.
No direct measured train carriage, wind exposure, repetitive mowing,
genetic adaptation or independent matched RAIL-vs-nonrail flora.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path

import analyze_Slovenia_original_guild_vs_species_microhabitat_15section_predict_v172 as base
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx

SOURCE_MAIN_REFERENCE={
  "all_observed_three_state":{
       "M0_global":1.0982431044178191,
       "M1_source_lifespan_by_diaspore_guild":1.0855888865022834,
       "M2_original_source_species_identity_shrunk_to_guild":1.0346230109451737,
  },
  "exclusive_BS_or_BT_only":{
       "M0_global":0.69205610135632,
       "M1_source_lifespan_by_diaspore_guild":0.6777698225101966,
       "M2_original_source_species_identity_shrunk_to_guild":0.6482258245641337,
  },
}
FREQ_DROP=(1,5,10)
MIN_ROUTE_SECTION_SUPPORT=(3,5)
SPECIES_SHRINKAGE_PSEUDOS=(3.,30.)

def verify_original_v172_source(rows):
    if len(rows)!=2729 or len({r["source_taxon_EPPO"] for r in rows})!=261:
        raise ValueError("Authentic SI 2729 89field×261taxa incidence source study not recovered")
    if len({r["internal_original_site"] for r in rows})!=89 or len({r["author_original_section"] for r in rows})!=15:
        raise ValueError("Authentic original 15 independent source sections and 89 vegetation sites missing")
    for mode,expected in SOURCE_MAIN_REFERENCE.items():
        current=base.evaluate_study(rows,mode)
        actual=current["original_source_study_holdout_scores"][
          "out_of_original_source_section_log_loss"]["original_source_section_equal_site_equal_section"]
        for m,val in expected.items():
            if not math.isclose(actual[m],val,rel_tol=2e-10,abs_tol=1e-11):
                raise ValueError("Original source v17.2 15-section heldout true-SI flora prediction checksum drift: "+mode+" "+m)
    return True

def statistics(rows):
    by=Counter(r["source_taxon_EPPO"] for r in rows)
    sections=defaultdict(set)
    named={}
    for r in rows:
        taxon=r["source_taxon_EPPO"]
        sections[taxon].add(r["author_original_section"])
        named[taxon]=r["original_taxon_name_binomial_if_codebook_unambiguous"]
    return by,sections,named

def scenario_score(rows):
    result={}
    for mode in base.MODES:
        try:score=base.evaluate_study(rows,mode)
        except ValueError as err:
            result[mode]={"identified":False,"source_original_holdout_not_estimable_reason":str(err)}
            continue
        obs=score["original_source_study_holdout_scores"][
            "out_of_original_source_section_log_loss"][
            "original_source_section_equal_site_equal_section"]
        results=score["source_equal_section_improvement_percents_LOGLOSS"]
        result[mode]={
            "identified":True,
            "real_source_original_site_x_taxon_events":score["n_real_original_observed_source_site_x_species_model_events"],
            "n_original_railway_field_sites":score["n_original_unique_sites_source"],
            "n_original_route_sections_holdout":score["n_original_source_sections"],
            "n_original_trait_known_species":score["n_original_authors_distinct_species_original_observed"],
            "source_matched_microhabitat_state_events":
               score["n_original_BS_only_BT_only_BOTH_events"],
            "source_equal_section_logloss_M0_M1_M2":obs,
            "source_guild_and_species_vs_prior_MSE_not_produced_only_logloss_gain_pct":results,
            "n_source_sections_M1_guild_beats_M0":score["source_n_sections_M1_guild_beats_global_M0"],
            "n_source_sections_M2_identity_beats_M1":score["source_n_sections_M2_species_beats_guild_M1"],
            "source_v172_same_source_section_score_fold_dropping_but_NOT_refitting":
                {k:v for k,v in score["source_scored_fold_omission_only_M2_beyond_guild_influence"].items()
                 if k!="source_remaining_section_details"},
            "no_real_genetic_or_plasticity_trait_measured":True
        }
    return result

def study_sensitivity(original):
    rows,provenance=base.author_source_rows(original)
    verify_original_v172_source(rows)
    freq,sections,named=statistics(rows)
    ranked=sorted(freq,key=lambda code:(-freq[code],code))
    scenario={}
    scenario["original_FROZEN_species_pseudocount10_same_261_taxa"]=scenario_score(rows)
    for n in FREQ_DROP:
        remove=set(ranked[:n])
        retained=[r for r in rows if r["source_taxon_EPPO"] not in remove]
        scenario["retrain_after_REMOVING_"+str(n)+"_most_original_source_site_FREQUENT_species"]={
            "selection_USED_observed_site_frequency_not_habitat_effect":True,
            "source_removed_original_taxa":[{
                "original_EPPO":code,
                "author_source_codebook_binomial":named[code],
                "n_original_source_sites":freq[code],
                "n_original_sections":len(sections[code]),
            } for code in ranked[:n]],
            "two_original_response_definition_prediction":scenario_score(retained),
        }
    for minimum in MIN_ROUTE_SECTION_SUPPORT:
        retained=[r for r in rows if len(sections[r["source_taxon_EPPO"]])>=minimum]
        scenario["retrain_only_species_with_original_presence_in_at_least_"+str(minimum)+"_LINE_SECTIONS"]={
            "frequency_filter_CHOSEN_from_all_original_sections_not_a_separate_clean_holdout":True,
            "n_original_taxa_remaining":len({r["source_taxon_EPPO"] for r in retained}),
            "two_original_response_definition_prediction":scenario_score(retained)
        }
    preserved=base.SPECIES_PSEUDO
    if not math.isclose(preserved,10.):
        raise ValueError("Original frozen v17.2 species posterior prior pseudo-count !=10")
    for pseudo in SPECIES_SHRINKAGE_PSEUDOS:
        try:
            base.SPECIES_PSEUDO=pseudo
            scenario["retrain_full_source_with_species_identity_shrinkage_prior_"+str(int(pseudo))]={
                "original_main_v172_smoothed_prior_was10_NOT_tuned":True,
                "new_source_species_prior_pseudocount":pseudo,
                "two_original_response_definition_prediction":scenario_score(rows)
            }
        finally:
            base.SPECIES_PSEUDO=preserved
    return {
        "complete":True,
        "version":"v17.3_REAL_SI_original_BSBT_species_identity_cross_15route_common_taxon_leaveout",
        "source_original_source_study_DOI":"10.1007/s44473-026-00155-z",
        "source_original_publisher_XLSX_SHA":"26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f",
        "source_17_2_main_frozen_predicted_loss_reproduced":True,
        "source_original_89field_15section_261taxa_2729site_species_provenance":provenance,
        "n_original_source_study_taxa_removed_frequency_tests":list(FREQ_DROP),
        "n_source_rare_taxon_filt_original_route_section_tests":list(MIN_ROUTE_SECTION_SUPPORT),
        "original_source_model_unmodified_guild_prior_pseudocount10":True,
        "source_retrained_true_original_study_sensitivity_results":scenario,
        "no_3new_independent_original_field_studies":True,
        "limits":[
          "All sensitivity scenarios reuse ONE Slovenia original author railway BS vs BT study, 89 sites nested in 15 sections.",
          "Dropping species because they are globally frequent in the observed source uses heldout site presence data: exploratory influence diagnostic not a preregistered independent predictive test.",
          "This tests guild-based source microhabitat prediction vs species identity only and not actual selection mechanisms, railway vs road, seed movement or individual phenotypes.",
          "Site×species observations share train sections and individuals; no independent 2729 event p-values or meta effect.",
          "BS and BT surveyed area/cross-zone detection probabilities may vary; original authors' selected locations are not a randomized rail trial.",
          "Shrunk species-specific transfer across original rail sections is observational across plants and regions; source true evolutionary adaptation remains unidentifiable."
        ],
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-SI-original-XLSX",required=True)
    p.add_argument("--original-author-audit-v102",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.source_SI_original_XLSX,a.original_author_audit_v102)
    result=study_sensitivity(original)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    compact={}
    for k,v in result["source_retrained_true_original_study_sensitivity_results"].items():
        inp=v.get("two_original_response_definition_prediction",v)
        three=inp["all_observed_three_state"]
        excl=inp["exclusive_BS_or_BT_only"]
        compact[k]={
          "source_3_state_events":three.get("real_source_original_site_x_taxon_events"),
          "source_3_state_taxa":three.get("n_original_trait_known_species"),
          "3state_M1_guild_vs_M0_gain_pct":three.get("source_guild_and_species_vs_prior_MSE_not_produced_only_logloss_gain_pct",{}).get("guild_M1_vs_global_M0"),
          "3state_M2_species_vs_M1_gain_pct":three.get("source_guild_and_species_vs_prior_MSE_not_produced_only_logloss_gain_pct",{}).get("species_M2_vs_guild_M1"),
          "3state_n_source_sections_species_beats_guild":three.get("n_source_sections_M2_identity_beats_M1"),
          "exclusive_M2_species_vs_M1_gain_pct":excl.get("source_guild_and_species_vs_prior_MSE_not_produced_only_logloss_gain_pct",{}).get("species_M2_vs_guild_M1"),
          "exclusive_n_source_sections_species_beats_guild":excl.get("n_source_sections_M2_identity_beats_M1")
        }
    print(json.dumps({
         "version":result["version"],
         "main_original_source_v172_loss_reproduced":True,
         "source_common_species_omission_and_shrinkage":compact,
         "no_actual_plant_ITV_or_causal_mowing_wind_adaptation":True
    },ensure_ascii=False))

if __name__=="__main__":main()
