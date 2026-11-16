"""v16.3 synthetic model-score and original-source metadata consistency tests.
No synthetic row is an observation of botanical field traits or train seeds.
"""
import copy,math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_two_original_countries_four_guild_railroad_sensitivity_v163 import (
    STUDY_GATES,GROUP_NAMES,RICHNESS_MODES,V162,run,evaluate_cv,evaluate_study
)

def source_fake(study):
    d=STUDY_GATES[study]
    labels=["annual_wind","annual_nonwind","nonannual_wind","nonannual_nonwind"]
    original_taxa=[10,51,43,157] if study=="SI" else [19,157,70,287]
    original_incidences=[180,471,438,1640] if study=="SI" else [253,1332,703,2201]
    groups={}
    for idx,name in enumerate(labels):
        values={}
        for mode in RICHNESS_MODES:
            # synthetic positivity switch, explicitly NOT actual original results
            is_pos=(idx==0 and mode=="guild_specific_log1p")
            rail=0.9 if is_pos else 1.1
            fold=[{"original_author_stretch":str(i),
                   "source_n_test_dependent_site_pairs":100,
                   "rail_increment_after_road_source_fold_MSE_reduction":1-rail}
                  for i in range(d["segments"])]
            values[mode]={
                "source_comparison_eligibility":{
                    "source_full_original_comparable_site_pairs":d["pairs"],
                    "source_pairs_excluded_because_both_sites_have_zero_original_guild_species":0,
                    "source_guild_pair_comparisons_with_at_least_one_original_species":d["pairs"],
                },
                "original_source_heldout_rail_vs_road_prediction":{
                    "identified":True,
                    "source_original_section_holdout_MSE_M1_road":1.,
                    "source_original_section_holdout_MSE_M2_road_and_rail":rail,
                    "source_original_section_equal_weight_rail_after_road_MSE_gain_pct":(1-rail)*100,
                    "n_original_source_sections_rail_addition_improved":d["segments"] if is_pos else 0,
                    "n_original_source_sections_eligible":d["segments"],
                    "original_source_section_heldout_guild_fold_score":fold,
                }
            }
        groups[name]={
            "source_distinct_observed_taxa":original_taxa[idx],
            "source_guild_site_taxon_incidences":original_incidences[idx],
            "source_sites_with_at_least_one_group_taxon":d["sites"],
            "richness_source_sensitivity":values,
        }
    return {
        "complete":True,"version":V162,"source_original_country":study,
        "source_reference_assemblage_road_vs_rail_CV_reproduced":{
           "author_original_DOI":d["doi"],
           "original_field_source_2023_or_2024_OSM_SHA":d["pbf"],
           "n_original_study_sites":d["sites"],
           "n_original_line_sections":d["segments"],
           "n_original_comparable_rail_and_road_site_pairs":d["pairs"],
           "no_source_assemblage_reference_change":True,
        },
        "original_source_fully_known_4guild_partition":groups,
    }


class TwoSource4GuildAuditMathTests(unittest.TestCase):
    def test_positive_primary_but_negative_full_richness_not_robust(self):
        o=run(source_fake("SI"),source_fake("CZ"))
        self.assertEqual(o["n_original_study_guild_pairs"],8)
        self.assertEqual(o["full_original_source_country_guild_richness_evaluated"],16)
        self.assertEqual(o["n_country_study_guild_combinations_positive_for_primary_specific_richness"],2)
        self.assertEqual(o["n_country_study_guild_combinations_positive_under_BOTH_richness_controls"],0)
        self.assertEqual(o["n_group_types_positive_in_both_independent_original_studies"],0)
        self.assertTrue(o["NO_pooled_global_railway_dispersal_causal_effect"])

    def test_all_source_artifact_doi_and_geofabrik_sha_validation(self):
        c=source_fake("CZ")
        c["source_reference_assemblage_road_vs_rail_CV_reproduced"][
           "original_field_source_2023_or_2024_OSM_SHA"]="fabricated"
        with self.assertRaisesRegex(ValueError,"DOI, historic OSM SHA"):
            evaluate_study("CZ",c)

    def test_source_four_guild_partition_taxa_and_incidence_conservation(self):
        si=source_fake("SI")
        si["original_source_fully_known_4guild_partition"]["annual_wind"]["source_distinct_observed_taxa"]+=1
        with self.assertRaisesRegex(ValueError,"taxonomy/observed plant"):
            evaluate_study("SI",si)

    def test_original_double_zero_group_pairs_cannot_be_imputed_as_zero(self):
        si=source_fake("SI")
        x=si["original_source_fully_known_4guild_partition"]["annual_wind"]["richness_source_sensitivity"]["guild_specific_log1p"]["source_comparison_eligibility"]
        x["source_pairs_excluded_because_both_sites_have_zero_original_guild_species"]=1
        with self.assertRaisesRegex(ValueError,"double-zero"):
            evaluate_study("SI",si)

    def test_score_fold_mean_must_equal_aggregate_MSE_gain(self):
        fake=source_fake("CZ")
        cv=fake["original_source_fully_known_4guild_partition"]["annual_wind"]["richness_source_sensitivity"]["guild_specific_log1p"]["original_source_heldout_rail_vs_road_prediction"]
        cv["source_original_section_holdout_MSE_M2_road_and_rail"]=.8
        with self.assertRaisesRegex(ValueError,"inconsistent"):
            evaluate_cv(cv,5)

    def test_fold_drop_sign_detects_single_large_outlier(self):
        cv={
           "identified":True,
           "source_original_section_holdout_MSE_M1_road":1.,
           "source_original_section_holdout_MSE_M2_road_and_rail":1.04,
           "source_original_section_equal_weight_rail_after_road_MSE_gain_pct":-4.,
           "n_original_source_sections_rail_addition_improved":4,
           "n_original_source_sections_eligible":5,
           "original_source_section_heldout_guild_fold_score":[
               {"original_author_stretch":str(i),
                "source_n_test_dependent_site_pairs":5,
                "rail_increment_after_road_source_fold_MSE_reduction":-1. if i==0 else .2}
               for i in range(5)
           ],
        }
        d=evaluate_cv(cv,5)
        self.assertFalse(d["rail_increment_positive"])
        self.assertEqual(d["n_scored_fold_deletions_retaining_positive_mean_rail_gain"],1)
        self.assertEqual(d["n_scored_fold_deletions_reversing_original_gain_direction"],1)
        self.assertTrue(d["dependent_fold_errors_NO_independent_pvalues"])

if __name__=="__main__":unittest.main()
