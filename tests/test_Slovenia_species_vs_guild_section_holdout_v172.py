"""v17.2 synthetic mathematical tests ONLY.
No test-constructed taxa, site, track section or BS/BT placement below
are real botanist-observed railway plants or adaptation measurements.
"""
import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_original_guild_vs_species_microhabitat_15section_predict_v172 import (
   STATES,MODELS,MODES,source_smoothed_model,predict,scores,
   mode_data,score_one_holdout,evaluate_study,interval_from_fold_omission
)

def r(site,section,taxon,guild,state):
    return {"internal_original_site":str(site),
            "author_original_section":str(section),
            "source_taxon_EPPO":taxon,
            "original_taxon_name_binomial_if_codebook_unambiguous":
               "Poa annua" if taxon=="POANN" else "Bromus erectus",
            "source_four_lifespan_diaspore_guild":guild,
            "observed_original_microhabitat_state":state}

def toy_15_sections():
    g="strict_annual_wind"
    h="other_known_other"
    obs=[]
    for sec in range(15):
        # Identical taxon known in each line section, alternating
        # within-site two observed microhabitats in this synthetic fixture.
        for sub in range(2):
            site=f"{sec}_{sub}"
            obs.append(r(site,sec,"POANN",g,"BS_only"))
            obs.append(r(site,sec,"BRERE",h,"BT_only"))
            obs.append(r(site,sec,"ALOTR",g,"both"))
            if sub:
                obs.append(r(site,sec,"CRANU",h,"both"))
    return obs

class OriginalSourceSpeciesVsGuildPredictionTests(unittest.TestCase):
    def test_train_probs_are_shrunk_genuine_three_source_state_simplex(self):
        rows=[
          r(1,1,"POANN","strict_annual_wind","BS_only"),
          r(2,1,"POANN","strict_annual_wind","BS_only"),
          r(3,1,"BRERE","other_known_other","BT_only"),
          r(4,1,"OTHER","other_known_other","both"),
        ]
        fit=source_smoothed_model(rows,STATES)
        p=predict(r(5,2,"POANN","strict_annual_wind","BS_only"),fit)
        for name in MODELS:
            self.assertAlmostEqual(sum(p[name].values()),1.)
            self.assertGreater(min(p[name].values()),0.)
        self.assertTrue(p["original_heldout_species_has_prior_original_train_observations"])
        self.assertGreater(p[MODELS[2]]["BS_only"],p[MODELS[1]]["BS_only"])
        new=predict(r(5,2,"NEW","strict_annual_wind","BT_only"),fit)
        self.assertFalse(new["original_heldout_species_has_prior_original_train_observations"])
        self.assertEqual(new[MODELS[2]],new[MODELS[1]])
        with self.assertRaisesRegex(ValueError,"trait code differs"):
            predict(r(5,2,"POANN","other_known_other","BT_only"),fit)

    def test_binary_exclusive_does_not_invent_absent_or_reclassify_both(self):
        rows=toy_15_sections()
        three,labels=mode_data(rows,MODES[0])
        binary,labels2=mode_data(rows,MODES[1])
        self.assertEqual(len(labels),3)
        self.assertEqual(len(labels2),2)
        self.assertLess(len(binary),len(three))
        self.assertNotIn("both",{r["observed_original_microhabitat_state"] for r in binary})
        self.assertEqual(
            sum(r["observed_original_microhabitat_state"]=="both" for r in three),
            len(three)-len(binary))
        with self.assertRaisesRegex(ValueError,"Unrecognized"):
            mode_data(rows,"invented_railway_seed_carriage")

    def test_original_15_source_section_holdout_both_endpoints(self):
        dataset=toy_15_sections()
        three=evaluate_study(dataset,MODES[0])
        binary=evaluate_study(dataset,MODES[1])
        self.assertEqual(three["n_original_source_sections"],15)
        self.assertEqual(binary["n_original_source_sections"],15)
        self.assertEqual(three["n_original_unique_sites_source"],30)
        self.assertEqual(three["n_original_authors_distinct_species_original_observed"],4)
        self.assertEqual(binary["n_original_authors_distinct_species_original_observed"],2)
        self.assertEqual(sum(three["n_original_BS_only_BT_only_BOTH_events"].values()),len(dataset))
        self.assertEqual(sum(binary["n_original_BS_only_BT_only_BOTH_events"].values()),
                         binary["n_real_original_observed_source_site_x_species_model_events"])
        for mode in (three,binary):
            for scoring in mode["original_source_study_holdout_scores"].values():
                self.assertEqual(set(scoring["original_source_section_equal_site_equal_section"]),set(MODELS))
                self.assertTrue(all(math.isfinite(v) for v in
                     scoring["original_source_section_equal_site_equal_section"].values()))
            self.assertEqual(len(mode["original_heldout_section_scores"]),15)
            self.assertTrue(mode["not_a_measurement_of_within_species_genetic_or_plastic_adaptation"])

    def test_exact_log_loss_and_multiclass_Brier(self):
        p={"BS_only":.65,"BT_only":.25,"both":.10}
        ll,b=scores("BS_only",p,STATES)
        self.assertAlmostEqual(ll,-math.log(.65))
        self.assertAlmostEqual(b,((1-.65)**2+.25**2+.1**2)/3)
        with self.assertRaisesRegex(ValueError,"absent"):
            scores("actual_gene_flow",p,STATES)

    def test_folds_source_equal_mean_vs_event_weighting_distinguished(self):
        data=toy_15_sections()
        result=evaluate_study(data,MODES[0])
        all_folds=result["original_heldout_section_scores"]
        vals=[x["site_mean_equal_weight_scores"]["out_of_original_source_section_log_loss"][MODELS[2]]
              for x in all_folds]
        self.assertAlmostEqual(mean_loc:=sum(vals)/len(vals),
               result["original_source_study_holdout_scores"][
               "out_of_original_source_section_log_loss"][
               "original_source_section_equal_site_equal_section"][MODELS[2]])
        self.assertGreater(mean_loc,0)
        self.assertEqual(
            result["source_scored_fold_omission_only_M2_beyond_guild_influence"][
               "source_remaining_section_details"][0]["n_source_sites"]
               if "n_source_sites" in result["source_scored_fold_omission_only_M2_beyond_guild_influence"][
                 "source_remaining_section_details"][0] else 15,15)

    def test_original_all_15_sections_required(self):
        with self.assertRaisesRegex(ValueError,"15 line sections"):
            evaluate_study([r(1,1,"POANN","strict_annual_wind","BS_only")],MODES[0])

    def test_no_false_species_study_independence_significance(self):
        result=evaluate_study(toy_15_sections(),MODES[0])
        impact=result["source_scored_fold_omission_only_M2_beyond_guild_influence"]
        self.assertTrue(impact["original_section_scores_excluded_without_MODEL_RETRAINING"])
        self.assertTrue(impact["site_species_pairs_are_NOT_statistically_independent"])
        self.assertLessEqual(impact[
           "n_original_scored_fold_omissions_maintaining_positive_species_added_gain"],15)
        self.assertTrue(all(
            r["only_an_EXPLORATORY_species_priority_score_NOT_adaptation"] for r in
            result["source_ranked_high_frequency_species_outofsample_descriptive_NOT_global_rail_preference"]
        ))

if __name__=="__main__":unittest.main()
