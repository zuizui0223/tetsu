"""v17.3 source sensitivity tests using SYNTHETIC biological site×species rows
only. No author field evidence is synthesized or presented as real results.
"""
import sys,math,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
import analyze_Slovenia_original_guild_vs_species_microhabitat_15section_predict_v172 as base
from test_Slovenia_species_vs_guild_section_holdout_v172 import toy_15_sections
from audit_Slovenia_species_guild_rare_common_15section_CV_v173 import (
  statistics,scenario_score,verify_original_v172_source,
  SOURCE_MAIN_REFERENCE,FREQ_DROP,MIN_ROUTE_SECTION_SUPPORT,
  SPECIES_SHRINKAGE_PSEUDOS,
)

class SlovenianOriginalTaxonInfluenceSourceQualityTests(unittest.TestCase):
    def test_original_author_source_baseline_v172_fixed_predictive_reference(self):
        self.assertAlmostEqual(
          SOURCE_MAIN_REFERENCE["all_observed_three_state"]["M2_original_source_species_identity_shrunk_to_guild"],
          1.0346230109451737)
        self.assertEqual(FREQ_DROP,(1,5,10))
        self.assertEqual(MIN_ROUTE_SECTION_SUPPORT,(3,5))
        self.assertEqual(SPECIES_SHRINKAGE_PSEUDOS,(3.,30.))

    def test_real_source_not_swapped_with_synthetic_15_section_material(self):
        with self.assertRaisesRegex(ValueError,"2729"):
            verify_original_v172_source(toy_15_sections())

    def test_frequency_selection_only_on_site_occurrences_NOT_species_BS_preference(self):
        x=toy_15_sections()
        frequency,sections,names=statistics(x)
        self.assertEqual(frequency["POANN"],30)
        self.assertEqual(frequency["BRERE"],30)
        self.assertEqual(len(sections["POANN"]),15)
        self.assertEqual(names["POANN"],"Poa annua")
        self.assertEqual(names["BRERE"],"Bromus erectus")
        self.assertEqual(sorted(frequency,key=lambda code:(-frequency[code],code))[0],"ALOTR")

    def test_conditional_tristate_and_exclusive_are_both_scored_after_removing_common_taxon(self):
        data=toy_15_sections()
        without=[r for r in data if r["source_taxon_EPPO"]!="POANN"]
        scenarios=scenario_score(without)
        self.assertTrue(scenarios["all_observed_three_state"]["identified"])
        self.assertTrue(scenarios["exclusive_BS_or_BT_only"]["identified"])
        self.assertEqual(scenarios["all_observed_three_state"]["n_original_route_sections_holdout"],15)
        self.assertEqual(scenarios["exclusive_BS_or_BT_only"]["n_original_route_sections_holdout"],15)
        self.assertEqual(scenarios["all_observed_three_state"]["n_original_trait_known_species"],3)
        self.assertEqual(scenarios["exclusive_BS_or_BT_only"]["n_original_trait_known_species"],1)
        self.assertTrue(all(math.isfinite(scenarios["all_observed_three_state"]["source_equal_section_logloss_M0_M1_M2"][m])
            for m in base.MODELS))

    def test_species_shrinkage_prior_changes_source_prediction_without_changing_original_guild_model(self):
        data=toy_15_sections()
        old=base.SPECIES_PSEUDO
        try:
            base.SPECIES_PSEUDO=3.
            x=scenario_score(data)["all_observed_three_state"]
            base.SPECIES_PSEUDO=30.
            y=scenario_score(data)["all_observed_three_state"]
        finally:
            base.SPECIES_PSEUDO=old
        gx=x["source_equal_section_logloss_M0_M1_M2"]
        gy=y["source_equal_section_logloss_M0_M1_M2"]
        self.assertAlmostEqual(gx[base.MODELS[0]],gy[base.MODELS[0]])
        self.assertAlmostEqual(gx[base.MODELS[1]],gy[base.MODELS[1]])
        self.assertNotAlmostEqual(gx[base.MODELS[2]],gy[base.MODELS[2]])
        self.assertEqual(base.SPECIES_PSEUDO,10.)

    def test_source_section_holdout_dependency_is_not_independent_pvalue(self):
        data=toy_15_sections()
        result=scenario_score(data)["all_observed_three_state"]
        self.assertTrue(result["no_real_genetic_or_plasticity_trait_measured"])
        self.assertTrue(result[
          "source_v172_same_source_section_score_fold_dropping_but_NOT_refitting"][
          "site_species_pairs_are_NOT_statistically_independent"])

if __name__=="__main__":unittest.main()
