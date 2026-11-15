"""v16.9 SYNTHETIC pure maths tests: paired observed within-species
trait contrast vs unpaired-only species is NOT the same estimand.
No real track plants, iNaturalist measurements, train seed transport or
heritable phenotype data are included in synthetic arrays.
"""
import sys,math,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_turnover_itv_accounting import decompose_rail_vs_control

class ExactPairedSpeciesItvAndUnpairedTurnoverTests(unittest.TestCase):
    def test_no_shared_species_must_have_ZERO_direct_observed_ITV(self):
        x=decompose_rail_vs_control(
            {"Poa annua":(1,7.)},{"Bromus erectus":(1,2.)},
            {"Poa annua":5.,"Bromus erectus":3.})
        self.assertAlmostEqual(x["observed_rail_minus_control_CWM"],5.)
        self.assertAlmostEqual(x["species_turnover_component_relative_to_reference"],2.)
        self.assertAlmostEqual(x["STRICT_paired_only_shared_species_ITV_component"],0.)
        self.assertAlmostEqual(x["shared_species_reference_composition_by_trait_interaction"],0.)
        self.assertAlmostEqual(x["species_only_one_habitat_reference_deviation_NOT_observed_ITV"],3.)
        self.assertEqual(x["n_jointly_observed_with_real_within_species_comparison"],0)
        self.assertEqual(x["n_species_only_rail"],1)
        self.assertEqual(x["n_species_only_control"],1)
        self.assertTrue(x["old_within_species_component_INCLUDES_reference_deviations_of_species_found_in_ONLY_ONE_habitat_NOT_pure_paired_ITV"])

    def test_only_true_paired_species_shift_when_weights_identical(self):
        x=decompose_rail_vs_control(
            {"A":(1,8.),"B":(1,4.)},{"A":(1,6.),"B":(1,4.)},
            {"A":6.,"B":4.})
        self.assertAlmostEqual(x["STRICT_paired_only_shared_species_ITV_component"],1.)
        self.assertAlmostEqual(x["shared_species_reference_composition_by_trait_interaction"],0.)
        self.assertAlmostEqual(x["species_only_one_habitat_reference_deviation_NOT_observed_ITV"],0.)
        self.assertAlmostEqual(x["strict_four_term_reconstruction"],1.)

    def test_all_shared_species_composition_reference_interaction_can_be_nonzero(self):
        x=decompose_rail_vs_control(
          {"A":(3,10.),"B":(1,2.)},
          {"A":(1,6.),"B":(3,4.)},
          {"A":5.,"B":5.})
        self.assertAlmostEqual(x["observed_rail_minus_control_CWM"],3.5)
        self.assertAlmostEqual(x["STRICT_paired_only_shared_species_ITV_component"],1.)
        self.assertAlmostEqual(x["shared_species_reference_composition_by_trait_interaction"],2.5)
        self.assertAlmostEqual(x["species_only_one_habitat_reference_deviation_NOT_observed_ITV"],0.)

    def test_nonmatched_reference_residual_does_not_invent_control_trait(self):
        x=decompose_rail_vs_control(
          {"A":(3,9.),"B":(2,4.)},
          {"A":(1,3.),"C":(2,7.)},
          {"A":6.,"B":5.,"C":9.})
        self.assertAlmostEqual(x["STRICT_paired_only_shared_species_ITV_component"],
                               ((3/5+1/3)/2)*(9-3))
        self.assertAlmostEqual(x["species_only_one_habitat_reference_deviation_NOT_observed_ITV"],
                               (2/5)*(4-5)-(2/3)*(7-9))
        self.assertAlmostEqual(x["strict_four_term_reconstruction"],
                               x["observed_rail_minus_control_CWM"])
        self.assertEqual(set(x["source_jointly_observed_species_trait_shifts"]),{"A"})
        self.assertTrue(x["no_trait_mean_for_absent_species_imputed"])

    def test_change_external_trait_reference_shifts_components_not_exact_total(self):
        rail={"A":(3,9.),"B":(2,4.)}
        control={"A":(1,3.),"C":(2,7.)}
        a=decompose_rail_vs_control(rail,control,{"A":6.,"B":4.,"C":7.})
        b=decompose_rail_vs_control(rail,control,{"A":0.,"B":100.,"C":-3.})
        self.assertAlmostEqual(a["observed_rail_minus_control_CWM"],
                               b["observed_rail_minus_control_CWM"])
        self.assertNotAlmostEqual(a["species_turnover_component_relative_to_reference"],
                                  b["species_turnover_component_relative_to_reference"])
        self.assertAlmostEqual(a["STRICT_paired_only_shared_species_ITV_component"],
                               b["STRICT_paired_only_shared_species_ITV_component"])
        self.assertAlmostEqual(a["strict_four_term_reconstruction"],
                               b["strict_four_term_reconstruction"])

    def test_missing_original_taxon_reference_refuses_artificial_ITV(self):
        with self.assertRaisesRegex(ValueError,"Cannot impute"):
            decompose_rail_vs_control(
               {"A":(1,4.)},{"B":(1,8.)},{"A":4.})

if __name__=="__main__":unittest.main()
