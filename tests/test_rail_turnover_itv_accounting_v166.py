"""v16.6 mathematics only; all arrays synthetic, NOT real rail-plants or photos."""
import sys,math,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_turnover_itv_accounting import decompose_rail_vs_control

class ExactCommunityTurnoverAndITVIdentity(unittest.TestCase):
    def test_only_composition_turnover_with_no_within_species_shift(self):
        rail={"A":(3,10.),"B":(1,2.)}
        control={"A":(1,10.),"B":(3,2.)}
        ref={"A":10.,"B":2.}
        z=decompose_rail_vs_control(rail,control,ref)
        self.assertAlmostEqual(z["observed_rail_minus_control_CWM"],4.)
        self.assertAlmostEqual(z["species_turnover_component_relative_to_reference"],4.)
        self.assertAlmostEqual(z["within_species_component_rail_minus_control"],0.)

    def test_only_within_species_shift_when_composition_unchanged(self):
        rail={"A":(1,8.),"B":(1,4.)}
        control={"A":(1,6.),"B":(1,4.)}
        r=decompose_rail_vs_control(rail,control,{"A":6.,"B":4.})
        self.assertAlmostEqual(r["observed_rail_minus_control_CWM"],1.)
        self.assertAlmostEqual(r["species_turnover_component_relative_to_reference"],0.)
        self.assertAlmostEqual(r["within_species_component_rail_minus_control"],1.)
        self.assertEqual(r["source_jointly_observed_species_trait_shifts"]["A"]["source_rail_minus_nonrail_within_species_trait_shift"],2.)

    def test_exact_reference_sensitive_decomposition_but_total_invariant(self):
        rail={"A":(3,9.),"B":(2,4.)}
        control={"A":(1,3.),"C":(2,7.)}
        a=decompose_rail_vs_control(rail,control,{"A":6.,"B":4.,"C":7.})
        b=decompose_rail_vs_control(rail,control,{"A":3.,"B":100.,"C":-3.})
        self.assertAlmostEqual(a["observed_rail_minus_control_CWM"],b["observed_rail_minus_control_CWM"])
        self.assertNotAlmostEqual(a["species_turnover_component_relative_to_reference"],
                                  b["species_turnover_component_relative_to_reference"])
        self.assertEqual(a["n_jointly_observed_with_real_within_species_comparison"],1)
        self.assertNotIn("B",a["source_jointly_observed_species_trait_shifts"])

    def test_missing_reference_in_single_habitat_species_cannot_be_imputed(self):
        with self.assertRaisesRegex(ValueError,"Cannot impute"):
            decompose_rail_vs_control({"A":(2,3.)},{"B":(2,5.)},{"A":3.})

    def test_invalid_nonfinite_and_absence_as_zero_declined(self):
        with self.assertRaisesRegex(ValueError,"positive finite"):
            decompose_rail_vs_control({"A":(0,4.)},{"A":(1,4.)},{"A":4.})
        with self.assertRaisesRegex(ValueError,"positive finite"):
            decompose_rail_vs_control({"A":(1,float("nan"))},{"A":(1,4.)},{"A":4.})
        with self.assertRaisesRegex(ValueError,"Independent species-level"):
            decompose_rail_vs_control({"A":(1,1.)},{"A":(1,1.)},None)

    def test_rail_CWM_does_not_measure_seed_delivery_or_genetic_adaptation(self):
        x=decompose_rail_vs_control({"Poa annua":(1,12.)},
                                    {"Poa annua":(1,10.)},{"Poa annua":10.})
        self.assertTrue(x["plasticity_vs_genetic_selection_NOT_identified_from_CWM"])
        self.assertTrue(x["no_rail_transport_or_photophenology_effect_inferred"])

if __name__=="__main__":unittest.main()
