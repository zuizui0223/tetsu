"""Numerical and estimand guards for the v315 thermal comparison."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
    import analyze_thermal_affinity_reference_v315 as m
except ModuleNotFoundError as error:
    if error.name!="analyze_thermal_affinity_reference_v315": raise
    m=None


def trait(values,states=None,limits=(1.,12.)):
    return {"values":values,"statuses":states or {key:"numeric" for key in values},
            "limits":list(limits),"audit":{}}


class ThermalReferenceTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(m,"v315 numerical implementation not present")

    def test_indifference_and_missing_do_not_enter_numeric_mean_as_zero(self):
        t=trait({"A":2.,"B":8.},{"A":"numeric","B":"numeric","X":"indifferent","U":"missing"})
        out=m.summarize({"A":.6,"B":.2,"X":.1,"U":.1},t,min_species=2)
        self.assertAlmostEqual(out["mean"],3.5)
        self.assertAlmostEqual(out["numeric_frequency_coverage"],.8)
        self.assertAlmostEqual(out["status_frequency_fraction"]["indifferent"],.1)
        self.assertAlmostEqual(out["indifferent_fraction_among_known_status"],1/9)

    def test_species_weighting_and_frequency_rescaling_have_distinct_effects(self):
        t=trait({"A":2.,"B":8.})
        a=m.summarize({"A":.6,"B":.2},t,min_species=2)
        b=m.summarize({"A":.3,"B":.1},t,min_species=2)
        c=m.summarize({"A":.6,"B":.2},t,min_species=2,species_equal=True)
        self.assertAlmostEqual(a["mean"],3.5)
        self.assertAlmostEqual(a["mean"],b["mean"])
        self.assertAlmostEqual(c["mean"],5.)

    def test_same_species_has_same_value_and_habitat_exchange_flips_difference(self):
        t=trait({"A":2.,"B":8.})
        self.assertEqual(m.contrast({"A":1.},{"A":.5},t,min_species=1)["difference"],0.)
        a=m.contrast({"A":1.},{"B":1.},t,min_species=1)
        b=m.contrast({"B":1.},{"A":1.},t,min_species=1)
        self.assertEqual(a["difference"],-6.)
        self.assertEqual(b["difference"],6.)

    def test_taxonomic_unknown_frequency_masks_both_habitats(self):
        t=trait({"A":2.,"B":8.,"C":6.})
        out=m.contrast({"A":1.,"B":None},{"A":1.,"B":1.,"C":1.},t,
                       exclude={"C"},min_species=1)
        self.assertEqual(out["rail"]["n_named"],1)
        self.assertEqual(out["reference"]["n_named"],1)
        self.assertEqual(out["difference"],0.)

    def test_joint_unknown_envelope_cancels_a_shared_species(self):
        out=m.shared_unknown_envelope({"A":.5,"U":.5},{"A":1.,"U":1.},{"A":5.},(0.,10.))
        self.assertAlmostEqual(out["lower"],0.)
        self.assertAlmostEqual(out["upper"],0.)
        self.assertAlmostEqual(out["unknown_coefficient_l1"],0.)

    def test_joint_unknown_envelope_has_hand_calculated_endpoints(self):
        out=m.shared_unknown_envelope({"A":.8,"U":.2},{"A":.6,"U":.4},{"A":5.},(0.,10.))
        self.assertAlmostEqual(out["known_component"],1.)
        self.assertAlmostEqual(out["lower"],-1.)
        self.assertAlmostEqual(out["upper"],1.)
        self.assertAlmostEqual(out["unknown_coefficient_l1"],.2)

    def test_zero_is_valid_in_eive_but_out_of_scale_and_low_support_stop(self):
        t=trait({"A":0.,"B":10.},limits=(0.,10.))
        self.assertEqual(m.summarize({"A":1.,"B":1.},t,min_species=2)["mean"],5.)
        self.assertIsNone(m.summarize({"A":1.,"B":1.},t)["mean"])
        with self.assertRaises(ValueError):
            m.summarize({"A":1.},trait({"A":12.},limits=(0.,10.)),min_species=1)

    def test_common_numeric_universe_does_not_require_presence_in_both_habitats(self):
        a=trait({"RailOnly":8.,"ReferenceOnly":2.,"Extra":6.})
        b=trait({"RailOnly":7.,"ReferenceOnly":3.},limits=(0.,10.))
        selected=m.common_numeric_traits(a,b)
        out=m.contrast({"RailOnly":1.,"Extra":1.},{"ReferenceOnly":1.},selected[0],min_species=1)
        self.assertEqual(out["difference"],6.)
        self.assertEqual(out["rail"]["n_numeric"],1)
        self.assertEqual(out["rail"]["n_named"],2)

    def test_paired_summary_weights_stations_equally_and_reports_exclusions(self):
        t=trait({"A":2.,"B":8.})
        pairs=[({"B"},{"A"}),({"A"},{"B"}),({"Unknown"},{"A"})]
        out=m.paired_summary(pairs,t,min_species=1,min_pairs=2)
        self.assertEqual(out["n_candidates"],3)
        self.assertEqual(out["n_pairs"],2)
        self.assertEqual(out["mean_difference"],0.)
        self.assertEqual(out["n_positive"],1)
        self.assertEqual(out["n_negative"],1)

    def test_low_pair_support_keeps_coverage_but_withholds_outcome_summaries(self):
        out=m.paired_summary([({"B"},{"A"})],trait({"A":2.,"B":8.}),min_species=1)
        self.assertFalse(out["estimable"])
        self.assertEqual(out["n_pairs"],1)
        self.assertEqual(out["mean_rail_numeric_taxon_coverage"],1.)
        for key in ("mean_difference","median_difference","n_positive","n_negative","n_zero","difference_range"):
            self.assertIsNone(out[key],key)


if __name__=="__main__":
    unittest.main()
