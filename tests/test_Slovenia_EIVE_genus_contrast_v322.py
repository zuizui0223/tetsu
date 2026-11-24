"""v322 synthetic tests: check additive taxonomic thermal decomposition."""
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Slovenia_EIVE_genus_contrast_v322 as m

NAMES={"A":"genus alpha","B":"genus beta",
       "C":"other gamma","D":"other delta"}
TRAITS={
 "genus alpha":{"position":8.0,"width":2.0,"systems":3},
 "genus beta":{"position":6.0,"width":4.0,"systems":3},
 "other gamma":{"position":3.0,"width":3.0,"systems":3},
 "other delta":{"position":5.0,"width":5.0,"systems":3},
}
def pair(bs,bt,site="s1",section="r1",side="L"):
    return {"site":site,"section":section,"side":side,
            "BS":set(bs),"BT":set(bt)}

class GenusThermalDecomposition(unittest.TestCase):
    def test_species_residual_can_explain_total_without_genus_replacement(self):
        p=[pair(["A","D"],["B","C"])]
        res=m.contrasts(p,NAMES,TRAITS,"position",3,0)
        self.assertAlmostEqual(res["total_location_mean_delta"],2)
        self.assertAlmostEqual(res["between_genus_location_mean_delta"],0)
        self.assertAlmostEqual(res["within_genus_trait_residual_location_mean_delta"],2)
    def test_genus_contrast_can_explain_total_without_species_residual(self):
        p=[pair(["A","B"],["C","D"])]
        res=m.contrasts(p,NAMES,TRAITS,"position",3,0)
        self.assertAlmostEqual(res["total_location_mean_delta"],3)
        self.assertAlmostEqual(res["between_genus_location_mean_delta"],3)
        self.assertAlmostEqual(res["within_genus_trait_residual_location_mean_delta"],0)
    def test_two_track_sides_count_one_field_site(self):
        p=[pair(["A","D"],["B","C"],side="L"),
           pair(["A","D"],["B","C"],side="R"),
           pair(["B","C"],["A","D"],site="s2",section="r2")]
        res=m.contrasts(p,NAMES,TRAITS,"position",3,0)
        self.assertEqual(res["eligible_field_sites"],2)
        self.assertEqual(res["eligible_track_sides"],3)
        self.assertAlmostEqual(res["total_location_mean_delta"],0)
    def test_within_genus_permutation_reproduces_observed_stat(self):
        p=[pair(["A","D"],["B","C"])]
        res=m.contrasts(p,NAMES,TRAITS,"position",3,0)
        perm=m.within_genus_permutation(res,n_draws=19,seed=12)
        self.assertAlmostEqual(perm["actual_within_genus_component"],2)
        self.assertEqual(perm["n_draws"],19)
        self.assertEqual(perm["n_genera_permuted"],2)
        self.assertGreaterEqual(perm["conditional_two_sided_fraction_with_plus_one"],1/20)
    def test_coverage_gate_prevents_nonmatched_sample_leakage(self):
        p=[pair(["A","D","MISSING"],["B","C"])]
        ok=m.contrasts(p,NAMES,TRAITS,"position",3,0)
        self.assertEqual(ok["eligible_field_sites"],1)
        with self.assertRaisesRegex(ValueError,"No eligible"):
            m.contrasts(p,NAMES,TRAITS,"position",3,.7)
    def test_position_width_separate_compatible_decisions(self):
        p=[pair(["A","D"],["B","C"])]
        position=m.contrasts(p,NAMES,TRAITS,"position",3,0)
        width=m.contrasts(p,NAMES,TRAITS,"width",3,0)
        self.assertAlmostEqual(position["total_location_mean_delta"],2)
        self.assertAlmostEqual(width["total_location_mean_delta"],0)
        self.assertEqual(position["represented_genera"],2)
        self.assertEqual(position["species_in_multi_species_genera"],4)

if __name__=="__main__": unittest.main()
