"""Baselga 2012 Jaccard decomposition and original-source gates; synthetic math only."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_rail_beta_turnover_nestedness_v74 import (
    baselga_jaccard_components,pair_beta,METRICS,FAMILIES,K,one_country)

class JaccardDecompositionTests(unittest.TestCase):
    def test_identical_sets_have_zero_beta(self):
        v=baselga_jaccard_components({"a":1,"b":1},{"a":1,"b":1})
        self.assertEqual(tuple(v[k] for k in METRICS),(0,0,0))
    def test_all_turnover_species_replacement(self):
        v=baselga_jaccard_components({"a":1,"b":1},{"a":1,"c":1})
        self.assertAlmostEqual(v["beta_jac"],2/3)
        self.assertAlmostEqual(v["beta_jtu"],2/3)
        self.assertAlmostEqual(v["beta_jne"],0)
    def test_species_list_nestedness(self):
        v=baselga_jaccard_components({"a":1},{"a":1,"b":1})
        self.assertAlmostEqual(v["beta_jac"],0.5)
        self.assertAlmostEqual(v["beta_jtu"],0)
        self.assertAlmostEqual(v["beta_jne"],0.5)
    def test_both_turnover_and_nestedness(self):
        v=baselga_jaccard_components({"a":1,"b":1,"c":1},{"a":1,"d":1})
        self.assertAlmostEqual(v["beta_jac"],.75)
        self.assertAlmostEqual(v["beta_jtu"],2/3)
        self.assertAlmostEqual(v["beta_jne"],1/12)
    def test_undefined_without_shared_or_replacement_is_not_named_nestedness(self):
        self.assertIsNone(baselga_jaccard_components({"a":1},{}))
    def test_not_expectation_of_beta(self):
        a={"a":15,"b":12};b={"a":10,"c":30}
        p=pair_beta(a,b)
        self.assertEqual(set(p),{"raw","rarefied20"})
        for source in FAMILIES:
            self.assertAlmostEqual(p[source]["beta_jac"],
               p[source]["beta_jtu"]+p[source]["beta_jne"])
        self.assertEqual(K,20)
    def test_invalid_probabilities_do_not_produce_ecological_claim(self):
        with self.assertRaisesRegex(ValueError,"Species inclusion"):
            baselga_jaccard_components({"a":-1},{"a":1})
    def test_frozen_country_cell_pairs_not_significance_replicates(self):
        result=one_country([],0)
        self.assertFalse(result["source_support"])
        self.assertIsNone(result["coefficients"]["full"]["woody"]["raw"]
                          ["partial_rho_rail_ratio_vs_dissimilarity"]["beta_jac"])
        self.assertEqual(result["n_dependent_common_50km_cell_pairs"],0)

if __name__=="__main__":unittest.main()
