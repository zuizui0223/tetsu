"""v333 synthetic regression checks for source-fixed geographic pairs and trait-label controls."""
import unittest,sys,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_EIVE_matched_station_taxon_null_v333 as m

SOURCE={"hot":8.,"mid":6.,"cold":4.}
PAIR={"BT":{"hot","mid"},"BI":{"cold","mid"},
      "ST":{"cold","mid"},"SI":{"cold","mid"}}
class TaxonRobustness(unittest.TestCase):
    def test_matched_CWM_exact_coefficients(self):
        c=m.coefficients([PAIR])
        self.assertAlmostEqual(sum(c[s]*SOURCE[s] for s in SOURCE),2)
        self.assertAlmostEqual(m.paired_mean([PAIR],SOURCE),2)
        self.assertAlmostEqual(m.numeric_observed([PAIR],SOURCE),2)
    def test_contrast_changes_when_outcome_taxon_removed(self):
        result=m.leave_one_taxon_out([PAIR],SOURCE)
        self.assertEqual(result["n_taxa_seen_in_matched_source_plots"],3)
        self.assertEqual(result["taxa_preserving_same_all_original_pairs_after_removal"],3)
        self.assertTrue(result["number_of_full_pair_taxon_deletions_making_difference_nonpositive"]>=1)
    def test_technical_equal_known_species_rarefaction_is_not_field_bootstrap(self):
        pairs=[PAIR]*3
        r=m.rarefy_known_taxa(pairs,SOURCE,draws=49,seed=32)
        self.assertEqual(r["n_species_equalized_diagnostic_draws"],49)
        self.assertEqual(r["n_source_original_pairs"],3)
        self.assertAlmostEqual(r["mean_across_technical_draws"],2)
    def test_global_null_preserves_fixed_source_incidence(self):
        c=m.coefficients([PAIR])
        r=m.shuffled_value_distribution(c,SOURCE,"global",29,seed=3)
        self.assertEqual(r["n_draws"],29)
        self.assertEqual(r["n_original_source_mapped_taxa_in_matched_comparison"],3)
        self.assertTrue(0<r["technical_two_sided_fraction_plus_one"]<=1)
    def test_genus_null_singletons_are_fixed(self):
        p={"BT":{"g alpha","h beta"},"BI":{"g beta","h beta"},
           "ST":{"g beta","h beta"},"SI":{"g beta","h beta"}}
        values={"g alpha":8.,"g beta":4.,"h beta":5.}
        c=m.coefficients([p])
        r=m.shuffled_value_distribution(c,values,"within_genus",19,seed=5)
        self.assertEqual(r["n_genera_with_at_least_two_species"],1)
        self.assertEqual(r["n_taxa_in_permutable_genera"],2)
    def test_rejected_empty_plot(self):
        with self.assertRaises(ValueError):
            m.coefficients([{"BT":set(),"BI":{"cold"},"ST":{"mid"},"SI":{"hot"}}])
    def test_no_observed_unmatched_source_sites_emitted(self):
        self.assertEqual(len(m.PANELS),6)
        self.assertEqual(m.REFERENCE,"Interior Plot 2")
    def test_quantiles_run_no_sample_intervals_are_field_CIs(self):
        r=m.rarefy_known_taxa([PAIR],SOURCE,draws=39,seed=10)
        self.assertTrue(r["technical_percentiles_are_NOT_independent_source_field_confidence_intervals"])
if __name__=="__main__":unittest.main()
