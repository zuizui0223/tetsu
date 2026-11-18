"""v25.9 synthetic source-cover ordinal rank vs later-record persistence QA."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Kashmir_2014_species_cover_persistence_2017_v259 as m


class KashmirCoverPersistenceTests(unittest.TestCase):
    def test_cover_auc_is_rank_invariant_under_unit_scaling(self):
        self.assertEqual(m.auc_of_early_covers([2,3],[1,1]),1.)
        self.assertEqual(m.auc_of_early_covers([20,30],[10,10]),1.)
        self.assertEqual(m.auc_of_early_covers([1],[2]),0.)
        self.assertEqual(m.auc_of_early_covers([2],[2]),.5)
        self.assertEqual(m.auc_of_early_covers([2,4],[3,5]),.25)
        with self.assertRaises(ValueError):
            m.auc_of_early_covers([],[3])

    def test_risk_set_excludes_new_taxa_from_cover_predictor(self):
        x=m.site_early_persistence({1:2.,2:20.,3:1.},{2,4})
        self.assertEqual(x["n_original_2014_taxa"],3)
        self.assertEqual(x["n_persistent_recorded_2017"],1)
        self.assertEqual(x["n_not_recorded_2017"],2)
        self.assertEqual(x["n_newly_recorded_2017"],1)
        self.assertEqual(x["source_2014_cover_AUC"],1.)
        self.assertEqual(x["raw_ranks_all"],[1.,2.,20.])

    def test_not_eligible_if_all_lost_or_all_retained(self):
        self.assertFalse(m.site_early_persistence({1:1.},{2})["eligible_auc"])
        self.assertFalse(m.site_early_persistence({1:1.},{1,2})["eligible_auc"])

    def test_site_equal_not_taxon_event_pseudoreplication(self):
        early=[
            ({1:1.,2:2.,3:5.},{3}),
            ({1:1.,2:2.,3:5.},{3}),
            ({1:1.,2:2.,3:5.},{3}),
            ({1:1.,2:2.,3:5.},{3}),
            ({1:1.,2:2.,3:5.},{3}),
        ]
        rec={str(i):m.site_early_persistence(a,b) for i,(a,b) in enumerate(early)}
        out=m.source_stratum_analysis(rec,499,42)
        self.assertTrue(out["identified"])
        self.assertEqual(out["n_source_spatial_sites"],5)
        self.assertEqual(out["n_auc_eligible_original_sites"],5)
        self.assertEqual(out["n_2014_recorded_species_site_events"],15)
        self.assertAlmostEqual(out["site_equal_source_2014_cover_retained_vs_lost_AUC"],1.)
        self.assertAlmostEqual(out["null_label_within_source_plot_mean"],.5,delta=.07)
        self.assertAlmostEqual(out["leave_one_original_source_site_out"]["min"],1.)

    def test_no_four_valid_sites_no_auc_claim(self):
        r={str(i):m.site_early_persistence({1:1.},{1} if i<3 else {2}) for i in range(6)}
        v=m.source_stratum_analysis(r,99,1)
        self.assertFalse(v["identified"])
        self.assertNotIn("site_equal_source_2014_cover_retained_vs_lost_AUC",v)


if __name__=="__main__":unittest.main()
