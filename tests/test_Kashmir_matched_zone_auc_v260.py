"""Synthetic whole-site paired-plot AUC contrasts for v26.0."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Kashmir_same_sites_early_cover_AUC_track_vs_adjacent_v260 as m


class PairedRailwayCoverAUCTests(unittest.TestCase):
    @staticmethod
    def make(n=12,track=.8,reference=.6,noneligible=0):
        res={}
        for i in range(n):
            res[str(i)]={}
            for zone in m.ZONES:
                eligible=(i>=noneligible)
                score=track if zone=="Track" else reference
                res[str(i)][zone]={
                    "eligible_auc":eligible,
                    "source_2014_cover_AUC":score,
                    "n_original_2014_taxa":5,
                    "n_persistent_recorded_2017":2
                }
        return res

    def test_same_original_site_paired_difference(self):
        c=self.make()
        r=m.paired_summary(c,199,42)
        for zone in m.COMPARATORS:
            self.assertTrue(r[zone]["identified"])
            self.assertEqual(r[zone]["paired_original_site_count"],12)
            self.assertAlmostEqual(r[zone]["source_2014_cover_AUC_track_minus_reference_paired_mean"],.2)
            self.assertAlmostEqual(r[zone]["site_cluster_resampling_technical"]["q025_not_randomized_field_CI"],.2)
            self.assertAlmostEqual(r[zone]["leave_one_original_site_out"]["min"],.2)
            self.assertEqual(r[zone]["n_sites_track_auc_gt_reference_auc"],12)

    def test_equal_areas_are_not_assumed_equal_species_risk_sets(self):
        r=m.paired_summary(self.make(n=12,track=.5,reference=.5),199)
        self.assertEqual(r["Embankment"]["source_2014_cover_AUC_track_minus_reference_paired_mean"],0)
        self.assertEqual(r["Embankment"]["n_sites_track_auc_eq_reference_auc"],12)
        self.assertTrue(r["Embankment"]["no_randomized_zone_assignment_or_shared_plant_risk_set"])

    def test_insufficient_both_zone_site_support_marked_unknown(self):
        raw=self.make(n=12,noneligible=8)
        r=m.paired_summary(raw,199)
        self.assertFalse(r["Interior Plot 2"]["identified"])
        self.assertEqual(r["Interior Plot 2"]["paired_original_site_count"],4)
        self.assertNotIn("source_2014_cover_AUC_track_minus_reference_paired_mean",r["Interior Plot 2"])

    def test_site_cluster_not_taxon_event_bootstrap(self):
        raw=self.make(n=12)
        # One site's score reversed, original site remains one contribution.
        raw["0"]["Track"]["source_2014_cover_AUC"]=.2
        r=m.paired_summary(raw,199,1021)["Embankment"]
        self.assertAlmostEqual(r["source_2014_cover_AUC_track_minus_reference_paired_mean"],(11*.2-.4)/12)
        self.assertEqual(r["n_sites_track_auc_lt_reference_auc"],1)
        self.assertEqual(r["n_sites_track_auc_gt_reference_auc"],11)


if __name__=="__main__":unittest.main()
