#!/usr/bin/env python3
"""Synthetic invariants, not biological source observations."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_Kashmir_Cynodon_transition_v340 import (
    CYNODON,zone_counts,technical_site_did,exact_mcnemar_double_baseline,
    source_species_amounts,source_sites,transition,quantile
)


class V340Test(unittest.TestCase):
    def setUp(self):
        self.names=[CYNODON,"poa annua","sorghum halepense"]
        self.rows=[
            ("A", {
               (2014,"Track"):{0:8.0,1:2.0},
               (2017,"Track"):{1:3.0},
               (2014,"Interior Plot 2"):{0:4.0,1:1.0},
               (2017,"Interior Plot 2"):{0:2.0,1:1.0},
            }),
            ("B",{
               (2014,"Track"):{0:5.0,1:5.0},
               (2017,"Track"):{0:10.0,1:2.0},
               (2014,"Interior Plot 2"):{1:1.0},
               (2017,"Interior Plot 2"):{0:3.0,1:1.0},
            }),
            ("Between C",{
               (2014,"Track"):{1:1.0},
               (2017,"Track"):{0:1.0,1:1.0},
               (2014,"Interior Plot 2"):{1:1.0},
               (2017,"Interior Plot 2"):{1:2.0},
            }),
            ("Between D",{
               (2014,"Track"):{0:1.0,1:2.0},
               (2017,"Track"):{0:1.0,1:1.0},
               (2014,"Interior Plot 2"):{0:1.0,1:2.0},
               (2017,"Interior Plot 2"):{0:2.0,1:1.0},
            })
        ]

    def test_reobservation_counts_distinguish_01_and_10(self):
        z=zone_counts(self.rows,"Track",self.names,{CYNODON})
        self.assertEqual(z["transition_2014_to_2017"],{"00":0,"01":1,"10":1,"11":2})
        self.assertAlmostEqual(z["2014_presence_to_2017_reobserved_fraction"],2/3)
        self.assertAlmostEqual(z["2014_absence_to_2017_new_record_fraction"],1)

    def test_both2014_positive_paired_mcnemar_correct(self):
        p=exact_mcnemar_double_baseline(
            self.rows,("Track","Interior Plot 2"),self.names,{CYNODON})
        self.assertEqual(p["n_original_sites_occupied_both_habitats_in_2014"],2)
        self.assertEqual(p["2017_interior_only_reobserved"],1)
        self.assertEqual(p["2017_both_reobserved"],1)
        self.assertTrue(p["not_estimable_direct_zone_difference_if_both2014_positive_sites_under5"])

    def test_site_did_repeats_same_sites_not_species_rows(self):
        res=technical_site_did(self.rows,("Track","Interior Plot 2"),
                    self.names,{CYNODON},draws=99,boot=99)
        self.assertEqual(res["n_same_original_sites"],4)
        # A Track -1 vs interior 0 -> -1; B Track0 vs interior+1 -> -1
        # C Track+1 vs interior0 -> +1; D both 0. mean -0.25
        self.assertAlmostEqual(res["observed_change_in_track_minus_interior_occupancy_fraction"],-0.25)
        self.assertEqual(res["n_negative_site_changes"],2)
        self.assertEqual(res["n_positive_site_changes"],1)

    def test_numeric_temporal_ratios_only_positive_both(self):
        res=source_species_amounts(self.rows,("Track","Interior Plot 2"),self.names)
        t=res["Track"]["same_source_plot_Cynodon_positive_both_years"]
        self.assertEqual(t["n"],2)
        self.assertAlmostEqual(t["geometric_mean_ratio_2017_div_2014"],math.sqrt(2))
        self.assertEqual(t["n_source_values_increased"],1)
        self.assertEqual(t["n_source_values_unchanged"],1)
        self.assertEqual(res["Track"]["year_stats"]["2014"]["n_positive_Cynodon_source_values"],3)
        self.assertEqual(res["Track"]["year_stats"]["2017"]["n_positive_Cynodon_source_values"],3)

    def test_missing_source_plot_is_excluded_not_absence(self):
        abund={}
        for i in range(31):
            site="Between "+str(i) if i<15 else "Station "+str(i)
            for y in (2014,2017):
                for z in ("Track","Interior Plot 2"):
                    abund[site,y,z]={0:1}
        # change one source site field to missing, excludes entire site
        abund["Between 0",2017,"Track"]=None
        sites,why=source_sites(abund,("Track","Interior Plot 2"))
        self.assertEqual(len(sites),30)
        self.assertEqual(why["missing_at_least_one_original_site_year_zone"],1)

    def test_defined_transition_and_percentiles(self):
        self.assertEqual(transition(False,True),"01")
        self.assertEqual(quantile([0,1,2],.5),1)


if __name__=="__main__":unittest.main()
