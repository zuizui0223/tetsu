#!/usr/bin/env python3
import random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Kashmir_EIVE_species_reobservation_v341 import (
    cohort,site_trait_split,panel_stats,shuffle_traits,null_diagnostic,PRIMARY,CYNODON
)


class GenuineLongitudinalTraitTests(unittest.TestCase):
    def setUp(self):
        self.eive={CYNODON:8.0,"poa annua":2.0,"setaria viridis":6.0,
                   "avena sativa":3.0,"sorghum halepense":7.0}
        self.source=[
          {"Track":[(CYNODON,False,10.0),("poa annua",True,1.0),("setaria viridis",False,3.0)],
           "Interior Plot 2":[(CYNODON,True,12.0),("poa annua",False,1.0),("setaria viridis",True,9.0)]},
          {"Track":[(CYNODON,False,8.0),("avena sativa",True,2.0),("sorghum halepense",False,4.0)],
           "Interior Plot 2":[(CYNODON,True,6.0),("avena sativa",False,1.0),("sorghum halepense",True,3.0)]}
        ]

    def test_retained_minus_lost_sign_within_same_plot(self):
        a=site_trait_split(self.source[0]["Track"],self.eive)
        b=site_trait_split(self.source[0]["Interior Plot 2"],self.eive)
        self.assertAlmostEqual(a["observed_thermal_alignment_kept_minus_lost"],2-7)
        self.assertAlmostEqual(b["observed_thermal_alignment_kept_minus_lost"],7-2)
        self.assertLess(a["observed_thermal_alignment_kept_minus_lost"]-
                        b["observed_thermal_alignment_kept_minus_lost"],0)

    def test_2014_observed_species_not_2017_newcomers(self):
        z=site_trait_split([(CYNODON,False,10.0),("poa annua",True,2.0)],
                           self.eive)
        self.assertEqual(z["n_2014_EIVE_matched"],2)
        self.assertEqual(z["n_2017_reobserved"],1)
        self.assertEqual(z["n_2017_not_reobserved"],1)

    def test_exclusion_and_missing_traits_yield_none_instead_of_zero(self):
        z=site_trait_split([(CYNODON,False,10.0),("poa annua",True,2.0)],
                           self.eive,{CYNODON})
        self.assertIsNone(z["observed_thermal_alignment_kept_minus_lost"])
        self.assertEqual(z["n_2014_EIVE_matched"],1)

    def test_same_original_sites_pairing(self):
        obs=panel_stats(self.source,PRIMARY,self.eive,include_extra=True)
        self.assertEqual(obs["source_complete_site_count"],2)
        self.assertEqual(obs["n_informative_temperature_alignment_site_pairs"],2)
        self.assertLess(obs["mean_matched_track_minus_interior_thermal_persistence_alignment"],0)
        self.assertEqual(obs["zone_stats"]["Track"]["2014_species_site_record_count"],6)
        self.assertEqual(obs["zone_stats"]["Track"]["same_species_same_site_reobservation_count"],2)

    def test_shuffles_source_taxon_trait_values_without_changing_multiset(self):
        p=shuffle_traits(self.eive,"all_source_matched_taxa",random.Random(341))
        self.assertEqual(sorted(p.values()),sorted(self.eive.values()))
        g=shuffle_traits(self.eive,"within_genus_only",random.Random(341))
        self.assertEqual(g["poa annua"],self.eive["poa annua"])
        self.assertEqual(g[CYNODON],self.eive[CYNODON])

    def test_under10_site_matched_gate_stops_label_test(self):
        obs=panel_stats(self.source,PRIMARY,self.eive)
        out=null_diagnostic(self.source,PRIMARY,self.eive,frozenset(),obs,draws=19)
        self.assertEqual(out["status"],"under10_matched_informative_original_sites_not_inferential")


if __name__=="__main__":unittest.main()
