#!/usr/bin/env python3
"""Synthetic label/temporal gates; no ecological effects asserted."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_CzechVeg_strict_HABITAT_EIVE_v344 import (
    subset_by_era,safe_effect)
from analyze_CzechVeg_original_study_EIVE_thermal_v343 import (
    frozen_frame,trait_data_to_positions
)

class V344StrictSourceTests(unittest.TestCase):
    def setUp(self):
        self.blocks={
          ("paper1","region1",1990,"4_16"):
             {"rail":["r1","r2"],"road":["d1","d2"]},
          ("paper2","region2",2010,"4_16"):
             {"rail":["r3","r4"],"road":["d3","d4"]},
          ("paper3","region3",2000,"4_16"):
             {"rail":["r5","r6"],"road":["d5","d6"]},
        }
        self.traits={"cynodon dactylon":(8.0,3),"poa annua":(3.0,3),
                     "bromus japonicus":(2.0,3),"setaria viridis":(7.0,3),
                     "avena sativa":(4.0,3)}
        self.plots={}
        for i in range(1,7):
            self.plots["r"+str(i)]={"source_total_taxa":5,"exact_binomial_taxa":set(self.traits)}
            self.plots["d"+str(i)]={"source_total_taxa":5,"exact_binomial_taxa":set(self.traits)}

    def test_era_splits_unchanged_observation_eligibility(self):
        x=subset_by_era(self.blocks,"before_2000")
        y=subset_by_era(self.blocks,"year_2000_onward")
        z=subset_by_era(self.blocks,"year_2010_onward")
        self.assertEqual(len(x),1)
        self.assertEqual(len(y),2)
        self.assertEqual(len(z),1)
        self.assertEqual(set(x).union(y),set(self.blocks))

    def test_no_study_means_when_source_outcomes_unavailable(self):
        a=safe_effect({}, {}, self.traits, frozenset())
        self.assertEqual(a["n_original_studies"],0)
        self.assertIn("inference_gate",a)

    def test_identical_source_arms_have_zero_outcome(self):
        selected,score,qa=frozen_frame(self.blocks,self.plots,self.traits,3,.5)
        full=safe_effect(selected,score,self.traits,frozenset())
        self.assertEqual(full["n_original_studies"],3)
        self.assertEqual(full["n_original_matched_blocks"],3)
        self.assertAlmostEqual(full["mean_study_equal_RAIL_minus_ROAD_EIVE"],0.0)
        self.assertIn("inference_gate",full)

    def test_C4_deletion_does_not_rematch_source_papers(self):
        a,b,qa=frozen_frame(self.blocks,self.plots,self.traits,3,.5)
        full=safe_effect(a,b,self.traits,frozenset())
        no_c4=safe_effect(a,b,self.traits,{"cynodon dactylon"})
        self.assertEqual(full["n_original_matched_blocks"],no_c4["n_original_matched_blocks"])

    def test_very_strict_coverage_can_fail_without_imputing_traits(self):
        insufficient={"r1":{"source_total_taxa":5,
                          "exact_binomial_taxa":{"poa annua"}},
                      "r2":{"source_total_taxa":5,
                          "exact_binomial_taxa":{"poa annua"}},
                      "d1":{"source_total_taxa":5,
                          "exact_binomial_taxa":{"poa annua"}},
                      "d2":{"source_total_taxa":5,
                          "exact_binomial_taxa":{"poa annua"}}}
        s,b,q=frozen_frame({next(iter(self.blocks)):self.blocks[next(iter(self.blocks))]},
                           insufficient,self.traits,3,.5)
        self.assertEqual(q["n_selected_original_matched_blocks"],0)

if __name__=="__main__":unittest.main()
