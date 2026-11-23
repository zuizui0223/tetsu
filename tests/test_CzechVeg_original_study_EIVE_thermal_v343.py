#!/usr/bin/env python3
"""Synthetic v343 safeguards only; not field effects."""
import random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_CzechVeg_original_study_EIVE_thermal_v343 import (
  normalized_binomial,plot_score,frozen_frame,report_stats,trait_null
)

class V343OriginalStudySourceTests(unittest.TestCase):
    def setUp(self):
        self.names=["cynodon dactylon","poa annua","setaria viridis",
                    "bromus japonicus","avenas sativa"]
        self.traits={n:(float(i+3),3) for i,n in enumerate(self.names)}
        self.sites={k:{"source_total_taxa":5,
            "exact_binomial_taxa":set(self.names)} for k in ("r1","r2","d1","d2")}
        self.blocks={("paper1","region1",2010,"area_band"):
            {"rail":["r1","r2"],"road":["d1","d2"]}}

    def test_strict_exact_binomial_normalization(self):
        self.assertEqual(normalized_binomial(" Cynodon   dactylon  "),"cynodon dactylon")
        self.assertIsNone(normalized_binomial("Cynodon dactylon subsp. dactylon"))
        self.assertIsNone(normalized_binomial("Cynodon"))
        self.assertIsNone(normalized_binomial(""))

    def test_plot_coverage_requires_three_known_and_half_source(self):
        score=plot_score(self.sites["r1"],self.traits,3,.5)
        self.assertEqual(score["matched_n"],5)
        self.assertAlmostEqual(score["coverage"],1)
        fail=plot_score(self.sites["r1"],
                        {"cynodon dactylon":(8,3),"poa annua":(2,3)},3,.5)
        self.assertIsNone(fail)

    def test_preserve_original_study_blocks_and_minimum_arm(self):
        selected,scored,qa=frozen_frame(self.blocks,self.sites,self.traits,3,.5)
        self.assertEqual(qa["n_selected_original_field_studies"],1)
        self.assertEqual(qa["n_selected_original_matched_blocks"],1)
        self.assertEqual(qa["n_selected_original_RAIL_text_plots"],2)
        self.assertEqual(qa["n_selected_original_ROAD_text_plots"],2)
        b2={("paper1","region1",2010,"area_band"):
            {"rail":["r1"],"road":["d1","d2"]}}
        _,_,q=frozen_frame(b2,self.sites,self.traits,3,.5)
        self.assertEqual(q["n_selected_original_matched_blocks"],0)

    def test_study_equal_original_difference_from_fixed_plots(self):
        selected,scored,qa=frozen_frame(self.blocks,self.sites,self.traits,3,.5)
        floats={k:v[0] for k,v in self.traits.items()}
        observed,paper=report_stats(selected,scored,floats)
        self.assertEqual(observed["n_original_studies"],1)
        self.assertAlmostEqual(observed["mean_study_equal_RAIL_minus_ROAD_EIVE"],0)
        self.assertEqual(observed["n_original_matched_blocks"],1)

    def test_complete_source_empty_trait_null_is_not_field_effect(self):
        _,scored,qa=frozen_frame(self.blocks,self.sites,self.traits,3,.5)
        selected={("paper1","region1",2010,"area_band"):{"rail":["r1","r2"],"road":["d1","d2"]}}
        vals={k:v[0] for k,v in self.traits.items()}
        observed,_=report_stats(selected,scored,vals)
        out=trait_null(selected,scored,vals,observed,"global",draws=11)
        self.assertEqual(out["draws"],11)
        self.assertEqual(out["two_sided_fraction_centered_plus_one"],1)
        self.assertTrue(out["original_source_incidence_and_trait_known_masks_fixed"])

if __name__=="__main__":unittest.main()
