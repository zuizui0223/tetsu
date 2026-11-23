#!/usr/bin/env python3
"""Synthetic v345 study/block/identical plot source-selection tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_CzechVeg_common_paper_label_thermal_v345 import (
    restrict_papers,restrict_keys,paper_comparison,same_original_plots)
from analyze_CzechVeg_original_study_EIVE_thermal_v343 import frozen_frame


class V345SourceOverlapTests(unittest.TestCase):
    def setUp(self):
        self.key1=("paper1","reg1",1990,"4_16")
        self.key2=("paper2","reg2",2010,"4_16")
        self.strict={self.key1:{"rail":["r1","r2"],"road":["d1","d2"]}}
        self.broad={
           self.key1:{"rail":["r1","r2","r3"],"road":["d1","d2","d3"]},
           self.key2:{"rail":["r4","r5"],"road":["d4","d5"]}}
        self.traits={
           "cynodon dactylon":(8.0,3),"poa annua":(3.0,3),
           "bromus japonicus":(2.0,3),"setaria viridis":(7.0,3),
           "avena sativa":(4.0,3)}
        keys=set(x for v in self.broad.values() for side in ("rail","road") for x in v[side])
        self.plots={x:{"source_total_taxa":5,"exact_binomial_taxa":set(self.traits)}
                    for x in keys}

    def test_paper_subset_preserves_original_study_blocks(self):
        d=restrict_papers(self.broad,{"paper1"})
        self.assertEqual(set(d),{self.key1})

    def test_shared_original_block_keys_only(self):
        d=restrict_keys(self.broad,set(self.strict))
        self.assertEqual(set(d),{self.key1})

    def test_matched_study_comparison_does_not_pseudoreplicate_blocks(self):
        ss,sc,_=frozen_frame(self.strict,self.plots,self.traits,3,.5)
        bs,bc,_=frozen_frame(self.broad,self.plots,self.traits,3,.5)
        x=paper_comparison(ss,sc,restrict_papers(bs,{"paper1"}),bc,self.traits)
        self.assertEqual(x["source_papers_present_in_both_estimates"],1)
        self.assertAlmostEqual(x["same_paper_equal_mean_broad_minus_strict"],0)

    def test_same_source_plots_means_are_identical(self):
        ss,sc,_=frozen_frame(self.strict,self.plots,self.traits,3,.5)
        bs,bc,_=frozen_frame(self.broad,self.plots,self.traits,3,.5)
        v=same_original_plots(ss,bs,sc,bc,self.traits)
        self.assertEqual(v["source_S2_blocks_both_labels_have_identical_RAIL_ROAD_plot_IDs"],1)
        self.assertTrue(v["same_actual_plant_plot_temperature_estimates_exactly_identical"])
        self.assertEqual(v["n_original_source_plots_RAIL"],2)
        self.assertEqual(v["n_original_source_plots_ROAD"],2)

    def test_no_shared_study_results_missing_not_zero(self):
        ss,sc,_=frozen_frame(self.strict,self.plots,self.traits,3,.5)
        bs,bc,_=frozen_frame(self.broad,self.plots,self.traits,3,.5)
        e=paper_comparison(ss,sc,restrict_papers(bs,{"nothing"}),bc,self.traits)
        self.assertIsNone(e["same_paper_equal_mean_broad_minus_strict"])
        self.assertEqual(e["source_papers_present_in_both_estimates"],0)


if __name__=="__main__":
    unittest.main()
