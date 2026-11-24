#!/usr/bin/env python3
"""v365 source Heldout original line network-distance incremental predictive maths."""
import importlib.util,math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from evaluate_Slovenia_network_distance_line_heldout_v365 import (
    MODELS,MODELS_NO_RICH,LAM,
    source_degree_pair_weights,features,fit_predict,line_holdout_fit
)
HAS_NUMPY=importlib.util.find_spec("numpy") is not None

@unittest.skipUnless(HAS_NUMPY,"numerical source predicted heldout Jaccard tests require numpy")
class V365NetworkVsGeographyTests(unittest.TestCase):
    def record(self,i,j,net,geo,y,section_i="A",section_j="B",sameyear=True):
        return {"i":i,"j":j,
                "source_original_geodetic_km":geo,
                "source_true_rail_GERS_network_km":net,
                "same_original_survey_year":sameyear,
                "same_original_track_layout":True,
                "same_original_author_section":section_i==section_j,
                "original_author_section_pair":(section_i,section_j),
                "site_source_field_richness":(10,10),
                "response":{"BS_and_BT":y,"BS":y,"BT":y}}
    def test_true_network_length_improves_holdout_when_real_signal_is_network(self):
        import numpy as np
        train=[self.record(i%16,(i+1)%16,net=1.+i,
               geo=3.+(i%5),y=max(0.05,0.9-.012*i))
               for i in range(60)]
        test=[self.record(i%16,(i+1)%16,net=1.+i,
               geo=3.+(i%5),y=max(0.05,0.9-.012*i))
               for i in range(60,80)]
        pred0=fit_predict(train,test,"BS_and_BT",MODELS["M0_geo_year_layout_richness"])
        pred1=fit_predict(train,test,"BS_and_BT",MODELS["M1_M0_plus_TRUE_network_distance"])
        target=np.asarray([x["response"]["BS_and_BT"] for x in test])
        self.assertTrue(np.all(np.isfinite(pred0)) and np.all(np.isfinite(pred1)))
        self.assertTrue(np.all(pred1>=0) and np.all(pred1<=1))
        self.assertLess(np.mean((target-pred1)**2),np.mean((target-pred0)**2))

    def test_source_original_pair_weights_do_not_sum_to_pairs_as_independent_replicates(self):
        rows=[self.record(0,i,net=1.+i,geo=2.,y=0.4) for i in range(1,9)]
        w=source_degree_pair_weights(rows)
        self.assertAlmostEqual(sum(w),1.)
        self.assertTrue(all(x>0 for x in w))

    def test_same_frozen_models_m0_m1_m2_m3_and_no_site_richness_source_sensitivity(self):
        self.assertEqual(LAM,.05)
        self.assertEqual(len(MODELS),4)
        self.assertEqual(len(MODELS_NO_RICH),4)
        self.assertEqual(list(MODELS)[0],"M0_geo_year_layout_richness")
        self.assertIn("log_physical_network",MODELS["M1_M0_plus_TRUE_network_distance"])
        self.assertNotIn("log_physical_network",MODELS["M0_geo_year_layout_richness"])
        self.assertIn("same_author_section",MODELS["M3_M2_plus_TRUE_network_distance"])
        self.assertNotIn("mean_log_site_richness",MODELS_NO_RICH["M0_geo_year_layout"])

    def test_source_botany_taxonomic_occupancy_not_covariate_when_sample_richness_omitted(self):
        p=self.record(0,1,net=24,geo=12,y=.3)
        f=features(p)
        self.assertIn("log_physical_network",f)
        self.assertIn("log_geodesic_squared",f)
        self.assertAlmostEqual(f["log_geodesic"],math.log1p(12))
        self.assertTrue(set(MODELS_NO_RICH["M1_geo_plus_TRUE_network_distance"])<=set(f))

if __name__=="__main__":unittest.main()
