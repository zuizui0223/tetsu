import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_original_GloNAF_matrix_v384 import compare
from audit_v382_heldout_continent_influence_v384b import summarize
class TestV384(unittest.TestCase):
    def test_set_comparison(self):
        v=compare({"a","b"},{"b","c","d"})
        self.assertEqual(v["shared"],1);self.assertEqual(v["direct_only"],1)
        self.assertEqual(v["matrix_only"],2);self.assertEqual(v["Jaccard"],.25)
    def test_fold_score_deletion(self):
        f=[{"heldout_continent":c,"delta_MSE_positive_means_rail_added_prediction_helped":d,
            "baseline_test_blockweighted_MSE":1.} for c,d in [("Asia",-.1),("Africa",.4),("Europe",.2)]]
        v={"base_road__to__both":{"source_independent_continent_holdouts":f,"mean_holdout_continent_delta_MSE":1/6}}
        p=dict(width_m=100000,min_land_fraction=.5,n_matched=3,n_pairs=2,richness_prediction=v,jaccard_prediction=v)
        out=summarize(dict(version="v382-original-world-plant-richness-and-Jaccard-road-tree-falsification",panels=[p]))["panels"][0]["richness_prediction"]
        self.assertAlmostEqual(out["leave_one_continent_score_out_min"],.05)
        self.assertAlmostEqual(out["leave_one_continent_score_out_max"],.3)
if __name__=="__main__":unittest.main()
