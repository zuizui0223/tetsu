"""v15.9 synthetic fold-summary integrity; actual results read only from authenticated artifacts."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Czech_Slovenia_original_railroad_CV_fold_influence_v159 import fold_summary,run

class OriginalAuthorCVFoldInfluenceMath(unittest.TestCase):
    def test_one_bad_railway_stretch_can_reverse_equal_stretch_source_average(self):
        # SYNTHETIC numbers, not Czech vegetation observations:
        fake=[{"section":str(i),"mse_plus_road":1.,
               "mse_plus_road_and_rail":2. if i==0 else .8,
               "n_holdout":10} for i in range(5)]
        x=fold_summary(fake,"mse_plus_road","mse_plus_road_and_rail",
                       "n_holdout","section",5,-4.0,"synthetic")
        self.assertAlmostEqual(x["source_section_equal_weight_rail_increment_MSE_gain_pct"],-4.)
        loo=x["source_leave_one_PREDICTION_FOLD_OUT_influence"]
        self.assertEqual(loo["n_source_fold_deletions"],5)
        self.assertEqual(loo["n_original_fold_deletions_with_positive_rail_increment"],1)
        self.assertAlmostEqual(loo["maximum_rail_gain_pct"],20.)
        self.assertTrue(loo["NOT_original_retraining_after_removing_source_study_section"])

    def test_reproduce_weighted_alternative_but_not_invent_independent_pair_count(self):
        fake=[{"section":"A","mse_plus_road":1.,
               "mse_plus_road_and_rail":.5,"n_holdout":20},
              {"section":"B","mse_plus_road":1.,
               "mse_plus_road_and_rail":1.5,"n_holdout":2}]
        x=fold_summary(fake,"mse_plus_road","mse_plus_road_and_rail",
                       "n_holdout","section",2,0.,"synthetic")
        self.assertAlmostEqual(x["source_section_equal_weight_rail_increment_MSE_gain_pct"],0.)
        self.assertGreater(x["source_test_pair_weighted_rail_increment_MSE_gain_pct"],0)
        self.assertTrue(x["folds_and_pairs_NOT_independent_field_studies"])
        self.assertEqual(x["source_original_test_pair_weights_sum_DUPLICATED_ACROSS_FOLDS"],22)

    def test_nonunique_heldout_stretches_fail_closed(self):
        fake=[{"section":"A","mse_plus_road":1.,
               "mse_plus_road_and_rail":.9,"n_holdout":2}]*2
        with self.assertRaisesRegex(ValueError,"fold count"):
            fold_summary(fake,"mse_plus_road","mse_plus_road_and_rail",
                         "n_holdout","section",2,10.,"synthetic")

    def test_wrong_source_publisher_report_versions_fail(self):
        with self.assertRaisesRegex(ValueError,"Slovenia original"):
            run({"version":"fake"},{"version":"fake"})

if __name__=="__main__":unittest.main()
