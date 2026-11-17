import hashlib,math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from audit_Nov17_country_rail_vs_domestic_reference_robustness_v251 import (
    evaluate_comparison,git_blob,station_proximity
)

def comparison(delta=.3,extended=.1,maebetter=True,genus_lower=.1,recorded=65):
    return {
      "common":{"n":80,"recorded":recorded,
        "foreign":.6,"domestic":.6-delta,"delta":delta,
        "mae_foreign":.18,"mae_domestic":.23 if maebetter else .15,
        "genus_interval":[genus_lower,.5]},
      "extended":{"n":200,"recorded":105,
        "foreign":.5,"domestic":.5-extended,"delta":extended,
        "mae_foreign":.15,"mae_domestic":.18 if maebetter else .12},
      "delta_excluding_training_study_fixed_panel":{"CZ":.25,"SI":.2,"AD":.09}
    }

class TestNov17CountryHeldoutReliabilityAudit(unittest.TestCase):
    def test_strict_all_eight_exploratory_conditions(self):
        r=evaluate_comparison("IT_road",comparison())
        self.assertTrue(r["all_eight_conditions_met"])
        self.assertEqual(len(r["eight_exploratory_descriptive_reliability_conditions"]),8)
        self.assertEqual(r["n_source_original_training_paper_removal_scenarios"],3)
    def test_positive_rank_can_be_bad_frequency_calibration(self):
        r=evaluate_comparison("SI_urban",comparison(maebetter=False))
        self.assertFalse(r["all_eight_conditions_met"])
        self.assertTrue(r["eight_exploratory_descriptive_reliability_conditions"][
           "foreign_rank_beats_domestic_in_selected_common_taxa"])
        self.assertFalse(r["eight_exploratory_descriptive_reliability_conditions"][
           "foreign_absolute_frequency_MAE_better_in_common_taxa"])
    def test_missing_genus_support_and_low_observed_species_fail_without_false_positive(self):
        r=evaluate_comparison("BA_urban",comparison(genus_lower=-.1,recorded=18))
        self.assertFalse(r["all_eight_conditions_met"])
        self.assertFalse(r["eight_exploratory_descriptive_reliability_conditions"][
           "at_least_half_source_species_recorded_in_holdout"])
    def test_extended_country_target_can_reverse_rank_sign(self):
        r=evaluate_comparison("CZ_urban",comparison(delta=.4,extended=-.3))
        self.assertFalse(r["all_eight_conditions_met"])
        self.assertFalse(r["eight_exploratory_descriptive_reliability_conditions"][
           "foreign_rank_beats_domestic_in_extended_support"])
    def test_original_score_is_enforced_not_manual_fabricated_delta(self):
        row=comparison();row["common"]["delta"]=.4
        with self.assertRaises(ValueError):evaluate_comparison("IT_road",row)
    def test_git_blob_sha_matches_GitHub_original_blob_definition(self):
        b=b"hello\n"
        expected=hashlib.sha1(b"blob 6\0hello\n").hexdigest()
        self.assertEqual(git_blob(b),expected)
    def test_250m_station_proximity_not_a_zero_railway_homogenization_effect(self):
        base={"source_plots":95104,"source_station_records":66856,"bands_km":{
            "0.25":{"supported":[],"matched_stations":7},
            "1.0":{"supported":[{"country":"Z","stations":4,"cwm":.2221,"leda":{"n":0}}],"matched_stations":4},
            "2.0":{"supported":[{"country":"Z","stations":5,"cwm":0.0,"leda":{"n":5,"equal":.1}}],"matched_stations":5},
            "1.0_strict":{"supported":[],"matched_stations":1}}}
        r=station_proximity(base)
        self.assertEqual(r["0.25"]["n_countries_with_at_least_four_true_station_nearby_ref_pairs"],0)
        self.assertTrue(r["0.25"]["no_cwm_estimate_not_zero_effect"])
        self.assertAlmostEqual(r["1.0"]["TRY_cover_weighted_plantheight_seedmass_near_station_minus_other_mean_equal_country"],.2221)
    def test_source_support_requires_actual_site_units_not_only_species(self):
        with self.assertRaises(ValueError):evaluate_comparison("SK_urban",{
           **comparison(),"common":{**comparison()["common"],"n":3,"recorded":5}})
if __name__=="__main__":unittest.main()
