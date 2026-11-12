"""v8.4 tests for separately sourced FAF crop-vs-resource economic OD
comparison. Synthetic only; no botanical causation or 'seed-free' label.
"""
import math,sys,unittest
from hashlib import sha256
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_FAF6_crop_vs_resource_rail_floristic_v84 import (
  COMMODITIES,JOINT_CONTROLS,
  extract_OD,symmetric_OD,verify_matrix_marginals,
  augment_base_pairs,one_cohort,same_pair_cohort,partial_rank_pair,
  assert_original_FAF_integrity)

def network():
    return {
      ("crop_sctg02_03","RAIL","01","02"):10.,
      ("crop_sctg02_03","RAIL","02","01"):3.,
      ("crop_sctg02_03","TRUCK","01","02"):22.,
      ("resource_sctg10_12","RAIL","01","02"):5.,
      ("resource_sctg10_12","TRUCK","02","01"):9.,
      ("resource_sctg10_12","TRUCK","01","02"):6.,
    }

def source_marginals():
    prev={"data":{"national_trade_type_x_commodity_group_2022":{"domestic":{
      "crop_codes_02_03":{"rail_thousand_short_tons":13.,
        "truck_thousand_short_tons":22.},
      "resource_codes_10_11_12":{"rail_thousand_short_tons":5.,
        "truck_thousand_short_tons":15.}}}}}
    fa80={"mode_networks":{
      "crop_sctg02_03":{
          "RAIL":{"total_2022_estimated_thousand_short_tons":13.},
          "TRUCK":{"total_2022_estimated_thousand_short_tons":22.}},
      "resource_sctg10_12":{
          "RAIL":{"total_2022_estimated_thousand_short_tons":5.},
          "TRUCK":{"total_2022_estimated_thousand_short_tons":15.}}},
      "raw_official_state_OD_2022_crop_estimated_tonnage":{
         "RAIL":[{"o":"01","d":"02","estimated_thousand_short_tons":10},
                 {"o":"02","d":"01","estimated_thousand_short_tons":3}],
         "TRUCK":[{"o":"01","d":"02","estimated_thousand_short_tons":22}]},
    }
    return prev,fa80

class FAFResourceAsComparison(unittest.TestCase):
    def test_comparator_is_separate_cargo_category_not_sterile(self):
        self.assertEqual(COMMODITIES,("crop_sctg02_03","resource_sctg10_12"))
        self.assertIn("log1p_crop_truck_OD",JOINT_CONTROLS)
        self.assertIn("log1p_resource_truck_OD",JOINT_CONTROLS)
        self.assertIn("log1p_resource_rail_strength_geomean",JOINT_CONTROLS)
    def test_actual_OD_separated_by_mode_and_commodity(self):
        flow=network()
        crop=extract_OD(flow,"crop_sctg02_03","RAIL")
        resource=extract_OD(flow,"resource_sctg10_12","RAIL")
        self.assertEqual(symmetric_OD(crop,"01","02"),13)
        self.assertEqual(symmetric_OD(resource,"01","02"),5)
        self.assertEqual(symmetric_OD(resource,"02","02"),0)
        with self.assertRaisesRegex(ValueError,"predeclared"):
            extract_OD(flow,"viable_weed_seeds","RAIL")
    def test_authentic_2022_resource_and_crop_marginals_must_both_conserve(self):
        f=network()
        mat={g:{mode:extract_OD(f,g,mode) for mode in ("RAIL","TRUCK")}
             for g in COMMODITIES}
        a,b=source_marginals()
        self.assertTrue(verify_matrix_marginals(mat,a,b))
        mat["resource_sctg10_12"]["RAIL"]["01","02"]+=1
        with self.assertRaisesRegex(ValueError,"marginals"):
            verify_matrix_marginals(mat,a,b)
    def test_symmetric_floristic_pair_uses_two_direction_tonnage_but_not_route(self):
        flow=network()
        mats={g:{mode:extract_OD(flow,g,mode) for mode in ("RAIL","TRUCK")}
              for g in COMMODITIES}
        base=[{"a":"01","b":"02",
                "rail_estimated_thousand_short_tons":13,
                "truck_estimated_thousand_short_tons":22,
                "responses":{"rarefied20_expected_inclusion_overlap":.35,
                             "raw_recorded_GRIIS_Jaccard":.2},
                "covars":{"log_state_centroid_geo_distance":3.,
                          "log_allplant_record_geomean":8.,
                          "log_allplant_record_ratio":0.,
                          "log_state_land_area_geomean":12.,
                          "log_state_land_area_ratio":.1,
                          "same_census_region":1,
                          "same_census_division":1}}],
        out=augment_base_pairs(base,mats)
        self.assertEqual(out[0]["crop_rail"],13)
        self.assertEqual(out[0]["resource_rail"],5)
        self.assertEqual(out[0]["resource_truck"],15)
        self.assertIn("log1p_resource_truck_OD",out[0]["covars"])
        self.assertIn("log1p_crop_rail_strength_geomean",out[0]["covars"])
    def test_per_cargo_cohort_selected_without_botanical_response(self):
        rr=[{"a":"01","b":"02","crop_rail":1.,"crop_truck":2.,"resource_rail":0,
             "responses":{}},
            {"a":"01","b":"03","crop_rail":0.,"crop_truck":2.,"resource_rail":7,
             "responses":{}},
            {"a":"02","b":"03","crop_rail":3.,"crop_truck":0.,"resource_rail":5,
             "responses":{}}]
        self.assertEqual(len(same_pair_cohort(rr,"all_49_state_pairs")),3)
        self.assertEqual(len(same_pair_cohort(rr,"crop_rail_positive")),2)
        self.assertEqual(len(same_pair_cohort(rr,"crop_rail_and_truck_positive")),1)
        self.assertEqual(len(same_pair_cohort(rr,"both_crop_and_resource_rail_positive")),1)
        with self.assertRaisesRegex(ValueError,"Unfrozen"):
            same_pair_cohort(rr,"favorable_sign_only")
    def test_not_estimable_on_small_dependent_state_pairs(self):
        rr=[{"a":"01","b":"02","crop_rail":1,
             "resource_rail":2,"responses":{"rarefied20_expected_inclusion_overlap":.5},
             "covars":{}}]
        self.assertIsNone(partial_rank_pair(rr,"rarefied20_expected_inclusion_overlap",
                                             "crop_rail",JOINT_CONTROLS))
        c=one_cohort(rr)
        self.assertFalse(c["sufficient_independent_state_source_opportunity_not_independent_ecological_tests"])
        self.assertIsNone(c["outcome_source"]["rarefied20_expected_inclusion_overlap"]
                          ["crop_rail"]["joint_adjusted_partial_rank_rho"])
    def test_original_official_hash_and_state_universe_validated(self):
        raw=b"FAF6-2022-mock-bytes"
        sha=sha256(raw).hexdigest()
        old={"complete":True,"official_published_CSV_ZIP_sha256":sha,
             "data":{"n_original_official_FAF6_state_OD_commodity_mode_rows":599529}}
        fa80={"complete":True,"source_hash_verified_vs_prior_v79":True,
              "source_trace":{"previous_verified_v79_total_and_modes_reconciled":True,
                              "original_source_rows":599529}}
        self.assertTrue(assert_original_FAF_integrity(raw,old,fa80,599529,51))
        with self.assertRaisesRegex(ValueError,"drift"):
            assert_original_FAF_integrity(raw+b"x",old,fa80,599529,51)

if __name__=="__main__":unittest.main()
