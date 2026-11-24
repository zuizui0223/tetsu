"""v376c independent worldwide GHSL settlement-adjusted railkm plant macro model tests."""
import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_global_original_railkm_biodiversity_GHSL_control_v376c as m

def testentry(key,cont="Europe",built=0):
 return {"key":key,"continent":cont,"log_species":3,
         "log_observations":5,"log_observers":2,"bio1":12,"bio1_sq":144,
         "abs_lat":40,"abs_lat_sq":1600,
         "block":(1,2),
         "rail_km":2,"log_rail_km":math.log1p(2),
         "species":19,"observations":100,"observers":12}
class SourceTrueGlobalBuiltup(unittest.TestCase):
 def test_source_DOI_and_units(self):
  self.assertTrue(m.BUILT_DOI.startswith("10.2905/"))
  self.assertEqual(m.RICH_BUILT,"log1p_1000_times_original_2020_builtup_fraction")
 def test_GHSL_positive_join_to_source_botanical_land_cell(self):
  rows=[testentry((1,2)),testentry((2,2))]
  a,qa=m.inject_built(rows,{(1,2):0.1,(2,2):0.0})
  self.assertEqual(len(a),2)
  self.assertAlmostEqual(a[0][m.RICH_BUILT],math.log1p(100))
  self.assertEqual(a[1][m.RICH_BUILT],0.)
  self.assertEqual(qa["source_qualified_true_zero_builtup"],1)
 def test_GHSL_missing_not_zero_source(self):
  rows=[testentry((1,2)),testentry((2,2))]
  a,qa=m.inject_built(rows,{(1,2):0.2})
  self.assertEqual(len(a),1)
  self.assertEqual(qa["source_quality_ok_plant_cell_which_lacks_original_GHSL2020_data"],1)
 def test_impossible_built_area_fails_closed(self):
  with self.assertRaises(ValueError):
   m.inject_built([testentry((1,2))],{(1,2):-1})
 def test_GHSL_100km_200km_same_parent_plant_frame(self):
  self.assertEqual(m.GRID_WIDTHS,(100000,200000))
  self.assertEqual(len(m.BETA_BUILT),2)
 def test_original_source_v371_corrected_before_controls(self):
  self.assertAlmostEqual(m.ACTUAL_V371_CORRECTED[100000]["richness"],0.004463745492835951)
  self.assertAlmostEqual(m.ACTUAL_V371_CORRECTED[200000]["jaccard"],-0.00006120015031440586)
 def test_result_source_mean_not_claim_ecological_causal_evidence(self):
  a={"source_independent_continent_holdouts":[
   {"heldout_continent":"Europe","baseline_test_blockweighted_MSE":2.0,
    "railkm_extended_test_blockweighted_MSE":1.2}]}
  b={"source_independent_continent_holdouts":[
   {"heldout_continent":"Europe","baseline_test_blockweighted_MSE":1.0,
    "railkm_extended_test_blockweighted_MSE":.9,
    "n_original_heldout_source_records":42}]}
  out=m.model_compare(a,b)
  self.assertAlmostEqual(out["mean_six_continent_equal_source_MSE"]["GERS_rail_lift_after_JRC_settlement"],.1)
  self.assertEqual(out["n_continents_where_railkm_still_helps_after_builtup"],1)
  self.assertTrue(out["no_causal_rail_or_plant_source_claim"])
 def test_urban_ridge_added_only_once_and_not_again_in_rail_increment(self):
  from unittest.mock import patch
  called=[]
  def fake_holdout(rows,base,added,*args):
   called.append((tuple(base),tuple(added)))
   return {"n_original_continents_eligible_as_heldout_folds":6}
  with patch.object(m,"continent_holdout",side_effect=fake_holdout):
   m.fit_same_study_continent_cv([],[])
  self.assertEqual(len(called),4)
  self.assertEqual(called[1][0].count(m.RICH_BUILT),1)
  self.assertNotIn(m.RICH_BUILT,called[1][1])
  for name in m.BETA_BUILT:
   self.assertEqual(called[3][0].count(name),1)
   self.assertNotIn(name,called[3][1])
  for base,added in called:
   self.assertFalse(set(base)&set(added))
 def test_full_world_source_not_a_single_city(self):
  self.assertEqual(len(m.GRID_WIDTHS),2)
  self.assertEqual(m.MIN_CONTINENTS,4)
if __name__=="__main__":unittest.main()
