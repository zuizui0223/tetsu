"""v335 synthetic source numeric-weight gradient and taxon omission tests."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Kashmir_original_weight_sign_reversal_v335 as m
TRAIT={"hot":{"position":8.},"cold":{"position":4.},"mid":{"position":6.}}
PAIR={"BT":{"hot":10.,"cold":1.},"BI":{"cold":4.,"mid":1.},
      "ST":{"hot":1.,"cold":10.},"SI":{"cold":4.,"mid":1.}}
class SourceWeightControl(unittest.TestCase):
 def test_presence_cwm_diff_zero_and_source_weight_diff_positive(self):
  p=m.weighted_panel([PAIR],TRAIT,"position",0.,None)
  raw=m.weighted_panel([PAIR],TRAIT,"position",1.,None)
  self.assertAlmostEqual(p["mean_between_minus_station_track_relative_interior"],0)
  self.assertGreater(raw["mean_between_minus_station_track_relative_interior"],3)
 def test_sqrt_power_intermediate_for_synthetic_dominant_hot(self):
  a=m.weighted_panel([PAIR],TRAIT,"position",0.0,None)
  b=m.weighted_panel([PAIR],TRAIT,"position",0.5,None)
  c=m.weighted_panel([PAIR],TRAIT,"position",1.0,None)
  self.assertLess(a["mean_between_minus_station_track_relative_interior"],
                  b["mean_between_minus_station_track_relative_interior"])
  self.assertLess(b["mean_between_minus_station_track_relative_interior"],
                  c["mean_between_minus_station_track_relative_interior"])
 def test_original_positive_source_cap_limits_extreme(self):
  p=dict(PAIR, BT={"hot":103.6,"cold":1.})
  original=m.weighted_panel([p],TRAIT,"position",1,None)
  cap=m.weighted_panel([p],TRAIT,"position",1,10)
  self.assertGreater(original["mean_between_minus_station_track_relative_interior"],
                     cap["mean_between_minus_station_track_relative_interior"])
 def test_count_above100_no_source_site_ids(self):
  p=dict(PAIR,BT={"hot":103.6,"cold":1})
  out=m.panel_source_outlier_trace([p])
  self.assertEqual(out["BT"],1)
  self.assertEqual(out["BI"],0)
 def test_single_taxon_omission_tracks_eligible_pairs(self):
  r=m.single_taxon_impact([PAIR],TRAIT)
  self.assertEqual(r["n_source_EIVE_taxa_in_source_matched_plots"],3)
  self.assertEqual(r["n_evaluable_one_taxon_omissions"],3)
  self.assertEqual(r["n_preserving_all_paired_geographic_sites"],3)
 def test_invalid_empty_source_record_returns_not_estimable(self):
  self.assertIsNone(m.community_mean({},TRAIT,"position",1,None))
 def test_source_negative_number_fails_closed(self):
  with self.assertRaises(ValueError):m.community_mean({"hot":-1.},TRAIT,"position",1,None)
 def test_predefined_panels_and_exponents(self):
  self.assertEqual(len(m.PANELS),4)
  self.assertEqual(m.EXPONENTS,(0.,.25,.5,.75,1.))
  self.assertEqual(m.CAPS_AT_ONE,(None,100,20,10,5))
if __name__=="__main__":unittest.main()
