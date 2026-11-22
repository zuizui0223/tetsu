"""v334 original source numerical weighting synthetic tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_cover_weighted_matched_temperature_v334 as m

NAMES=["poa annua","galium aparine","lotus corniculatus"]
EIVE={"poa annua":{"position":7.0,"width":2.0,"systems":3},
      "galium aparine":{"position":5.0,"width":4.0,"systems":3},
      "lotus corniculatus":{"position":6.0,"width":3.0,"systems":3}}
ORIG={
  ("Between 1",2014,"Track"):{0:10,1:1},
  ("Between 1",2014,"Interior Plot 2"):{1:10,2:1},
  ("Station 1",2014,"Track"):{0:1,1:10},
  ("Station 1",2014,"Interior Plot 2"):{1:9,2:1}
}
PAIRS=[{"between_source_site":"Between 1","station_source_site":"Station 1"}]
class SourceCover(unittest.TestCase):
 def test_presence_cwm_unweighted_same_two_taxa(self):
  x=m.weighted_community({0:100,1:1},NAMES,EIVE,3,"position","presence")
  self.assertEqual(x["cwm"],6)
  self.assertEqual(x["effective_trait_species_number"],2)
 def test_original_positive_source_weight_is_not_percentage(self):
  x=m.weighted_community({0:100,1:1},NAMES,EIVE,3,"position","original_positive_source_value")
  self.assertGreater(x["cwm"],6.9)
  self.assertGreater(x["largest_single_species_weight_fraction"],.99-1e-3)
 def test_sqrt_weight_less_extreme_than_raw(self):
  a=m.weighted_community({0:100,1:1},NAMES,EIVE,3,"position","original_positive_source_value")
  b=m.weighted_community({0:100,1:1},NAMES,EIVE,3,"position","sqrt_original_source_value")
  self.assertGreater(a["cwm"],b["cwm"])
  self.assertGreater(b["cwm"],6)
 def test_unknown_weight_method_fails_closed(self):
  with self.assertRaises(ValueError):
   m.weighted_community({0:2},NAMES,EIVE,3,"position","linear_imputed_temperature")
 def test_missing_species_observations_not_zero(self):
  self.assertIsNone(m.weighted_community(None,NAMES,EIVE,3,"position","presence"))
 def test_four_source_plot_comparison_preserves_additivity(self):
  one=m.group_comparison(ORIG,PAIRS,NAMES,EIVE,2014,
      "strict_3systems_5species_50pct","presence","position")
  two=m.group_comparison(ORIG,PAIRS,NAMES,EIVE,2014,
      "strict_3systems_5species_50pct","original_positive_source_value","position")
  self.assertAlmostEqual(one["mean_between_minus_station_track_relative_interior"],0)
  self.assertGreater(two["mean_between_minus_station_track_relative_interior"],0)
  self.assertAlmostEqual(
   two["mean_between_minus_station_track_relative_interior"],
   two["mean_between_minus_station_track_only"]-
   two["mean_between_minus_station_interior_only"])
 def test_no_independent_signflip_claim_when_one_pair(self):
  x=m.group_comparison(ORIG,PAIRS,NAMES,EIVE,2014,
       "strict_3systems_5species_50pct","presence","position")
  self.assertIsNone(x["conditional_matched_signflip_diagnostic"])
 def test_fixed_panels_and_weight_methods(self):
  self.assertEqual(len(m.PANELS),4)
  self.assertEqual(len(m.WEIGHTINGS),3)
if __name__=="__main__":unittest.main()
