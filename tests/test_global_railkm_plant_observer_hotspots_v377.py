"""v377 high-observer-count 2025 iNaturalist source cell exclusion sensitivity tests."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_global_railkm_plant_observer_hotspots_v377 as m

def row(i,cont="Europe",n=20,rail=0):
 return {"continent":cont,"key":(i,0),"observers":n,
         "observations":100,"species":40,"rail_km":rail,
         "log_rail_km":math.log1p(rail)}
class GlobalObserverHotspotAudit(unittest.TestCase):
 def test_zero_fraction_retains_original_source_2025_cell_population(self):
  original=[row(i,n=10+i) for i in range(25)]
  kept,qa=m.select_unhot_cells(original,0.)
  self.assertEqual(len(kept),len(original))
  self.assertEqual(qa["Europe"]["n_excluded_source_public_plant_observer_hotspot_cells"],0)
 def test_top_five_percent_removed_by_observer_effort_only(self):
  vals=[row(i,n=i+5) for i in range(100)]
  vals[-1]["rail_km"]=10000
  vals[0]["rail_km"]=1000000
  kept,qa=m.select_unhot_cells(vals,.05)
  self.assertEqual(len(kept),95)
  self.assertEqual(qa["Europe"]["n_excluded_source_public_plant_observer_hotspot_cells"],5)
  self.assertNotIn((99,0),{r["key"] for r in kept})
  self.assertIn((0,0),{r["key"] for r in kept})
 def test_continent_separate_source_sampling_no_global_only_threshold(self):
  original=[row(i,cont="Europe",n=i+5) for i in range(100)]
  # Real source globally unique EPSG6933 cell ids cannot be reused
  # across two continents even in synthetic input fixtures.
  original += [row(i+100,cont="Asia",n=i+15) for i in range(20)]
  kept,qa=m.select_unhot_cells(original,.10)
  self.assertEqual(qa["Europe"]["n_excluded_source_public_plant_observer_hotspot_cells"],10)
  self.assertEqual(qa["Asia"]["n_excluded_source_public_plant_observer_hotspot_cells"],2)
  self.assertEqual(len(kept),108)
 def test_tied_observer_count_exclusion_deterministic_original_grid_key(self):
  original=[row(3,n=100),row(2,n=100),row(1,n=100),
            row(0,n=20)] + [row(i,n=10) for i in range(4,20)]
  kept,_=m.select_unhot_cells(original,.10)
  self.assertNotIn((1,0),{r["key"] for r in kept})
  self.assertNotIn((2,0),{r["key"] for r in kept})
  self.assertIn((3,0),{r["key"] for r in kept})
 def test_wrong_posthoc_threshold_fails_closed(self):
  with self.assertRaises(ValueError):m.select_unhot_cells([row(1)],.25)
 def test_source_continent_reference_not_inferred_from_rail_density(self):
  with self.assertRaises(ValueError):m.select_unhot_cells([row(1,cont="Atlantis")],.05)
 def test_original_frozen_source_regression_exact(self):
  self.assertEqual(m.FRACTIONS,(0.,.05,.10))
  self.assertAlmostEqual(m.CORRECTED_V371[100000]["richness"],.004463745492835951)
  self.assertAlmostEqual(m.CORRECTED_V371[200000]["jaccard"],-.00006120015031440586)

if __name__=="__main__":unittest.main()
