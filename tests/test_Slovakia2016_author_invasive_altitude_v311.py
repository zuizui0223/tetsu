"""v311 synthetic tests: public original 14-name author status and additive group null."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
 import fitz
 from analyze_Slovakia2016_invasive_shared_species_altitude_v311 import (
     PUBLISHED_NAMED_INVASIVES,GROUPS,match_original_author_names,dyad_components,
     design,fit,standard_condition,source_elevation_tertiles
 )
except ModuleNotFoundError as e:
 if e.name not in ("fitz","numpy"):raise
 fitz=None

@unittest.skipUnless(fitz is not None,"Source-specific PDF and numerical packages installed only in canonical Actions")
class TestIndependentSourceInvasiveSpatial(unittest.TestCase):
 def test_14_exact_author_explicit_source_names_nonduplicated(self):
  self.assertEqual(len(PUBLISHED_NAMED_INVASIVES),14)
  self.assertEqual(len(set(PUBLISHED_NAMED_INVASIVES)),14)
  self.assertIn("Conyza canadensis",PUBLISHED_NAMED_INVASIVES)
  self.assertNotIn("Erigeron canadensis",PUBLISHED_NAMED_INVASIVES)
 def test_strict_name_match_does_not_guess_fifteenth_or_synonym(self):
  obs={k:{1} for k in PUBLISHED_NAMED_INVASIVES}
  m=match_original_author_names(obs)
  self.assertEqual(len(m),14)
  broken=dict(obs);broken.pop("Conyza canadensis")
  with self.assertRaises(ValueError):match_original_author_names(broken)
 def test_original_all_shared_taxon_denominator_used_by_invasive(self):
  source={1:{"A","B","C"},2:{"A","C"},3:{"B","C"}}
  pairs=[(1,2),(1,3),(2,3)]
  c=dyad_components(source,{"A","B"},pairs)
  self.assertEqual(c.shape,(3,3))
  np.testing.assert_allclose(c[:,0],c[:,1]+c[:,2],atol=1e-12)
  self.assertAlmostEqual(float(c[0,1]),1/3)
  self.assertAlmostEqual(float(c[1,1]),1/3)
  self.assertAlmostEqual(float(c[2,1]),0.0)
 def test_joint_regression_coefficient_additivity(self):
  km=np.array([1.,2.,4.,7.,10.,13.])
  alt=np.array([.1,1.,2.,.25,.7,3.])
  X=design(km,alt)
  invasive=2-0.1*km+0.04*alt
  other=1-0.05*km+0.08*alt
  C=np.column_stack([invasive+other,invasive,other])
  b=fit(X,C)
  self.assertAlmostEqual(float(b[1,0]),float(b[1,1]+b[1,2]),places=10)
  cond,corr=standard_condition(km,alt)
  self.assertTrue(np.isfinite(cond))
  self.assertTrue(np.isfinite(corr))
 def test_elevation_tertiles_are_fixed_by_altitude_not_plant_results(self):
  c={i:(49.,20.,700-10*i,str(i)) for i in range(1,20)}
  g=source_elevation_tertiles(c)
  self.assertEqual(list(map(len,g)),[6,6,7])
  self.assertEqual(g[0],[1,2,3,4,5,6])
  self.assertEqual(g[1],[7,8,9,10,11,12])
 def test_source_error_not_silently_zeroed(self):
  pairs=[(1,2)]
  with self.assertRaises(ValueError):
   dyad_components({1:set(),2:set()},set(),pairs)
if __name__=="__main__":unittest.main()
