"""v376b exact worldwide original GHSL 2020 builtup m2 coarse cell rollup tests."""
import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import aggregate_original_GHSL2020_global_builtup_v376b as m

class GHSL2020Area(unittest.TestCase):
 def test_same_global_GERS_plant_grid_geographic_origin(self):
  self.assertEqual((m.X_START,m.Y_START),(-17400000,-7400000))
  self.assertEqual(m.ROWS,148)
  self.assertEqual(m.COLS,348)
 def test_source_JRC_2020_globe_mollweide_zip(self):
  self.assertIn("_E2020_",m.SOURCE_URL)
  self.assertIn("_54009_1000",m.SOURCE_URL)
  self.assertEqual(m.EXPECTED_ZIP_BYTES,152527564)
 def test_preserve_builtup_m2_under_100km_200km_sum(self):
  data=np.full((m.ROWS,m.COLS),m.NODATA)
  ok=np.zeros((m.ROWS,m.COLS),dtype=bool)
  # two north-origin 100km rows, ensure world source 200km cell has
  # four 100km subcells with exact m2 sum and source southern iy.
  yy=m.ROWS-2
  data[yy:yy+2,:2]=np.array([[1e8,2e8],[3e8,4e8]])
  ok[yy:yy+2,:2]=True
  a=m.rowwise_coarse_grid(data,ok,100000)
  b=m.rowwise_coarse_grid(data,ok,200000)
  self.assertEqual(len(a),4)
  self.assertEqual(len(b),1)
  self.assertAlmostEqual(b[0]["source_GHSL2020_total_builtup_m2_in_cell"],1e9)
  self.assertAlmostEqual(b[0]["source_GHSL2020_builtup_fraction_in_original_cell"],1e9/(200000**2))
  self.assertEqual(sorted(x["source_equal_area_grid_iy"] for x in a),[0,0,1,1])
 def test_builtup_zero_retains_truly_valid_zero(self):
  data=np.full((m.ROWS,m.COLS),m.NODATA)
  ok=np.zeros((m.ROWS,m.COLS),dtype=bool)
  data[0,0]=0.;ok[0,0]=True
  a=m.rowwise_coarse_grid(data,ok,100000)
  self.assertEqual(len(a),1)
  self.assertEqual(a[0]["source_GHSL2020_total_builtup_m2_in_cell"],0.)
 def test_unknown_ocean_not_invented_zero(self):
  data=np.full((m.ROWS,m.COLS),m.NODATA)
  ok=np.zeros((m.ROWS,m.COLS),dtype=bool)
  self.assertEqual(m.rowwise_coarse_grid(data,ok,100000),[])
 def test_ghsl_pixel_builtup_area_cannot_exceed_cell(self):
  data=np.full((m.ROWS,m.COLS),m.NODATA)
  ok=np.zeros((m.ROWS,m.COLS),dtype=bool)
  data[0,0]=1.1e10;ok[0,0]=True
  with self.assertRaises(ValueError):
   m.rowwise_coarse_grid(data,ok,100000)
 def test_original_2020_full_source_nodata_is_separate_from_valid_urban_zero(self):
  self.assertLess(m.NODATA,0)
  self.assertEqual(m.WIDTH,100000)
 def test_100km_and_200km_are_two_resolution_of_same_world(self):
  self.assertEqual(m.BANDS,(100000,200000))
  self.assertEqual(m.ROWS%2,0)
  self.assertEqual(m.COLS%2,0)
if __name__=="__main__":unittest.main()
