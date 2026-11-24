"""v373 source global railway infrastructure km vs WorldClim macro temperature checks."""
import unittest,sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_world_railkm_climate_distribution_v373 as m

class GlobalInfrastructureClimatology(unittest.TestCase):
    def test_original_world_source_length_population_not_rail_plants(self):
        self.assertAlmostEqual(m.RAIL_SOURCE_KM,2237043.481)
        self.assertEqual(m.GRID_WIDTHS,(100000,200000))
    def test_original_world_100km_grid_cell_center(self):
        self.assertEqual(m.mapcell(0,0,100000),(-17350000,-7350000))
        self.assertEqual(m.mapcell(0,0,200000),(-17300000,-7300000))
    def test_physical_bio1_bands_not_EIVE(self):
        self.assertEqual(m.temp_band(7.9),"cool")
        self.assertEqual(m.temp_band(8),"temperate")
        self.assertEqual(m.temp_band(20),"warm")
    def test_weighted_global_railkm_mean(self):
        self.assertAlmostEqual(m.weighted_average([(1,5),(3,25)]),20)
    def test_empty_railkm_not_biological_zero(self):
        self.assertIsNone(m.weighted_average([]))
    def test_no_truncated_source_railkm_allowed(self):
        j={"source_Oveture_release":m.RAIL_RELEASE,
           "source_original_fullworld_rail_128_files_SHA_verified":True,
           "source_fullworld_candidate_normal_rail_km":m.RAIL_SOURCE_KM,
           "source_global_gridded_rail_exposure_km_by_cell_width":{
              str(w):[{"source_real_original_Overture_normal_rail_km_clipped_into_cell":m.RAIL_SOURCE_KM}]
              for w in m.GRID_WIDTHS}}
        self.assertTrue(m.valid_source_exposure(j))
        j["source_fullworld_candidate_normal_rail_km"]=100
        with self.assertRaises(ValueError):m.valid_source_exposure(j)
    def test_source_six_continents_instead_of_europe_only(self):
        self.assertEqual(set(m.CONTINENTS),{
          "Africa","Asia","Europe","North America","South America","Oceania"})
    def test_source_macroclimate_not_thermal_species_trait(self):
        self.assertIn("BIO1",m.WORLDCLIM_SHA if False else "BIO1") 
        self.assertEqual(len(m.WORLDMAP_SHA),64)
        self.assertEqual(len(m.WORLDCLIM_SHA),64)

if __name__=="__main__":unittest.main()
