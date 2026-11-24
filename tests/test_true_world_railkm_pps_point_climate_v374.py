"""v374 global physically sampled railway-km climate vs coarse cell source QA tests."""
import sys, unittest, math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_true_world_railkm_pps_point_climate_v374 as m
from pyproj import Transformer
from shapely.geometry import box
from shapely.prepared import prep

def sample():
    out=[]
    for i in range(m.SAMPLE_K):
        out.append({
            "original_independent_rail_km_point_draw_index":i,
            "source_original_2026_WGS84_track_point_lon":float(10+i%60*.01),
            "source_original_2026_WGS84_track_point_lat":float(40+i%30*.01),
            "source_original_rail_class":"standard_gauge",
            "source_site_geocode_is_synthetic_rail_km_point_NOT_plant_or_observer_GPS":True,
            "source_overpass_or_tunnel_candidate":False
        })
    return {"source_global_rail_sample_meta":{
       "sample_seed_SHA256":m.SEED_SHA,
       "n_independent_uniform_world_rail_km_draws_with_replacement":m.SAMPLE_K,
       "original_128_global_source_parquet_complete_verified":True,
       "sampling_frame_original_fullworld_candidate_normal_rail_metres":m.EXPECTED_CANDIDATE_KM*1000,
       "source_sampled_original_rail_class_counts":{"standard_gauge":m.SAMPLE_K}},
       "original_public_rail_point_draws_no_botanist_coordinates":out}

class AuthenticGlobalLineClimate(unittest.TestCase):
    def test_true_2048_original_draws_pass(self):
        pts,_=m.valid_source_sample(sample())
        self.assertEqual(len(pts),2048)
    def test_truncated_original_sample_rejected(self):
        d=sample();d["original_public_rail_point_draws_no_botanist_coordinates"].pop()
        with self.assertRaises(ValueError):m.valid_source_sample(d)
    def test_wrong_source_salt_rejected(self):
        d=sample();d["source_global_rail_sample_meta"]["sample_seed_SHA256"]="0"*64
        with self.assertRaises(ValueError):m.valid_source_sample(d)
    def test_real_wgs84_100km_200km_grid_origin_consistency(self):
        f=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
        inv=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
        a=m.grid_from_true_point(15,48,100000,f,inv)
        b=m.grid_from_true_point(15,48,200000,f,inv)
        self.assertEqual(a[0][0]//2,b[0][0])
        self.assertEqual(a[0][1]//2,b[0][1])
    def test_macro_temperature_cutoffs_are_exactly_8_and_20(self):
        self.assertEqual(m.safe_band(7.99),"cool")
        self.assertEqual(m.safe_band(8.0),"temperate")
        self.assertEqual(m.safe_band(19.99),"temperate")
        self.assertEqual(m.safe_band(20.0),"warm")
        self.assertIsNone(m.safe_band(None))
    def test_land_geographic_region_not_guessed_if_ocean(self):
        land={c:prep(box(13,47,16,50)) if c=="Europe"
              else prep(box(-180,-90,-170,-85)) for c in m.CONTINENTS}
        self.assertEqual(m.land_source_continent(48,15,land),"Europe")
        self.assertEqual(m.land_source_continent(48,135,land),"unassigned_land_map_units")
    def test_band_fractions_do_not_treat_NA_as_cool(self):
        data=[{"band":"cool"},{"band":"temperate"},{"band":None},{"band":"warm"}]
        r=m.weighted_counts(data,"band")
        self.assertEqual(r["n_valid"],3)
        self.assertAlmostEqual(r["band_fractions_valid"]["cool"],1/3)
    def test_true_origin_mismatch_v351_not_used(self):
        self.assertEqual(m.X_START,-17400000)
        self.assertEqual(m.Y_START,-7400000)
    def test_v373_grid_comparison_requires_both_scales(self):
        d={"global_source_grid_rail_climate_results":[{"source_grid_width_m":100000,
               "source_full_global_climate_context_QA":{
                "original_full_global_rail_geodesic_km":m.EXPECTED_V373_KM},
               "source_each_of_six_continents_climate_band":{c:{} for c in m.CONTINENTS}}]}
        with self.assertRaises(ValueError):m.source_v373(d)
    def test_2048_rail_km_sample_not_plant_presence_data(self):
        self.assertEqual(m.SAMPLE_K,2048)
        self.assertEqual(len(m.BANDS),3)
        self.assertGreater(m.EXPECTED_CANDIDATE_KM,2000000)

if __name__=="__main__":unittest.main()
