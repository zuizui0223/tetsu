import unittest,math
from analyze_global_rail_height_road_contrast import road_lengths,fit_roads,ROAD_KM,ROAD_VALID,DENSITY
class HeightRoad(unittest.TestCase):
    def test_common_units(self):
        # 150 m/km2, 40,000 km2 total land, 98% valid => 6000 km vs 5880 km
        r={DENSITY:math.log1p(150)}
        s={"source_GRIP4_total_land_area_m2":4e10,"source_GRIP4_road_coverage_of_source_land":.98,
           "source_GRIP4_observed_road_density_m_per_km2_land":150,
           "source_GRIP4_implied_road_km_in_valid_land":5880}
        out=road_lengths(r,s)
        self.assertAlmostEqual(math.expm1(out[ROAD_KM]),6000)
        self.assertAlmostEqual(math.expm1(out[ROAD_VALID]),5880)
        with self.assertRaises(ValueError):road_lengths(r,dict(s,source_GRIP4_implied_road_km_in_valid_land=123))
    def test_no_support(self):
        self.assertNotEqual(fit_roads([])[0]["status"],"estimated")
if __name__=="__main__":unittest.main()
