import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from audit_Parma_rail_vs_Italy_MIREN_road_geographic_mismatch_v252 import (
  CITY_PARMA_UNGEGN,git_blob,haversine_km,min_to_lonlat_bbox,analyze
)
class TestTrueParmaVsDistantMIRENItalianRoadReference(unittest.TestCase):
    def source(self):
        return {
         "sites":{"IT":28},
         "controls":{"IT_road":{"roads":3,"plots":60,"years":[2022],
                     "lat":[42.149692,42.47389002],"lon":[12.89728797,14.13079]}},
         "comparisons":{"IT_road":{"common":{"delta":.30809816690971126,
                                  "genus_interval":[.1380,.4739]},
                                  "extended":{"delta":.12511153764367122}}}}
    def test_real_city_to_miren_bbox_distance(self):
        source=self.source()["controls"]["IT_road"]
        val,corner=min_to_lonlat_bbox(CITY_PARMA_UNGEGN,source["lat"],source["lon"])
        self.assertGreater(val,330)
        self.assertLess(val,333)
        self.assertAlmostEqual(corner[0],42.47389002)
        self.assertAlmostEqual(corner[1],12.89728797)
    def test_city_inside_bbox_is_zero_distance(self):
        d,_=min_to_lonlat_bbox((44.8,10.3),[44,45],[9,11])
        self.assertAlmostEqual(d,0)
    def test_reject_invalid_bbox(self):
        with self.assertRaises(ValueError):min_to_lonlat_bbox((44.8,10.3),[43,42],[9,11])
    def test_study_era_gap_and_no_all_station_geolocation_claim(self):
        a=analyze(self.source())
        self.assertEqual(a['source_original_rail_station_28_from_Parma_PROVINCE'],28)
        self.assertEqual(a['domestic_reference_MIREN_mountain_road_source'][
            'source_2022_minus_latest_original_main_rail_year'],9)
        self.assertTrue(a['not_claim_each_of_the_28_original_Parma_railway_stations_is_331km_from_road_reference'])
    def test_geography_source_drift_halts(self):
        s=self.source()
        s["controls"]["IT_road"]["plots"]=61
        with self.assertRaises(ValueError):analyze(s)
if __name__=="__main__":unittest.main()
