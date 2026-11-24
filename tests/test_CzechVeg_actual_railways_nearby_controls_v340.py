"""v340 real-source Czech railway substrate and spatial matching unit checks."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_CzechVeg_actual_railways_nearby_controls_v340 as m

class TrueRailwayPlotCertification(unittest.TestCase):
    def test_actual_railway_track_and_embankment_certified(self):
        self.assertEqual(m.rail_certainty("v kolejisti"),"railway_track_yard_or_structure")
        self.assertEqual(m.rail_certainty("železniční násep"),"railway_embankment_ditch_or_cut")
    def test_adjacent_not_spuriously_on_rail(self):
        self.assertEqual(m.rail_certainty("hromady hlíny cca 200 m od železniční zastávky"),
                         "railway_nearby_or_below_not_confirmed")
        self.assertEqual(m.rail_certainty("mokřad pod železniční tratí"),
                         "railway_nearby_or_below_not_confirmed")
    def test_station_property_not_confirmed_track(self):
        self.assertEqual(m.rail_certainty("nádraží"),
                         "station_property_not_necessarily_track")
    def test_source_geography_and_area(self):
        a=m.src_metadata({"RELEVE_NR":"1","DEG_LAT":"49.5","DEG_LON":"15.2",
                          "SURF_AREA":"4","PRECISION":"30","YEAR":"2001",
                          "REFERENCE":"10"})
        self.assertEqual(a["year"],2001)
        self.assertEqual(a["area"],4)
        self.assertEqual(a["precision"],30)
    def test_missing_year_not_an_imputed_1970(self):
        a=m.src_metadata({"RELEVE_NR":"2","YEAR":"","SURF_AREA":""})
        self.assertIsNone(a["year"])
        self.assertIsNone(a["area"])
    def test_matching_source_year_area_and_precision_gate(self):
        a={"lat":49.0,"lon":15.0,"year":2000,"area":4,"precision":30,
           "reference":"R"}
        b=dict(a,lat=49.001,lon=15.001)
        ok,d=m.in_bounds(a,b,m.MATCH_GATES["G1_same_source_high_precision"])
        self.assertTrue(ok)
        self.assertGreater(d,0)
        c=dict(b,reference="OTHER")
        self.assertFalse(m.in_bounds(a,c,m.MATCH_GATES["G1_same_source_high_precision"])[0])
        self.assertTrue(m.in_bounds(a,c,m.MATCH_GATES["G3_cross_source_high_precision"])[0])
    def test_precision_quality_disqualifies(self):
        a={"lat":49,"lon":15,"year":2000,"area":4,"precision":200,
           "reference":"R"}
        b=dict(a)
        self.assertFalse(m.in_bounds(a,b,m.MATCH_GATES["G1_same_source_high_precision"])[0])
    def test_same_year_not_enough_if_different_plot_area(self):
        a={"lat":49,"lon":15,"year":2000,"area":4,"precision":20,
           "reference":"R"}
        b=dict(a,area=100)
        self.assertFalse(m.in_bounds(a,b,m.MATCH_GATES["G2_same_source_extended"])[0])
    def test_author_duplicate_substudies_not_independent_climate_regions(self):
        self.assertEqual(len(m.MATCH_GATES),4)
        self.assertEqual(len(m.CERTIFIED),2)

if __name__=="__main__":unittest.main()
