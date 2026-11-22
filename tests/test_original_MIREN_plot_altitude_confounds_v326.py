"""v326 synthetic tests for MIREN plot elevation-displacement sensitivity."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_original_MIREN_plot_altitude_confounds_v326 as m
EIVE={"poa annua":{"position":7,"width":3,"systems":3},
      "galium aparine":{"position":5,"width":4,"systems":3}}
def record(region,road,site,alt):
    return {"region":region,"road":road,"site":site,"source_year":2018,
        "plot1_minus_plot3_alt_m":alt,
        "position_delta":2.0,"width_delta":-1.0}
class AltitudeMethods(unittest.TestCase):
    def test_strict_physical_altitude_missing_not_0(self):
        self.assertIsNone(m.num(""))
        self.assertIsNone(m.num("inf"))
        self.assertIsNone(m.num("NA"))
        self.assertEqual(m.num("200.5"),200.5)
    def test_road_region_hierarchy_with_equal_elevation(self):
        rows=[record("AUT",r,str(i),0.) for r in ("r1","r2")
              for i in range(4)]
        panel=m.summarize(rows,"position",0)
        self.assertEqual(panel["supported_regions"],1)
        self.assertEqual(panel["supported_roads"],2)
        self.assertEqual(panel["supported_sites"],8)
        self.assertEqual(panel["region_equal_mean_roadside_minus_interior_EIVE"],2)
    def test_outcome_blind_altitude_gate_reduces_sample(self):
        rows=[record("AUT",r,str(i),(0 if i<3 else 20))
              for r in ("r1","r2") for i in range(4)]
        p=m.summarize(rows,"position",5)
        self.assertEqual(p["all_pre_gate_matched_sites_within_altitude_band"],6)
        self.assertEqual(p["supported_regions"],0)
        self.assertIsNone(p["region_equal_mean_roadside_minus_interior_EIVE"])
    def test_no_imputation_of_site_with_missing_elev(self):
        sel={("AUT","r1"):[{"site":"s","year":2018,
            "roadside":{"poa annua"},"inland":{"galium aparine"}}]}
        out,q=m.make_records(sel,{},EIVE)
        self.assertEqual(out,[])
        self.assertEqual(q["missing_real_paired_source_elevation"],1)
    def test_outcome_sign_separate_from_elevation_sign(self):
        rows=[record("AUT",r,str(i),-5) for r in ("r1","r2")
              for i in range(4)]
        p=m.summarize(rows,"width",10)
        self.assertEqual(p["region_equal_mean_roadside_minus_interior_EIVE"],-1)
    def test_minimum_two_roads_per_region(self):
        rows=[record("NOR","r1",str(i),0) for i in range(4)]
        self.assertEqual(m.summarize(rows,"position",0)["supported_regions"],0)
if __name__=="__main__":unittest.main()
