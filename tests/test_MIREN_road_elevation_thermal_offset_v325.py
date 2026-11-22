"""v325 synthetic elevation gradient tests: no original MIREN plot rows here."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_road_elevation_thermal_offset_v325 as m

T={"poa annua":{"position":7.,"width":4.,"systems":3},
   "galium aparine":{"position":5.,"width":3.,"systems":3}}
def row(road,site,elev):
    return {"site":site,"year":2018,"elevation_m":elev,
            "roadside":{"poa annua"},"inland":{"galium aparine"}}

class SlopeTests(unittest.TestCase):
    def test_positive_slope_per_100m(self):
        self.assertAlmostEqual(m.road_slope([100,200,300,400],[0,1,2,3]),1.0)
    def test_negative_slope_per_100m(self):
        self.assertAlmostEqual(m.road_slope([100,200,300,400],[3,2,1,0]),-1.0)
    def test_no_road_elevation_variation(self):
        self.assertIsNone(m.road_slope([100]*4,[0,1,2,3]))
    def test_bad_elevation_missing_not_zero(self):
        self.assertIsNone(m.float_elev(""))
        self.assertIsNone(m.float_elev("NA"))
        self.assertIsNone(m.float_elev("inf"))
        self.assertEqual(m.float_elev("1200.5"),1200.5)
    def test_elevation_site_pair_requires_actual_source_year(self):
        selected={("AUT","r1"):[{"site":"a","year":2018,"roadside":{"poa annua"},
                    "inland":{"galium aparine"}} for _ in range(4)]}
        result,qa=m.build_cohort(selected,{})
        self.assertEqual(result,{})
        self.assertEqual(qa["selected_sites_without_same_year_plot_elevation"],4)
    def test_region_not_inferred_from_one_road(self):
        sel={("AUT","r1"):[row("r1",str(i),100+100*i) for i in range(4)]}
        d=m.model(sel,T,"position")
        self.assertEqual(d["eligible_original_region_codes"],0)
        self.assertIsNone(d["region_equal_mean_delta_EIVE_per_100m_altitude"])
    def test_two_road_four_site_region_unit(self):
        sel={("AUT",r):[row(r,str(i),100+100*i) for i in range(4)]
             for r in ("r1","r2")}
        d=m.model(sel,T,"position")
        self.assertEqual(d["eligible_original_region_codes"],1)
        self.assertEqual(d["eligible_original_road_ids"],2)
        self.assertEqual(d["eligible_site_count"],8)
        self.assertAlmostEqual(d["region_equal_mean_delta_EIVE_per_100m_altitude"],0)
if __name__=="__main__":unittest.main()
