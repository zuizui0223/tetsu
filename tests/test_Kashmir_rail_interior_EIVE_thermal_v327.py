"""v327 synthetic fixed-gate Kashmir railway vs same-site interior tests."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_rail_interior_EIVE_thermal_v327 as m

EIVE={
 "poa annua":{"position":7.0,"width":4.0,"systems":4},
 "galium aparine":{"position":5.0,"width":3.0,"systems":3},
 "lactuca serriola":{"position":6.0,"width":3.5,"systems":1}
}
NAMES=["poa annua","galium aparine","lactuca serriola","unknown"]
def fake_sites():
    sites=[f"Station {i}" for i in range(1,17)]+[f"Between {i}" for i in range(1,16)]
    plots={}
    for s in sites:
        for year in (2014,2017):
            plots[s,year,"Track"]={0,2}
            plots[s,year,"Interior Plot 1"]={1,2}
            plots[s,year,"Interior Plot 2"]={1,3}
            plots[s,year,"Embankment"]={0,1}
    return plots

class KashmiriThermal(unittest.TestCase):
    def test_missing_plot_fails_closed(self):
        self.assertIsNone(m.profile(None,NAMES,EIVE,"position",3))
        delta,qa=m.site_effect(fake_sites(),NAMES,EIVE,"Station 1",2014,
                              ("Track","Interior Plot 2"),"position",3,2,.5)
        self.assertIsNone(delta)
        self.assertEqual(qa["reason"],"trait_coverage_gate")
    def test_exact_matched_EIVE_index_denominator(self):
        o=m.profile({0,2,3},NAMES,EIVE,"position",3)
        self.assertEqual(o["n_source"],3)
        self.assertEqual(o["n_eive"],1)
        self.assertAlmostEqual(o["coverage"],1/3)
        self.assertEqual(o["cwm"],7)
    def test_min_source_systems_matters(self):
        one=m.profile({0,2,3},NAMES,EIVE,"position",1)
        three=m.profile({0,2,3},NAMES,EIVE,"position",3)
        self.assertEqual(one["n_eive"],2)
        self.assertEqual(three["n_eive"],1)
    def test_station_between_and_years_do_not_double_sample(self):
        plots=fake_sites()
        p=m.panel(plots,NAMES,EIVE,2014,("Track","Interior Plot 2"),
                  "position",1,1,.25)
        self.assertEqual(p["observed"]["all"]["sites"],31)
        self.assertEqual(p["observed"]["station"]["sites"],16)
        self.assertEqual(p["observed"]["between"]["sites"],15)
        self.assertEqual(p["observed"]["all"]["track_higher_sites"],31)
        self.assertEqual(p["observed"]["all"]["mean_track_minus_interior"],1.5)
    def test_missing_xlsx_source_site_not_a_false_zero(self):
        plots=fake_sites()
        plots["Station 1",2017,"Track"]=None
        p=m.panel(plots,NAMES,EIVE,2017,("Track","Interior Plot 2"),
                  "position",1,1,.25)
        self.assertEqual(p["observed"]["all"]["sites"],30)
        self.assertEqual(p["excluded_source_sites_by_reason"]["missing_source_record"],1)
    def test_repeated_year_temporal_is_site_paired(self):
        plots=fake_sites()
        a=m.panel(plots,NAMES,EIVE,2014,("Track","Interior Plot 2"),"position",1,1,.25)
        plots["Station 1",2017,"Track"]=None
        b=m.panel(plots,NAMES,EIVE,2017,("Track","Interior Plot 2"),"position",1,1,.25)
        t=m.temporal(a,b)
        self.assertEqual(t["n_common_spatial_sites"],30)
        self.assertEqual(t["within_site_temporal"]["all"]["mean_2017_minus_2014_of_track_interior_delta"],0)
    def test_strict_gate_marks_nonestimable_not_positive(self):
        p=m.panel(fake_sites(),NAMES,EIVE,2014,("Track","Interior Plot 2"),"position",3,5,.5)
        self.assertEqual(p["observed"]["all"]["sites"],0)
        self.assertIsNone(p["observed"]["all"]["mean_track_minus_interior"])
    def test_temperature_width_not_interchangeable(self):
        a=m.profile({0,1},NAMES,EIVE,"position",3)
        b=m.profile({0,1},NAMES,EIVE,"width",3)
        self.assertNotEqual(a["cwm"],b["cwm"])

if __name__=="__main__":unittest.main()
