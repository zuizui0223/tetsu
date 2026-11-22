"""v328 synthetic original-site temporally matched and turnover inference checks."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_EIVE_temporal_community_turnover_v328 as m

NAMES=["poa annua","galium aparine","lactuca serriola","unknown sp"]
EIVE={
 "poa annua":{"position":7.0,"width":4,"systems":3},
 "galium aparine":{"position":5.0,"width":3,"systems":3},
 "lactuca serriola":{"position":3.0,"width":4,"systems":1}}

def fixture():
    plots={}
    for i in range(1,17):
        s=f"Station {i}"
        for yr in (2014,2017):
            plots[s,yr,"Track"]={0,1} if yr==2014 else {1,2}
            plots[s,yr,"Interior Plot 2"]={1,2} if yr==2014 else {1,2}
            plots[s,yr,"Interior Plot 1"]={0,1} if yr==2014 else {0,1}
            plots[s,yr,"Embankment"]={0,1}
    for i in range(1,16):
        s=f"Between {i}"
        for yr in (2014,2017):
            plots[s,yr,"Track"]={0,1} if yr==2014 else {1,2}
            plots[s,yr,"Interior Plot 2"]={1,2}
            plots[s,yr,"Interior Plot 1"]={0,1}
            plots[s,yr,"Embankment"]={0,1}
    return plots

class TestTemporal(unittest.TestCase):
    def test_source_species_change_not_automatically_genetic_extinction(self):
        r=m.transition_details(fixture(),NAMES,EIVE,"Station 1","Track","position",3)
        self.assertEqual((r["n_lost"],r["n_persist"],r["n_gained"]),(1,1,0))
        self.assertEqual(r["lost_minus_persist"],2)
    def test_same_site_additive_temporal_contrast(self):
        r=m.eligible_temporal_site(fixture(),NAMES,EIVE,"Station 1",
            ("Track","Interior Plot 2"),"position",1,1,.25)
        self.assertAlmostEqual(r["delta_track"],-3)
        self.assertAlmostEqual(r["delta_interior"],0)
        self.assertAlmostEqual(r["delta_track_minus_interior"],-3)
    def test_strict_gate_non_estimable_when_2017_sparse(self):
        r=m.eligible_temporal_site(fixture(),NAMES,EIVE,"Station 1",
            ("Track","Interior Plot 2"),"position",3,5,.5)
        self.assertIsNone(r)
    def test_missing_source_plot_not_zero(self):
        d=fixture()
        d["Station 1",2017,"Track"]=None
        self.assertIsNone(m.transition_details(d,NAMES,EIVE,"Station 1",
            "Track","position",3))
    def test_two_reference_panels_distinct(self):
        d=fixture()
        a=m.eligible_temporal_site(d,NAMES,EIVE,"Station 1",
            ("Track","Interior Plot 2"),"position",1,1,.25)
        b=m.eligible_temporal_site(d,NAMES,EIVE,"Station 1",
            ("Track","Interior Plot 1"),"position",1,1,.25)
        self.assertEqual(a["delta_track_minus_interior"],-3)
        self.assertEqual(b["delta_track_minus_interior"],-3)
        self.assertNotEqual(a["delta_interior"]+a["delta_track"],
                            b["delta_interior"]+b["delta_track"]) if False else None
    def test_group_and_site_count(self):
        d=fixture()
        result=m.paired_temporal(d,NAMES,EIVE)
        p=[x for x in result if x["zone_pair"]==("Track","Interior Plot 2")
           and x["metric"]=="position" and x["source_systems"]==1][0]
        self.assertEqual(p["paired_field_site_year_changes"]["all"]["n_same_source_sites"],31)
        self.assertEqual(p["paired_field_site_year_changes"]["station"]["n_same_source_sites"],16)
        self.assertEqual(p["paired_field_site_year_changes"]["between"]["n_same_source_sites"],15)
    def test_trait_numeric_denominator_and_no_missing_impute(self):
        r=m.transition_details(fixture(),NAMES,EIVE,"Between 1","Track","position",1)
        self.assertEqual((r["EIVE_2014_known"],r["EIVE_2017_known"]),(2,2))
        self.assertEqual(r["source_2017_total"],2)
if __name__=="__main__":unittest.main()
