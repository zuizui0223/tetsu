"""v336 tests for MIREN source ordinal cover-class weighting."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_EIVE_source_cover_class_v336 as m

EIVE={
 "poa annua":{"position":7.0,"width":4.0,"systems":3},
 "galium aparine":{"position":5.0,"width":3.0,"systems":3},
 "lotus corniculatus":{"position":6.0,"width":2.0,"systems":3},
 "lactuca serriola":{"position":8.0,"width":5.0,"systems":3},
 "sonchus arvensis":{"position":4.0,"width":6.0,"systems":3},
}
def fixture():
 roads={}
 cover={}
 for road in ("r1","r2"):
  sites=[]
  for i in range(4):
   species=set(EIVE)
   site=str(i)
   sites.append({"site":site,"year":2018,
             "roadside":species,"inland":species})
   cover["AUT",road,site,2018,"1"]={x:(8 if x=="lactuca serriola" else 1) for x in species}
   cover["AUT",road,site,2018,"3"]={x:(8 if x=="sonchus arvensis" else 1) for x in species}
  roads["AUT",road]=sites
 return roads,cover
class MIRENOrdinalMethods(unittest.TestCase):
 def test_source_cover_missing_not_biological_zero(self):
  self.assertIsNone(m.numeric_class(""))
  self.assertIsNone(m.numeric_class("NA"))
 def test_cover_class_outside_published_1_to_8_rejected(self):
  for v in ("0","9","inf","bogus"):
   with self.assertRaises(ValueError):
    m.numeric_class(v)
 def test_presence_vs_cover_class_same_recorded_species_can_diverge(self):
  sp=set(EIVE)
  vals={s:(8 if s=="lactuca serriola" else 1) for s in sp}
  a=m.cwm(sp,vals,EIVE,"position","presence")
  b=m.cwm(sp,vals,EIVE,"position","cover_class")
  self.assertGreater(b["mean"],a["mean"])
  self.assertEqual(a["covered_trait_species"],5)
  self.assertEqual(a["cover_completeness"],1)
 def test_cover_class_sqrt_between_presence_and_full(self):
  sp=set(EIVE)
  vals={s:(8 if s=="lactuca serriola" else 1) for s in sp}
  a=m.cwm(sp,vals,EIVE,"position","presence")
  b=m.cwm(sp,vals,EIVE,"position","sqrt_cover_class")
  c=m.cwm(sp,vals,EIVE,"position","cover_class")
  self.assertTrue(a["mean"] < b["mean"] < c["mean"])
 def test_cover_missing_changes_coverage_denominator(self):
  sp=set(EIVE)
  data={s:1 for s in list(EIVE)[:3]}
  a=m.cwm(sp,data,EIVE,"position","presence")
  self.assertAlmostEqual(a["cover_completeness"],.6)
  self.assertEqual(a["covered_trait_species"],3)
 def test_source_road_region_minimum_support(self):
  roads,cover=fixture()
  cohorts,qa=m.collect_sitepanels(roads,cover,EIVE,.8)
  self.assertEqual(len(cohorts),2)
  self.assertEqual(sum(len(x) for x in cohorts.values()),8)
  panel=m.mean_by_region(cohorts,EIVE,"position","cover_class")
  self.assertEqual(panel["n_source_regions"],1)
  self.assertEqual(panel["n_source_roads"],2)
  self.assertEqual(panel["n_original_source_site_year_first_visits"],8)
  self.assertGreater(panel["region_equal_mean_roadside_minus_inland_EIVE"],0)
 def test_cover_known_strict_gate_blocks_incomplete_plot(self):
  roads,cover=fixture()
  for x in cover:
   if x[-1]=="1":
    cover[x].pop("poa annua")
  cohorts,qa=m.collect_sitepanels(roads,cover,EIVE,1.0)
  self.assertEqual(len(cohorts),0)
  self.assertEqual(qa["dropped_cover_missingness_or_min_species"],8)
 def test_identical_source_species_two_zones_presence_unchanged(self):
  roads,cover=fixture()
  sel,_=m.collect_sitepanels(roads,cover,EIVE,.8)
  results=m.mean_by_region(sel,EIVE,"position","presence")
  self.assertAlmostEqual(results["region_equal_mean_roadside_minus_inland_EIVE"],0)
if __name__=="__main__":unittest.main()
