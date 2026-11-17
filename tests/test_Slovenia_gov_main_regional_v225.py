import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from analyze_Slovenia_gov_official_main_regional_traits_v225 import (
  MAIN,REGIONAL,KNOWN,official_class,categorical_contrast,mean_by_section,aggregate
)
class TestGovernmentLineCodeMapping(unittest.TestCase):
    def test_all_15_source_codes_are_distinct_and_classified(self):
        self.assertEqual(len(KNOWN),15)
        self.assertEqual(len(set(KNOWN)),15)
        self.assertEqual({official_class(x) for x in KNOWN},{"main_G","regional_R"})
        self.assertEqual(sum(official_class(x)=="main_G" for x in KNOWN),6)
        self.assertEqual(sum(official_class(x)=="regional_R" for x in KNOWN),9)
    def test_historically_verified_G44_R70_not_line_speed_proxy(self):
        self.assertEqual(official_class("44"),"main_G")
        self.assertEqual(official_class("70"),"regional_R")
        self.assertEqual(official_class(64),"main_G")
        self.assertIsNone(official_class("unverified"))
    def test_section_equal_metric_not_site_weighted(self):
        data=[{"section":"10","m":1.0}]*5+[{"section":"20","m":3.0},
              {"section":"21","m":0.0},{"section":"31","m":1.0}]
        out=mean_by_section(data,"m")
        self.assertEqual(out["10"],1.)
        self.assertEqual(out["20"],3.)
        self.assertEqual(out["31"],1.)
        self.assertEqual(categorical_contrast(out)["main_minus_regional"],1.5)
        self.assertEqual(aggregate([1,2,3])["median"],2)
    def test_one_year_not_assumed_due_to_categorization(self):
        self.assertEqual(official_class("21"),"regional_R")
        self.assertEqual(official_class("41"),"main_G")
if __name__=="__main__":unittest.main()
