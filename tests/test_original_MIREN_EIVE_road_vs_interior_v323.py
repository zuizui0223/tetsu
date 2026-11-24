"""v323 synthetic checks for original MIREN pairs and EIVE temperature traits."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_original_MIREN_EIVE_road_vs_interior_v323 as m

EIVE={
  "poa annua":{"position":7.,"width":3.,"systems":3},
  "galium aparine":{"position":5.,"width":4.,"systems":3},
  "lactuca serriola":{"position":6.,"width":3.,"systems":1},
}
class RoadBenchmarkMethods(unittest.TestCase):
    def test_earliest_complete_year_not_latest(self):
        data={
         ("AUT","R","S",2011,"1"):{"poa annua"},
         ("AUT","R","S",2011,"3"):{"galium aparine"},
         ("AUT","R","S",2019,"1"):{"galium aparine"},
         ("AUT","R","S",2019,"3"):{"poa annua"},
         ("AUT","R","T",2018,"1"):{"poa annua"}}
        got=m.select_first_complete(data,2010)
        self.assertEqual(len(got),1)
        self.assertEqual(got[0]["year"],2011)
        after=m.select_first_complete(data,2017)
        self.assertEqual(len(after),1)
        self.assertEqual(after[0]["year"],2019)
    def test_source_system_minimum_prevents_trait_leakage(self):
        self.assertEqual(m.traitmean({"poa annua","lactuca serriola"},EIVE,"position",3),
                         (7,.5,1))
        self.assertEqual(m.traitmean({"poa annua","galium aparine"},EIVE,"width",3),
                         (3.5,1.0,2))
    def test_region_hierarchy_requires_two_roads_and_four_sites_each(self):
        pairs=[]
        for road in ("a","b"):
            for i in range(4):
                pairs.append({"region":"AUT","road":road,"site":str(i),
                  "year":2018,"roadside":{"poa annua"},
                  "inland":{"galium aparine"}})
        panel=m.analyze_panel(pairs,EIVE,"position",3,.5,1)
        self.assertEqual(panel["supported_region_codes"],1)
        self.assertEqual(panel["supported_road_ids"],2)
        self.assertEqual(panel["supported_first_survey_sites"],8)
        self.assertEqual(panel["region_equal_mean_roadside_minus_inland"],2.0)
        self.assertEqual(panel["regions"]["AUT"]["positive_roads"],2)
    def test_single_road_cannot_claim_region_replication(self):
        pairs=[{"region":"AUT","road":"a","site":str(i),"year":2018,
                "roadside":{"poa annua"},"inland":{"galium aparine"}}
                for i in range(4)]
        p=m.analyze_panel(pairs,EIVE,"position",3,.5,1)
        self.assertEqual(p["supported_region_codes"],0)
        self.assertIsNone(p["region_equal_mean_roadside_minus_inland"])
    def test_unmatched_road_side_does_not_generate_zero_plant_plot(self):
        input_={
          ("CHE","r1","A",2020,"1"):{"poa annua"},
          ("CHE","r1","B",2020,"3"):{"galium aparine"}}
        self.assertEqual(m.select_first_complete(input_,2010),[])
    def test_trait_width_not_measured_temperature_celsius(self):
        x=m.traitmean({"poa annua"},EIVE,"width",3)
        self.assertEqual(x[0],3.)
        self.assertNotEqual(x[0],EIVE["poa annua"]["position"])

if __name__=="__main__":unittest.main()
