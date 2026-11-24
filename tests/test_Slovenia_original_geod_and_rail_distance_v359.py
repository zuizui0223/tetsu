#!/usr/bin/env python3
"""v359: no source site identity or species can define physical match; fail closed."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_geod_and_rail_distance_matched_v359 import (
  GEO_CALIPER_KM,NET_CALIPER_KM,GEO_MAX_KM,
  joint_geo_and_network_triples,botanical_comparison,MATCH_MODES
)

class V359TrueNetworkMatchedRailBotanyTests(unittest.TestCase):
    def setUp(self):
        self.sites=[
          {"lat":46., "lon":15.+i*.001,"section":"A" if i<3 else "B",
           "year":2022,"layout":"SINGLE",
           "BS":({"poa","bromus","trifolium"} if i<3 else {"poa","saponaria"}),
           "BT":{"galium","poa"} if i<3 else {"poa","rumex"}}
          for i in range(6)]
        self.d={(i,j):abs(i-j)*.35+1 for i in range(6) for j in range(6) if i!=j}

    def test_network_and_geo_matching_source_fixed_limits(self):
        self.assertEqual(GEO_CALIPER_KM,2.5)
        self.assertEqual(NET_CALIPER_KM,5.0)
        self.assertEqual(GEO_MAX_KM,40.0)
        triads,qa=joint_geo_and_network_triples(self.sites,self.d,"geo_and_network")
        self.assertGreater(len(triads),0)
        for x in triads:
            self.assertNotEqual(x["same"],x["different"])
            self.assertEqual(self.sites[x["anchor"]]["section"],self.sites[x["same"]]["section"])
            self.assertNotEqual(self.sites[x["anchor"]]["section"],self.sites[x["different"]]["section"])
            self.assertLessEqual(x["source_abs_network_gap_km"],NET_CALIPER_KM)
            self.assertLessEqual(x["source_abs_geographic_gap_km"],GEO_CALIPER_KM)

    def test_no_network_path_is_not_a_distance_zero(self):
        triads,qa=joint_geo_and_network_triples(self.sites,{},"geo_and_network")
        self.assertEqual(triads,[])
        self.assertGreater(qa["no_source_Oveture_shortest_path_or_unmatched"],0)
        summary=botanical_comparison(self.sites,triads)
        self.assertFalse(summary["source_qualifies_predefined_minimum_12_anchors_4_sections"])
        self.assertEqual(summary["biological_result"],"NOT_IDENTIFIED_due_to_source_support")

    def test_site_year_layout_control(self):
        sites=[dict(s) for s in self.sites]
        sites[3]["year"]=2023
        triads,qa=joint_geo_and_network_triples(sites,self.d,"geo_network_year_layout")
        for x in triads:
            a=sites[x["anchor"]]
            for target in (x["same"],x["different"]):
                self.assertEqual(a["year"],sites[target]["year"])
                self.assertEqual(a["layout"],sites[target]["layout"])

    def test_source_plant_species_names_not_used_for_match(self):
        first,_=joint_geo_and_network_triples(self.sites,self.d,"geo_and_network")
        altered=[dict(s,BS={"totally","unrelated"},BT={"rare"}) for s in self.sites]
        second,_=joint_geo_and_network_triples(altered,self.d,"geo_and_network")
        self.assertEqual(first,second)

if __name__=="__main__":unittest.main()
