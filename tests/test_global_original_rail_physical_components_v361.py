#!/usr/bin/env python3
"""v361: exact GERS physical topology tests; crossings without connector IDs split."""
import sys,unittest,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from build_global_Overture_physical_connector_components_v361 import (
  Disjoint,original_connectors,eligible_panels,PANELS,NORMAL_CLASSES,
  SOURCE_RELEASE,EXPECTED_INVENTORY_SHA
)
class V361GenuineGlobalPhysicalRailGraphTests(unittest.TestCase):
    def test_frozen_world_source_release_and_two_graph_panels(self):
        self.assertEqual(SOURCE_RELEASE,"2026-11-18.0")
        self.assertEqual(len(EXPECTED_INVENTORY_SHA),64)
        self.assertIn("all_original_rail_subtypes",PANELS)
        self.assertIn("conventional_active_candidate_source_rails",PANELS)

    def test_only_shared_real_connector_UUID_connects_crossing_lines(self):
        uf=Disjoint()
        for a,b in [("a0","a1"),("b0","b1")]:uf.union(a,b)
        self.assertEqual(uf.components,2)
        uf.union("a1","b0")
        self.assertEqual(uf.components,1)

    def test_original_connector_chain_sorted_and_duplicates_removed(self):
        row={"connectors":[{"connector_id":"c2","at":1},
          {"connector_id":"c0","at":0},
          {"connector_id":"c1","at":.5},
          {"connector_id":"c1","at":.5},
          {"connector_id":None,"at":.66},
          {"connector_id":"bad","at":2.0}]}
        self.assertEqual(original_connectors(row),["c0","c1","c2"])

    def test_normal_rail_keeps_bridges_tunnels_but_not_abandoned(self):
        original={"subtype":"rail","class":"standard_gauge",
                  "rail_flags":[{"values":["is_bridge","is_tunnel"]}]}
        p=eligible_panels(original)
        self.assertTrue(p["all_original_rail_subtypes"])
        self.assertTrue(p["conventional_active_candidate_source_rails"])
        bad=copy.deepcopy(original);bad["rail_flags"]=[{"values":["is_abandoned"]}]
        self.assertTrue(eligible_panels(bad)["all_original_rail_subtypes"])
        self.assertFalse(eligible_panels(bad)["conventional_active_candidate_source_rails"])
        tram=copy.deepcopy(original);tram["class"]="tram"
        self.assertTrue(eligible_panels(tram)["all_original_rail_subtypes"])
        self.assertFalse(eligible_panels(tram)["conventional_active_candidate_source_rails"])
        self.assertEqual(NORMAL_CLASSES,frozenset(("standard_gauge","narrow_gauge","unknown")))

    def test_true_connector_names_not_regional_closeness(self):
        uf=Disjoint()
        for n in ["left","right","isolated"]:uf.add(n)
        uf.union("left","right")
        self.assertEqual(uf.find("left"),uf.find("right"))
        self.assertNotEqual(uf.find("isolated"),uf.find("left"))

if __name__=="__main__":unittest.main()
