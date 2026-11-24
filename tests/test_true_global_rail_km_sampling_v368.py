#!/usr/bin/env python3
"""v368: GLOBAL railway-km PPSWR sampling independent of plant observations."""
import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from sample_entire_global_GERS_railkm_WGS84_v368b import (
    SAMPLE_K,SEED,EXPECTED_CANDIDATE_SEGMENTS,EXPECTED_CANDIDATE_KM,
    fixed_sample_targets,candidate_segment
)

class GlobalRailLengthProbabilitySourceTests(unittest.TestCase):
    def test_population_is_all_global_active_candidate_normal_rail_km(self):
        self.assertEqual(SAMPLE_K,2048)
        self.assertEqual(EXPECTED_CANDIDATE_SEGMENTS,1536226)
        self.assertEqual(EXPECTED_CANDIDATE_KM,2237043.481)
        self.assertIn("rail-km-global-random",SEED)
    def test_genuine_world_rail_km_probability_samples_with_replacement(self):
        a=fixed_sample_targets(2_237_043_481.)
        b=fixed_sample_targets(2_237_043_481.)
        self.assertEqual(a,b)
        self.assertEqual(len(a),2048)
        self.assertEqual({i for _,i in a},set(range(2048)))
        self.assertTrue(all(0<=v<2_237_043_481. for v,_ in a))
        self.assertEqual(a,sorted(a))
    def test_not_named_lines_no_tunnel_filter_but_exclude_inactive(self):
        row={"subtype":"rail","class":"standard_gauge",
             "rail_flags":[{"values":["is_bridge","is_tunnel"]}]}
        self.assertTrue(candidate_segment(row))
        row["rail_flags"]=[{"values":["is_disused"]}]
        self.assertFalse(candidate_segment(row))
        row["rail_flags"]=[]
        row["class"]="tram"
        self.assertFalse(candidate_segment(row))
        row["class"]="unknown"
        self.assertTrue(candidate_segment(row))
if __name__=="__main__":unittest.main()
