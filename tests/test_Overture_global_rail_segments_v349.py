#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Overture_global_rail_segments_v349 import eligible_segment,ALLOWED_RAIL_CLASS,EXCLUDED_FLAGS

class V349OvertureRailSourceTests(unittest.TestCase):
    def test_frozen_intercity_rail_allowed(self):
        good,_=eligible_segment({"subtype":"rail","class":"standard_gauge",
               "rail_flags":[],"geometry":b"placeholder"})
        self.assertTrue(good)
    def test_tram_subway_excluded(self):
        for cl in ("tram","subway","light_rail","funicular","monorail"):
            ok,why=eligible_segment({"subtype":"rail","class":cl,
                      "geometry":b"a"})
            self.assertFalse(ok);self.assertEqual(why,"wrong_class")
    def test_no_abandoned_tunnels_or_bridges(self):
        for name in EXCLUDED_FLAGS:
            ok,why=eligible_segment({"subtype":"rail","class":"standard_gauge",
               "geometry":b"a","rail_flags":[{"values":[name],"between":None}]})
            self.assertFalse(ok);self.assertEqual(why,"rail_flag_excluded")
    def test_invalid_geometry_excluded_without_imputing(self):
        ok,why=eligible_segment({"subtype":"rail","class":"standard_gauge"})
        self.assertFalse(ok);self.assertEqual(why,"no_geometry")
    def test_nonrail_source_not_conflated(self):
        ok,why=eligible_segment({"subtype":"road","class":"standard_gauge","geometry":b"a"})
        self.assertFalse(ok);self.assertEqual(why,"not_rail")

if __name__=="__main__":unittest.main()
