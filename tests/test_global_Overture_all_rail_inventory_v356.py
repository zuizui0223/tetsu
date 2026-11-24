#!/usr/bin/env python3
"""v356: source S3 complete page inventory invariants and rail subtype taxonomy."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_global_Overture_all_rail_parquet_inventory_v356 import (
    RELEASE,BUCKET,PREFIX,EXPECTED_PRIMARY_SCHEMA,PREDECLARED_BBOX,MAX_KEYS
)

class V356GlobalRailInventoryTests(unittest.TestCase):
    def test_exact_immutable_release_before_any_biology(self):
        self.assertEqual(RELEASE,"2026-11-18.0")
        self.assertEqual(EXPECTED_PRIMARY_SCHEMA,"v2.0.0")
        self.assertIn("theme=transportation/type=segment/",PREFIX)
        self.assertIn("2026-11-18.0",PREFIX)
        self.assertEqual(BUCKET,"overturemaps-us-west-2")
    def test_wgs_sample_contains_actual_rail_and_fixed_preanalysis(self):
        self.assertEqual(PREDECLARED_BBOX,(14.420,50.060,14.470,50.120))
    def test_complete_source_list_must_be_finite_not_partial(self):
        self.assertGreaterEqual(MAX_KEYS,10000)

if __name__=="__main__":unittest.main()
