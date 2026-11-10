"""Offline source tests for geodesically clipped Overture road/rail covariates."""
import sys
from pathlib import Path
import unittest

from shapely.geometry import LineString,box,Point
from shapely.wkb import dumps

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_overture_transport_cells_v26 import (
    query_transport,geo_line_length_km,cell_area_km2,
    process_rows,build_cell_geography,run
)
from audit_overture_gbif_overlap_v19 import REGIONS


class TransportCellSourceTests(unittest.TestCase):
    def test_query_rail_and_road_not_only_rail(self):
        q=query_transport(REGIONS["jp_sendai"])
        self.assertIn("subtype='rail'",q)
        self.assertIn("subtype='road'",q)
        self.assertIn("ST_AsWKB(geometry)",q)
        self.assertIn("2026-10-21.0",q)
        self.assertNotIn("LIMIT 1",q)

    def test_frozen_nine_cell_geodesic_area(self):
        g,a=build_cell_geography("jp_sendai")
        self.assertEqual(len(a),9)
        self.assertTrue(all(0<x<10000 for x in a.values()))
        self.assertAlmostEqual(sum(a.values()),
            cell_area_km2(box(140.77,38.16,141.01,38.37)),delta=1.)
        with self.assertRaises(ValueError):
            build_cell_geography("unfrozen_city")

    def test_clipped_rail_and_road_splits_actual_length(self):
        # One long source line intersects two cells; clipped lengths must not
        # be counted twice as full source segment length.
        line=LineString([(140.8,38.19),(140.94,38.19)])
        rows=[
            ("road_1","road","primary",dumps(line)),
            ("rail_1","rail","standard_gauge",dumps(line)),
        ]
        r=process_rows("jp_sendai",rows)
        self.assertEqual(r["n_source_segment_ids"],2)
        self.assertEqual(r["source_segment_classes"],{"road":1,"rail":1})
        self.assertTrue(r["all_rail_length_km_within_target_cells"]>0)
        self.assertAlmostEqual(r["all_road_length_km_within_target_cells"],
                               r["all_rail_length_km_within_target_cells"],places=5)
        self.assertLess(r["all_road_length_km_within_target_cells"],
                        geo_line_length_km(line))
        self.assertTrue(all(v["area_calculated_geodesically_not_from_degrees"] for v in r["cells"].values()))

    def test_empty_and_point_intersections_zero(self):
        self.assertEqual(geo_line_length_km(Point(140.8,38.2)),0)
        self.assertEqual(geo_line_length_km(LineString()),0)
        self.assertGreater(geo_line_length_km(LineString([(140.8,38.2),(140.81,38.2)])),0)

    def test_malformed_or_duplicate_transport_source_is_rejected(self):
        g=dumps(LineString([(140.8,38.2),(140.9,38.2)]))
        with self.assertRaisesRegex(ValueError,"Duplicate segment"):
            process_rows("jp_sendai",[("id","rail","standard_gauge",g),
                                      ("id","road","primary",g)])
        with self.assertRaisesRegex(ValueError,"Unexpected source"):
            process_rows("jp_sendai",[("id","road","subway",g)])
        with self.assertRaisesRegex(ValueError,"invalid"):
            run(["not_real_city"],db=None)

    def test_missing_source_city_is_not_claimed_complete(self):
        class FakeDB:
            def execute(self,_):return self
            def fetchall(self):return []
            def close(self):pass
        out=run(["jp_sendai"],db=FakeDB())
        self.assertTrue(out["complete"])
        self.assertEqual(out["regions"]["jp_sendai"]["all_road_length_km_within_target_cells"],0)
        # A null spatial return has no proof that the actual city has no
        # mapped roads: API sources need independent expected counts in real run.


if __name__=="__main__":
    unittest.main()
