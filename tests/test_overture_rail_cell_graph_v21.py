"""Synthetic WKB and connector tests for Overture source-cell macro network."""
import sys
from pathlib import Path
import unittest
from shapely.geometry import LineString
from shapely.wkb import dumps
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_overture_rail_cell_graph_v21 import (
    rail_query,source_region,graph_component_ids,cell_for_point
)
from harvest_gbif_cell_species_facets_v20 import grid_cells


class FakeDB:
    def __init__(self,rows):self.rows=rows
    def execute(self,_query):return self
    def fetchall(self):return self.rows


class OvertureRailCellGraphTests(unittest.TestCase):
    def test_rail_source_query_is_pinned_and_class_filtered(self):
        from audit_overture_gbif_overlap_v19 import REGIONS
        q=rail_query(REGIONS["jp_sendai"])
        self.assertIn("ST_AsWKB(geometry)",q)
        self.assertIn("connectors",q)
        self.assertIn("2026-10-21.0",q)
        self.assertIn("subtype='rail'",q)
        self.assertNotIn("railway='passenger'",q)

    def test_3x3_city_bin_bounds_and_outside_null(self):
        cells=grid_cells("jp_sendai")
        self.assertEqual(cell_for_point(cells,140.8,38.19),"jp_sendai_x0_y0")
        self.assertIsNone(cell_for_point(cells,0.,0.))

    def test_source_graph_connects_two_distinct_cells_by_connector_edge(self):
        shape=LineString([(140.8,38.2),(140.90,38.28)])
        src=[("r1","standard_gauge",
              [{"connector_id":"a","at":0},{"connector_id":"b","at":1}],dumps(shape))]
        out=source_region(FakeDB(src),"jp_sendai")
        self.assertEqual(out["n_eligible_source_rail_segments"],1)
        self.assertEqual(out["n_connectors_in_source_graph"],2)
        self.assertEqual(out["n_distinct_components_in_source_window"],1)
        self.assertEqual(out["n_cells_with_eligible_connector"],2)
        self.assertEqual(out["n_cell_pairs_connected_in_clipped_source_graph"],1)
        self.assertEqual(out["n_cell_pairs"],36)

    def test_two_crossing_segments_not_automatically_joined(self):
        first=LineString([(140.8,38.20),(140.90,38.28)])
        second=LineString([(140.8,38.28),(140.90,38.20)])
        src=[
            ("r1","standard_gauge",
             [{"connector_id":"a","at":0},{"connector_id":"b","at":1}],dumps(first)),
            ("r2","standard_gauge",
             [{"connector_id":"c","at":0},{"connector_id":"d","at":1}],dumps(second))
        ]
        out=source_region(FakeDB(src),"jp_sendai")
        self.assertEqual(out["n_distinct_components_in_source_window"],2)

    def test_component_grouping_handles_separate_graph_pieces(self):
        x=graph_component_ids({"a":{"b":10},"b":{"a":10},"z":{"x":20},"x":{"z":20}})
        self.assertEqual(sorted(len(a) for a in x),[2,2])

    def test_missing_connectors_fail_source_not_invented(self):
        src=[("r1","standard_gauge",[],dumps(LineString([(140.8,38.2),(140.9,38.28)])))]
        with self.assertRaisesRegex(ValueError,"invalid rail geospatial"):
            source_region(FakeDB(src),"jp_sendai")


if __name__=="__main__":
    unittest.main()
