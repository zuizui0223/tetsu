"""Synthetic 50km equal-area physical ROAD path and source QA."""
import sys
from pathlib import Path
from itertools import combinations
import unittest
from pyproj import Transformer
from shapely.geometry import LineString
from shapely.wkb import dumps
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_overture_50km_road_paths_v33 import (
    road_query,network_source,route_summary,run,_dijkstra
)
from audit_overture_50km_macro_rail_v30 import macro_geometry
from audit_gbif_equalarea_macroblocks_v29 import GRID_CRS


def road(block,cell_a,cell_b,sid,mid=False):
    _,cells=macro_geometry(block)
    w=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    p=cells[cell_a].centroid;q=cells[cell_b].centroid
    links=[{"connector_id":sid+"a","at":0.0}]
    if mid:
        links.extend([{"connector_id":sid+"p","at":0.5},
                      {"connector_id":sid+"q","at":0.5}])
    links.append({"connector_id":sid+"b","at":1.0})
    return (sid,"primary",links,dumps(LineString([
        w.transform(p.x,p.y),w.transform(q.x,q.y)
    ])))


class RoadMacroPathTests(unittest.TestCase):
    def test_source_sql_only_road_and_halo(self):
        q=road_query("us_midwest")
        self.assertIn("subtype='road'",q)
        self.assertIn("ST_AsWKB(geometry)",q)
        self.assertIn("2026-10-21.0",q)
        self.assertNotIn("subtype='rail'",q)

    def test_connected_two_cells_have_same_geodesic_network_route(self):
        block="us_midwest"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        g,anchors,meta=network_source(block,[road(block,a,b,"r")])
        result=route_summary(block,g,anchors)
        q=result["nearest_in_cell"]
        self.assertEqual(q["n_anchored_cells"],2)
        self.assertEqual(q["n_road_route_pairs"],1)
        self.assertAlmostEqual(q["median_road_to_geodesic_ratio"],1,delta=.02)
        self.assertEqual(meta["n_source_road_segments_within_bbox_plus_halo"],1)

    def test_equal_at_different_connector_ids_are_local_road_edges(self):
        block="fr_burgundy"
        one=road(block,f"{block}_x0_y0",f"{block}_x1_y0","same_at",True)
        graph,anchors,m=network_source(block,[one])
        self.assertEqual(m["n_road_equal_at_connector_links_epsilon"],1)
        self.assertEqual(route_summary(block,graph,anchors)["nearest_in_cell"]["n_road_route_pairs"],1)

    def test_separate_road_components_not_connected_by_map_crossings(self):
        block="fr_burgundy"
        a=road(block,f"{block}_x0_y0",f"{block}_x1_y0","first")
        b=road(block,f"{block}_x0_y1",f"{block}_x1_y1","second")
        graph,anchors,meta=network_source(block,[a,b])
        self.assertEqual(meta["n_source_graph_components"],2)
        self.assertEqual(route_summary(block,graph,anchors)["nearest_in_cell"]["n_road_route_pairs"],2)

    def test_dijkstra_and_unknown_source(self):
        g={"a":{"b":12.,"c":50.},"b":{"a":12.,"c":17.},
           "c":{"a":50.,"b":17.}}
        self.assertEqual(_dijkstra(g,"a",["c"])["c"],29.)
        with self.assertRaises(ValueError):
            _dijkstra(g,"missing",[])

    def test_incomplete_endpoint_source_rows_explicitly_skipped_only_under_two_percent(self):
        block="gb_yorkshire"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        actual=[road(block,a,b,str(i)) for i in range(100)]
        invalid=list(road(block,a,b,"incomplete"))
        invalid[2]=[
            {"connector_id":"internal_start","at":0.25},
            {"connector_id":"internal_end","at":1.0}
        ]
        g,anchors,meta=network_source(block,actual+[tuple(invalid)])
        self.assertEqual(meta["n_excluded_source_road_segments_missing_endpoints"],1)
        self.assertEqual(meta["n_original_source_road_segments_returned"],101)
        self.assertLess(meta["excluded_incomplete_road_fraction"],.02)
        self.assertEqual(meta["n_source_road_segments_within_bbox_plus_halo"],100)
        self.assertTrue(meta["incomplete_road_source_NOT_interpreted_as_biological_absence"])

    def test_incomplete_source_fraction_over_two_percent_stops_without_fabricated_links(self):
        block="au_southeast"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        valid=road(block,a,b,"valid")
        bad=list(road(block,a,b,"bad"))
        bad[2]=[
            {"connector_id":"only_inside","at":.2},
            {"connector_id":"real_end","at":1.}
        ]
        with self.assertRaisesRegex(ValueError,"exceeds fixed 2%"):
            network_source(block,[valid,tuple(bad)])

    def test_missing_or_duplicate_source_never_claimed_road_free(self):
        block="us_midwest";a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        one=road(block,a,b,"same")
        with self.assertRaisesRegex(ValueError,"Missing Overture"):
            network_source(block,[])
        with self.assertRaisesRegex(ValueError,"Duplicate road"):
            network_source(block,[one,one])
        with self.assertRaises(ValueError):
            run(["unknown_block"])


if __name__=="__main__":
    unittest.main()
