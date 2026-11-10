"""Synthetic topology and route-distance tests for Overture road negative control."""
import sys
from pathlib import Path
import unittest
from itertools import combinations
from shapely.geometry import LineString
from shapely.wkb import dumps

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_overture_road_network_v28 import (
    query_road_source,graph_from_source_rows,graph_components,
    shortest_paths_to_targets,network_for_city,haversine_m
)
from audit_overture_gbif_overlap_v19 import REGIONS
from harvest_gbif_cell_species_facets_v20 import grid_cells


def road(sid, lonlat, a, b):
    return (sid,"residential",[
        {"connector_id":a,"at":0.0},{"connector_id":b,"at":1.0}
    ],dumps(LineString(lonlat)))


class RoadNetworkTests(unittest.TestCase):
    def test_source_road_sql_pin_schema(self):
        q=query_road_source(REGIONS["jp_sendai"])
        self.assertIn("subtype='road'",q)
        self.assertIn("connectors",q)
        self.assertIn("ST_AsWKB(geometry)",q)
        self.assertIn("2026-10-21.0",q)
        self.assertNotIn("subtype='rail'",q)

    def test_road_share_connector_but_crossed_lines_not_connected(self):
        g,pos,meta=graph_from_source_rows([
            road("r1",[(140.79,38.18),(140.90,38.18)],"a","b"),
            road("r2",[(140.90,38.18),(140.96,38.18)],"b","c"),
            road("r3",[(140.88,38.13),(140.88,38.25)],"x","y")
        ])
        groups,sizes=graph_components(g)
        self.assertEqual(len(sizes),2)
        self.assertEqual(groups["a"],groups["c"])
        self.assertNotEqual(groups["a"],groups["y"])
        result=shortest_paths_to_targets(g,"a",["b","c"])
        self.assertGreater(result["c"],result["b"])
        self.assertEqual(meta["n_eligible_source_road_segments"],3)

    def test_nine_cell_grid_route_matrix_on_connected_synthetic_roads(self):
        cells=grid_cells("jp_sendai")
        centroids={}
        for c in cells:
            w,s,e,n=c["bbox"]
            centroids[c["cell_id"]]=((w+e)/2,(s+n)/2)
        roads=[]
        for a,b in combinations(cells,2):
            ia,ib=a["cell_id"],b["cell_id"]
            ax,ay=int(ia[-4]),int(ia[-1])
            bx,by=int(ib[-4]),int(ib[-1])
            if abs(ax-bx)+abs(ay-by)==1:
                roads.append(road("link_"+ia+"_"+ib,
                                  [centroids[ia],centroids[ib]],ia,ib))
        actual=network_for_city("jp_sendai",roads)
        self.assertEqual(actual["n_eligible_source_road_segments"],12)
        self.assertEqual(actual["n_source_cells_snapped_within_5km"],9)
        self.assertEqual(actual["n_source_cell_pairs_road_connected"],36)
        self.assertEqual(len(actual["cell_pair_road_paths"]),36)
        self.assertTrue(all(x["optimistic_road_connector_path_km"] is not None for x in actual["cell_pair_road_paths"]))
        self.assertTrue(all(x["optimistic_road_to_euclidean_ratio"]>=.9 for x in actual["cell_pair_road_paths"]))
        self.assertTrue(actual["window_clipped_unknown_outside_links"])

    def test_invalid_road_source_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"Duplicate Overture"):
            graph_from_source_rows([road("x",[(140.8,38.2),(140.9,38.2)],"a","b"),
                                    road("x",[(140.8,38.2),(140.9,38.2)],"a","b")])
        with self.assertRaisesRegex(ValueError,"connector"):
            graph_from_source_rows([
                ("x","residential",[{"connector_id":"a","at":0.0}],
                 dumps(LineString([(140.8,38.2),(140.9,38.2)])))
            ])
        with self.assertRaises(ValueError):
            graph_from_source_rows([])

    def test_impossible_paths_and_bad_weights(self):
        graph={"a":{"b":1.},"b":{"a":1.},"z":{"y":1.},"y":{"z":1.}}
        self.assertEqual(shortest_paths_to_targets(graph,"a",["b"])["b"],1.0)
        self.assertNotIn("y",shortest_paths_to_targets(graph,"a",["y"]))
        with self.assertRaises(ValueError):
            shortest_paths_to_targets(graph,"unknown",["a"])
        bad={"a":{"b":-1.},"b":{"a":-1.}}
        with self.assertRaises(ValueError):
            shortest_paths_to_targets(bad,"a",["b"])

    def test_near_centers_and_source_country_not_final_botanical_effect(self):
        self.assertEqual(haversine_m((140.,38.),(140.,38.)),0.)
        self.assertGreater(haversine_m((140.,38.),(141.,38.)),80000.)
        self.assertLess(haversine_m((140.,38.),(141.,38.)),90000.)


if __name__=="__main__":
    unittest.main()
