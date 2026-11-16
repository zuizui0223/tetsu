"""Offline unit tests for physically source-anchored 50km rail shortest routes."""
import sys
from pathlib import Path
from itertools import combinations
import unittest
from shapely.geometry import LineString
from shapely.wkb import dumps
from pyproj import Transformer

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_overture_50km_rail_paths_v32 import (
    source_network,route_diagnostics,analyze_block,path_m,run
)
from audit_overture_50km_macro_rail_v30 import macro_geometry
from audit_gbif_equalarea_macroblocks_v29 import GRID_CRS


def rail(block,a,b,id,at_extra=False):
    _,polys=macro_geometry(block)
    trans=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    ca=polys[a].centroid;cb=polys[b].centroid
    coords=[trans.transform(ca.x,ca.y),trans.transform(cb.x,cb.y)]
    ends=[{"connector_id":id+"a","at":0.0}]
    if at_extra:
        ends.extend([{"connector_id":id+"p","at":.5},
                     {"connector_id":id+"q","at":.5}])
    ends.append({"connector_id":id+"b","at":1.0})
    return (id,"standard_gauge",ends,dumps(LineString(coords)))


class MacroRailPathsTests(unittest.TestCase):
    def test_one_real_geographic_line_yields_route_at_approx_geodesic(self):
        block="us_midwest"
        line=rail(block,block+"_x0_y0",block+"_x1_y0","r")
        ans=analyze_block(block,[line])
        r=ans["routes"]["nearest_in_cell"]
        self.assertEqual(r["n_geographical_pairs"],120)
        self.assertEqual(r["n_anchored_cells"],2)
        self.assertEqual(r["n_pairs_reachable_in_local_source_graph"],1)
        self.assertEqual(r["n_route_ratios_below_point95_geodesic"],0)
        self.assertAlmostEqual(r["median_rail_to_geodesic_ratio"],1,delta=.02)
        self.assertTrue(ans["not_train_service_or_seed_transfer"])

    def test_equal_at_explicit_nodes_are_source_connected(self):
        block="fr_burgundy"
        item=rail(block,block+"_x0_y0",block+"_x1_y0","s",at_extra=True)
        out=analyze_block(block,[item])
        self.assertEqual(out["source"]["same_at_source_links_epsilon"],1)
        self.assertEqual(out["routes"]["nearest_in_cell"]["n_pairs_reachable_in_local_source_graph"],1)

    def test_two_nonshared_rail_paths_not_fused_at_crossings(self):
        block="jp_tohoku"
        r1=rail(block,block+"_x0_y0",block+"_x1_y0","rail1")
        r2=rail(block,block+"_x0_y1",block+"_x1_y1","rail2")
        out=analyze_block(block,[r1,r2])
        self.assertEqual(out["source"]["n_distinct_source_physical_components"],2)
        self.assertEqual(out["routes"]["nearest_in_cell"]["n_pairs_reachable_in_local_source_graph"],2)
        self.assertEqual(out["routes"]["nearest_in_cell"]["n_anchored_cells"],4)

    def test_graph_shortest_path_targets_and_no_imputed_reachable(self):
        g={"a":{"b":10.,"c":50.},"b":{"a":10.,"c":20.},
           "c":{"b":20.,"a":50.},"z":{"y":10.},"y":{"z":10.}}
        self.assertAlmostEqual(path_m(g,"a",["b","c","z"])["c"],30.)
        self.assertNotIn("z",path_m(g,"a",["z"]))
        with self.assertRaises(ValueError):
            path_m(g,"a",["unknown"])
        with self.assertRaises(ValueError):
            path_m(g,"unknown",[])

    def test_bad_rail_source_and_unknown_blocks_rejected(self):
        block="us_midwest"
        item=rail(block,block+"_x0_y0",block+"_x1_y0","x")
        with self.assertRaisesRegex(ValueError,"Duplicate Overture"):
            source_network(block,[item,item])
        with self.assertRaises(ValueError):
            source_network(block,[])
        with self.assertRaises(ValueError):
            run(["not_a_frozen_block"])


if __name__=="__main__":
    unittest.main()
