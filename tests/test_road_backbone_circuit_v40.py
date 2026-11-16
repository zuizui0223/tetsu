"""Offline major-road backbone resistor sensitivity on real 50km source frame."""
import sys
from pathlib import Path
import unittest
from pyproj import Transformer
from shapely.geometry import LineString
from shapely.wkb import dumps
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_overture_50km_road_backbone_circuit_v40 import (
    MAJOR_ROAD_CLASSES,filtered_query,major_road_circuit_graph,
    analyze_block,run
)
from audit_overture_50km_macro_rail_v30 import macro_geometry
from audit_gbif_equalarea_macroblocks_v29 import GRID_CRS
from tetsu.circuit_graph import pairwise_effective_resistance


def road(block,a,b,uid,source="a",target="b",klass="primary"):
    _,polys=macro_geometry(block)
    projection=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    p=polys[a].centroid
    q=polys[b].centroid
    return (uid,klass,[
        {"connector_id":source,"at":0.},
        {"connector_id":target,"at":1.}
    ],dumps(LineString([
        projection.transform(p.x,p.y),
        projection.transform(q.x,q.y)
    ])))


class MajorRoadCircuitTests(unittest.TestCase):
    def test_class_scope_is_explicit_not_full_road_network(self):
        query=filtered_query("us_midwest")
        self.assertIn("subtype='road'",query)
        for label in MAJOR_ROAD_CLASSES:
            self.assertIn(label,query)
        self.assertNotIn("'residential'",query)
        self.assertNotIn("'tertiary'",query)
        self.assertIn("COUNT(*)",filtered_query("us_midwest",count_only=True))
        with self.assertRaises(ValueError):
            filtered_query("us_midwest",classes=("residential",))

    def test_two_parallel_backbone_links_retained_in_circuit(self):
        block="us_midwest"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        rows=[road(block,a,b,"s1"),road(block,a,b,"s2")]
        graph,stats=major_road_circuit_graph(rows)
        self.assertEqual(stats["n_parallel_connector_pairs_in_source"],1)
        resistances,_=pairwise_effective_resistance(graph,["a","b"])
        self.assertAlmostEqual(resistances[("a","b")],1./graph["a"]["b"],places=6)
        out=analyze_block(block,rows)
        p=out["strategies"]["nearest_in_cell"]
        self.assertEqual(p["n_connected_source_pairs"],1)
        self.assertAlmostEqual(p["median_source_Q_redundancy"],.5,places=6)
        self.assertTrue(out["road_backbone_not_full_network"])

    def test_shortest_path_equals_effective_resistance_for_single_backbone(self):
        block="fr_burgundy"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        out=analyze_block(block,[road(block,a,b,"line")])
        p=out["strategies"]["nearest_in_cell"]
        self.assertEqual(p["n_connected_source_pairs"],1)
        self.assertAlmostEqual(p["median_source_Q_redundancy"],0,places=6)

    def test_source_rejects_outcome_selected_or_missing_regions(self):
        block="us_midwest"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        with self.assertRaises(ValueError):
            major_road_circuit_graph([road(block,a,b,"footway","a","b","residential")])
        with self.assertRaises(ValueError):
            major_road_circuit_graph([road(block,a,b,"dup"),road(block,a,b,"dup")])
        with self.assertRaises(ValueError):
            run(["not_predeclared"])
        with self.assertRaises(ValueError):
            run([])

    def test_no_fake_solvable_road_if_no_source(self):
        with self.assertRaises(ValueError):
            major_road_circuit_graph([])


if __name__=="__main__":
    unittest.main()
