"""Synthetic analytic circuit-theory validations, no plant effects or train claims."""
import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from tetsu.circuit_graph import (
    add_conductance,components,simplify_resistance_graph,
    pairwise_effective_resistance,resistor_edges,write_circuitscape_network
)
from audit_overture_50km_circuit_v37 import load_physical_multigraph_conductance,diagnostic
from audit_overture_50km_macro_rail_v30 import macro_geometry
from audit_gbif_equalarea_macroblocks_v29 import GRID_CRS
from pyproj import Transformer
from shapely.geometry import LineString
from shapely.wkb import dumps


def graph_from_resistors(*edges):
    graph={}
    for a,b,r in edges:
        add_conductance(graph,a,b,1./r)
    return graph


def source_line(block,a,b,sid,connector_start="start",connector_end="end"):
    _,boxes=macro_geometry(block)
    to_geo=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    p=boxes[a].centroid;q=boxes[b].centroid
    shape=LineString([to_geo.transform(p.x,p.y),to_geo.transform(q.x,q.y)])
    return (sid,"standard_gauge",[
        {"connector_id":connector_start,"at":0.0},
        {"connector_id":connector_end,"at":1.0}
    ],dumps(shape))


class NetworkCircuitscapeMathTests(unittest.TestCase):
    def test_series_tree_is_identical_to_shortest_path(self):
        g=graph_from_resistors(("a","b",1.),("b","c",2.))
        eff,meta=pairwise_effective_resistance(g,["a","c"])
        self.assertAlmostEqual(eff[("a","c")],3.,places=7)
        self.assertEqual(meta["n_disconnected_pairs"],0)
        self.assertEqual(meta["n_exact_grounded_solves"],1)
        self.assertEqual(meta["individual_components"][0]["n_reduced_nodes"],2)

    def test_triangle_has_nontrivial_path_redundancy(self):
        g=graph_from_resistors(("a","b",1.),("b","c",1.),("c","a",1.))
        eff,_=pairwise_effective_resistance(g,["a","b","c"])
        for v in eff.values():
            self.assertAlmostEqual(v,2./3.,places=7)
        self.assertEqual(len(resistor_edges(g)),3)

    def test_two_parallel_equal_physical_tracks_halves_resistance(self):
        g=graph_from_resistors(("a","b",1.),("a","b",1.))
        self.assertEqual(len(resistor_edges(g)),1)
        self.assertAlmostEqual(1./g["a"]["b"],.5)
        eff,_=pairwise_effective_resistance(g,["a","b"])
        self.assertAlmostEqual(eff[("a","b")],.5,places=8)

    def test_two_alternative_paths_diamond_reduces_effective_resistance(self):
        g=graph_from_resistors(
            ("a","b",1.),("b","d",1.),
            ("a","c",1.),("c","d",1.)
        )
        eff,_=pairwise_effective_resistance(g,["a","d"])
        self.assertAlmostEqual(eff[("a","d")],1.,places=8)
        simplified,meta=simplify_resistance_graph(g,["a","d"])
        self.assertEqual(len(simplified),2)
        self.assertAlmostEqual(1./simplified["a"]["d"],1.,places=8)
        self.assertEqual(meta["n_exact_series_eliminated"],2)

    def test_disconnected_terminals_are_none_not_ecological_zero(self):
        g=graph_from_resistors(("a","b",1.),("c","d",1.))
        res,meta=pairwise_effective_resistance(g,["a","b","c"])
        self.assertAlmostEqual(res[("a","b")],1.)
        self.assertIsNone(res[("a","c")])
        self.assertIsNone(res[("b","c")])
        self.assertEqual(meta["n_disconnected_pairs"],2)

    def test_no_genuine_source_connector_crossing_joins(self):
        g=graph_from_resistors(("west","east",1.),("north","south",1.))
        self.assertEqual(len(components(g)),2)
        with self.assertRaises(ValueError):
            add_conductance(g,"west","west",1.)
        with self.assertRaises(ValueError):
            add_conductance(g,"west","south",0.)
        with self.assertRaises(ValueError):
            pairwise_effective_resistance(g,["west","west"])

    def test_circuitscape_text_export_has_no_double_counted_undirected_links(self):
        g=graph_from_resistors(("a","b",1.),("b","c",2.),("a","c",3.))
        with tempfile.TemporaryDirectory() as tmp:
            result=write_circuitscape_network(g,["a","c"],tmp,stem="test")
            self.assertTrue(result["julia_circuitscape_not_automatically_executed"])
            self.assertEqual(len(Path(result["network"]).read_text().splitlines()),3)
            self.assertIn("data_type = network",Path(result["config"]).read_text())
            self.assertEqual(len(Path(result["focal"]).read_text().splitlines()),2)

    def test_real_source_loader_sums_two_parallel_routes_not_minimum(self):
        block="us_midwest";a=block+"_x0_y0";b=block+"_x1_y0"
        one=source_line(block,a,b,"s1","aa","bb")
        two=source_line(block,a,b,"s2","aa","bb")
        graph,stats=load_physical_multigraph_conductance([one,two])
        self.assertEqual(stats["n_parallel_connector_links"],1)
        self.assertEqual(stats["multigraph_cycle_rank"],1)
        eff,_=pairwise_effective_resistance(graph,["aa","bb"])
        single,_=load_physical_multigraph_conductance([one])
        self.assertAlmostEqual(eff[("aa","bb")],1./single["aa"]["bb"]/2.,places=6)
        out=diagnostic(block,[one,two])
        q=out["strategies"]["nearest_in_cell"]
        self.assertEqual(q["n_circuit_reachable_source_pairs"],1)
        self.assertEqual(q["n_pairs_with_at_least_5pct_independent_path_redundancy"],1)
        self.assertAlmostEqual(q["median_R_effective_over_shortest"],.5,places=5)
        self.assertFalse(out["actual_Julia_Circuitscape_executed"])

    def test_source_no_rails_is_error_not_effective_zero(self):
        with self.assertRaises(ValueError):
            load_physical_multigraph_conductance([])


if __name__=="__main__":
    unittest.main()
