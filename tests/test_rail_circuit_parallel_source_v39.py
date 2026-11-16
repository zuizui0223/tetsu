"""Check that source parallel track duplication can be separated from real loops."""
import sys
from pathlib import Path
import unittest
from pyproj import Transformer
from shapely.geometry import LineString
from shapely.wkb import dumps

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from audit_rail_circuit_parallel_source_v39 import audit_block,run
from audit_overture_50km_macro_rail_v30 import macro_geometry
from audit_overture_50km_circuit_v37 import load_physical_multigraph_conductance
from audit_gbif_equalarea_macroblocks_v29 import GRID_CRS
from tetsu.circuit_graph import pairwise_effective_resistance


def rail(block,a,b,uid,src="a",dst="b"):
    _,cells=macro_geometry(block)
    togeo=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    p=cells[a].centroid;q=cells[b].centroid
    return (uid,"standard_gauge",[
        {"connector_id":src,"at":0.},
        {"connector_id":dst,"at":1.},
    ],dumps(LineString([togeo.transform(p.x,p.y),togeo.transform(q.x,q.y)])))


class ParallelGraphAuditTests(unittest.TestCase):
    def test_same_source_endpoints_parallel_summed_vs_once_changes_circuit(self):
        block="us_midwest"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        src=[rail(block,a,b,"first"),rail(block,a,b,"second")]
        summed,sa=load_physical_multigraph_conductance(src,parallel_policy="sum")
        one,sb=load_physical_multigraph_conductance(src,parallel_policy="shortest_only")
        rsum,_=pairwise_effective_resistance(summed,["a","b"])
        rone,_=pairwise_effective_resistance(one,["a","b"])
        self.assertAlmostEqual(rsum[("a","b")]/rone[("a","b")],.5,delta=.0001)
        self.assertEqual(sa["n_parallel_connector_links"],1)
        self.assertTrue(sa["source_parallel_linkages_not_dropped_for_circuit"])
        self.assertFalse(sb["source_parallel_linkages_not_dropped_for_circuit"])
        summary=audit_block(block,src)
        result=summary["strategies"]["nearest_in_cell"]
        self.assertAlmostEqual(result["median_Q_with_parallel_sum"],.5,places=6)
        self.assertAlmostEqual(result["median_Q_with_parallel_collapsed"],0.,places=6)
        self.assertEqual(result["n_pairs_Q_collapsed_at_least_5pct"],0)

    def test_distinct_real_route_loops_persist_after_parallel_collapse(self):
        block="fr_burgundy"
        a=f"{block}_x0_y0";b=f"{block}_x1_y0"
        c=f"{block}_x0_y1";d=f"{block}_x1_y1"
        lines=[rail(block,a,b,"ab","a","b"),
               rail(block,b,d,"bd","b","d"),
               rail(block,a,c,"ac","a","c"),
               rail(block,c,d,"cd","c","d")]
        graph,stats=load_physical_multigraph_conductance(lines,parallel_policy="shortest_only")
        self.assertEqual(stats["n_parallel_connector_links"],0)
        resistances,_=pairwise_effective_resistance(graph,["a","d"])
        self.assertIsNotNone(resistances[("a","d")])
        self.assertGreater(resistances[("a","d")],0)
        summary=audit_block(block,lines)
        x=summary["strategies"]["nearest_in_cell"]
        self.assertGreater(x["n_pairs_Q_collapsed_at_least_5pct"],0)
        self.assertLessEqual(x["median_Q_with_parallel_collapsed"],x["median_Q_with_parallel_sum"]+1e-6)

    def test_invalid_policy_does_not_invent_map_evidence(self):
        with self.assertRaises(ValueError):
            load_physical_multigraph_conductance([],parallel_policy="permanent_truth")
        with self.assertRaises(ValueError):
            run(["unfrozen_block"])


if __name__=="__main__":
    unittest.main()
