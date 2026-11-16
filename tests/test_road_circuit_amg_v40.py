"""Optional AMG electrical solver parity test: runs in FR source workflow.

General repository smoke CI can skip this if pyamg is unavailable;
the FR workflow installs pyamg and executes it before actual source.
"""
import sys,importlib.util
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.circuit_graph import add_conductance,pairwise_effective_resistance


def graph(*edges):
    g={}
    for a,b,r in edges:
        add_conductance(g,a,b,1./r)
    return g


@unittest.skipUnless(importlib.util.find_spec("pyamg") is not None,
                     "Optional pyamg is only installed in FR major ROAD source workflow")
class CircuitAMGTests(unittest.TestCase):
    def test_direct_vs_amg_triangle(self):
        g=graph(("a","b",1.),("b","c",1.),("a","c",1.))
        direct,_=pairwise_effective_resistance(g,["a","b","c"],solver="direct")
        approx,meta=pairwise_effective_resistance(g,["a","b","c"],solver="amg")
        self.assertEqual(meta["numerical_solver"],"amg")
        self.assertTrue(meta["AMG_is_iterative_with_explicit_residual_checks"])
        for p in direct:
            self.assertAlmostEqual(approx[p],direct[p],delta=1e-6)

    def test_direct_vs_amg_diamond_tree_and_parallel(self):
        cases=[
            (graph(("a","b",1.),("b","d",1.),("a","c",1.),("c","d",1.)),["a","d"]),
            (graph(("a","b",1.),("b","c",2.)),["a","c"]),
            (graph(("a","b",1.),("a","b",1.)),["a","b"])
        ]
        for g,terminals in cases:
            d,_=pairwise_effective_resistance(g,terminals,solver="direct")
            a,_=pairwise_effective_resistance(g,terminals,solver="amg")
            for pair in d:
                self.assertAlmostEqual(a[pair],d[pair],delta=1e-6)

    def test_solver_does_not_silently_change(self):
        g=graph(("a","b",1.))
        with self.assertRaises(ValueError):
            pairwise_effective_resistance(g,["a","b"],solver="unsupported")


if __name__=="__main__":
    unittest.main()
