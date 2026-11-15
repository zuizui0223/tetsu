"""v15.0 synthetic topology tests; actual output requires publisher original+Geofabrik PBF."""
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_true_OSM_rail_road_paths_v150 import (
    SourceNetwork,snap_sites,pair_paths,MAIN_RAIL_MAX_SNAP_KM,
    MAIN_ROAD_MAX_SNAP_KM,ROAD_KINDS,RAIL_KINDS,
)

class RealGraphDesignTests(unittest.TestCase):
    def test_physical_crossing_without_shared_original_osm_node_remains_disconnected(self):
        g=SourceNetwork()
        g.append_way([(1,46.000,15.000),(2,46.010,15.010)])
        g.append_way([(3,46.010,15.000),(4,46.000,15.010)])
        graph,coords,components,qc=g.finalize()
        self.assertEqual(qc["original_osm_unique_node_count"],4)
        self.assertEqual(qc["original_osm_connected_components"],2)
        s=[{"lat":46.0,"lon":15.0},{"lat":46.010,"lon":15.0}]
        snap=snap_sites(coords,s,0.5)
        dist=pair_paths(graph,snap)
        self.assertTrue(math.isinf(dist[0,1]))
        self.assertNotEqual(components[snap[0]["node"]],components[snap[1]["node"]])

    def test_route_distance_exceeds_geodesic_and_uses_real_node_connectors(self):
        g=SourceNetwork()
        # A hairpin connected by node ID 2 only, unlike straight-line shortcut.
        g.append_way([(1,46.0,15.0),(2,46.015,15.015),(3,46.0,15.030)])
        graph,coords,components,qc=g.finalize()
        sites=[{"lat":46.0,"lon":15.0},{"lat":46.0,"lon":15.030}]
        snaps=snap_sites(coords,sites,0.5)
        dist=pair_paths(graph,snaps)
        self.assertTrue(math.isfinite(dist[0,1]))
        self.assertGreater(dist[0,1],2.3)
        self.assertAlmostEqual(dist[0,1],dist[1,0],places=9)
        self.assertAlmostEqual(dist[0,0],0.0,places=9)

    def test_two_original_ways_join_only_on_same_node_id(self):
        g=SourceNetwork()
        g.append_way([(1,46,15),(2,46,15.01)])
        g.append_way([(2,46,15.01),(3,46.01,15.01)])
        graph,_,components,qc=g.finalize()
        self.assertEqual(qc["original_osm_connected_components"],1)
        self.assertEqual(qc["original_osm_unique_undirected_edges"],2)

    def test_nearest_network_node_gate_does_not_invent_site_connector(self):
        g=SourceNetwork()
        g.append_way([(1,46,15),(2,46,15.01)])
        graph,coords,components,qc=g.finalize()
        snap=snap_sites(coords,[{"lat":46.02,"lon":15.02},{"lat":46,"lon":15}],0.5)
        self.assertIsNone(snap[0]["node"])
        self.assertTrue(math.isfinite(snap[0]["distance_km"]))
        self.assertTrue(snap[1]["on_supported_source_network"])
        self.assertTrue(math.isinf(pair_paths(graph,snap)[0,1]))

    def test_parallel_edge_duplicated_osm_way_preserves_shortest_length(self):
        g=SourceNetwork()
        g.append_way([(1,46,15),(2,46,15.01)])
        g.append_way([(2,46,15.01),(1,46,15)])
        graph,coords,_,qc=g.finalize()
        self.assertEqual(qc["original_osm_unique_undirected_edges"],1)
        s=snap_sites(coords,[{"lat":46,"lon":15},{"lat":46,"lon":15.01}],0.5)
        self.assertAlmostEqual(pair_paths(graph,s)[0,1],0.773,delta=0.04)

    def test_predeclared_source_tag_filters_reject_footpath_and_non_active_rail(self):
        self.assertIn("rail",RAIL_KINDS)
        self.assertNotIn("abandoned",RAIL_KINDS)
        self.assertIn("service",ROAD_KINDS)
        self.assertNotIn("footway",ROAD_KINDS)
        self.assertNotIn("construction",ROAD_KINDS)

if __name__=="__main__":unittest.main()
