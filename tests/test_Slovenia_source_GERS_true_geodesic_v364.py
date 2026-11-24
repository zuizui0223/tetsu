#!/usr/bin/env python3
"""v364: synthetic source unit-consistency tests, safe in general CI without GIS."""
import importlib.util,math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from recompute_Slovenia_GERS_true_geodesic_network_v364 import (
    SOURCE_RELEASE,ANALYSIS_NULL_KINDS,source_geodetic_line_m,
    snapped_WGS84_source_at,corrected_physical_graph,corrected_source_attachments
)
HAS_GIS=all(importlib.util.find_spec(x) is not None
            for x in ("networkx","shapely","pyproj"))

@unittest.skipUnless(HAS_GIS,"source true GERS geo-distance tests need GIS in original-source workflow")
class V364TrueGERSGeodeticDistances(unittest.TestCase):
    def setUp(self):
        from pyproj import Transformer,Geod
        from shapely.geometry import LineString
        self.geod=Geod(ellps="WGS84")
        self.tr=Transformer.from_crs("EPSG:4326","EPSG:32633",always_xy=True)
        self.line=LineString([(15.0,46.0),(15.02,46.0),(15.02,46.015)])

    def test_original_Overture_GERS_position_is_exact_true_geodesic_fraction(self):
        L,parts=source_geodetic_line_m(self.line,self.geod)
        x,y=self.line.coords[1]
        f,offset=snapped_WGS84_source_at(self.line,self.tr.transform(x,y),self.tr,self.geod)
        self.assertAlmostEqual(f,parts[0]/L,places=10)
        self.assertAlmostEqual(offset,0.,places=5)
        self.assertAlmostEqual(sum(parts),L,places=8)
        self.assertGreater(L,0)

    def test_nearest_WGS84_segment_geodetic_offset_never_use_simplistic_UTM_total(self):
        lon,lat=self.line.coords[0]
        self.assertAlmostEqual(snapped_WGS84_source_at(
           self.line,self.tr.transform(lon,lat),self.tr,self.geod)[0],0.)
        lon,lat=self.line.coords[-1]
        self.assertAlmostEqual(snapped_WGS84_source_at(
           self.line,self.tr.transform(lon,lat),self.tr,self.geod)[0],1.)

    def test_true_original_connector_graph_preserves_UUID_network_nodes(self):
        from shapely.geometry import LineString
        import networkx as nx
        rows=[{"id":"source-rail-1","subtype":"rail",
          "class":"standard_gauge","rail_flags":[],"geometry":self.line.wkb,
          "connectors":[
             {"connector_id":"a","at":0.0},
             {"connector_id":"b","at":0.45},
             {"connector_id":"c","at":1.0}],
          "names":None}]
        old,new,segments,QA,metric=corrected_physical_graph(rows,self.tr)
        self.assertEqual(old.number_of_nodes(),3)
        self.assertEqual(new.number_of_nodes(),3)
        self.assertEqual(old.number_of_edges(),2)
        self.assertEqual(new.number_of_edges(),2)
        true_length=segments[0]["source_geodetic_length_m"]
        self.assertAlmostEqual(nx.shortest_path_length(new,"a","c",weight="weight"),true_length,places=6)
        self.assertAlmostEqual(new.get_edge_data("a","b")["weight"],true_length*.45,places=6)

    def test_original_site_source_at_and_distances_rewired_to_geodesic(self):
        import networkx as nx
        rows=[{"id":"source-rail-1","subtype":"rail",
           "class":"standard_gauge","rail_flags":[],"geometry":self.line.wkb,
           "connectors":[{"connector_id":"a","at":0.},
                         {"connector_id":"b","at":1.}],
           "names":None}]
        old,new,segments,QA,metric=corrected_physical_graph(rows,self.tr)
        site=[{"lon":15.005,"lat":46.00001},
              {"lon":15.015,"lat":46.00001}]
        aug,attach,qa=corrected_source_attachments(new,segments,site,self.tr)
        self.assertEqual(qa["source_original_site_matches_recalibrated"],2)
        self.assertTrue(all(s["status"]=="matched_to_real_rail_connector_network" for s in attach))
        self.assertTrue(all(s["source_correct_geodetic_linear_referencing"] for s in attach))
        # Site nodes are attached through their two bounding connectors.
        # The downstream original physical_network_distances() includes the
        # shorter direct in-segment path when both sites share one GERS line.
        from analyze_Slovenia_geod_and_rail_distance_matched_v359 import physical_network_distances
        dist=physical_network_distances(aug,attach)[0,1]*1000
        estimated=(attach[1]["source_fraction"]-attach[0]["source_fraction"])*segments[0]["source_geodetic_length_m"]
        self.assertAlmostEqual(dist,estimated,places=5)
if __name__=="__main__":unittest.main()
