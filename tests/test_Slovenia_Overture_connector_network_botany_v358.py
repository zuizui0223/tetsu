#!/usr/bin/env python3
"""v358: explicit physical network connectors, NOT rail crossing geometry as edges."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_Overture_connector_network_botany_v358 import (
   source_eligible,adjacent_connector_chain,add_segment_edges,build_rail_graph
)
from shapely.geometry import LineString
from pyproj import Transformer
import networkx as nx

class V358OriginalPhysicalRailGraphTests(unittest.TestCase):
    def t(self,a,b):
        return {"id":a,"subtype":"rail","class":"standard_gauge",
                "rail_flags":[],"geometry":b[0].wkb,
                "connectors":b[1]}
    def test_original_rail_classes_not_all_same_line(self):
        self.assertTrue(source_eligible({"subtype":"rail","class":"standard_gauge",
                                        "rail_flags":[{"values":["is_tunnel"]}]}))
        self.assertTrue(source_eligible({"subtype":"rail","class":"standard_gauge",
                                        "rail_flags":[{"values":["is_bridge"]}]}))
        self.assertFalse(source_eligible({"subtype":"rail","class":"tram"}))
        self.assertFalse(source_eligible({"subtype":"rail","class":"standard_gauge",
                                        "rail_flags":[{"values":["is_disused"]}]}))

    def test_source_interior_connector_chain_sorted_and_weighted(self):
        row={"connectors":[{"connector_id":"C","at":1.0},
              {"connector_id":"A","at":0.0},{"connector_id":"B","at":0.5}]}
        chain=adjacent_connector_chain(row)
        self.assertEqual(chain,[(0.0,"A"),(.5,"B"),(1.0,"C")])
        g=nx.Graph()
        add_segment_edges(g,chain,1000.)
        self.assertAlmostEqual(nx.shortest_path_length(g,"A","C",weight="weight"),1000.0)

    def test_true_geo_crossings_without_shared_uuid_remain_disconnected(self):
        tran=Transformer.from_crs("EPSG:4326","EPSG:32633",always_xy=True)
        a=LineString([(14.999,46),(15.001,46)])
        b=LineString([(15,45.999),(15,46.001)])
        rows=[
           self.t("a",(a,[{"connector_id":"a0","at":0.0},{"connector_id":"a1","at":1.0}])),
           self.t("b",(b,[{"connector_id":"b0","at":0.0},{"connector_id":"b1","at":1.0}]))
        ]
        g,segments,qa=build_rail_graph(rows,tran)
        self.assertEqual(g.number_of_edges(),2)
        self.assertEqual(nx.number_connected_components(g),2)
        # Both geometries share exact WGS intersection; network remains split.
        self.assertTrue(a.intersects(b))

    def test_physical_connector_shared_UUID_joins_lines(self):
        tran=Transformer.from_crs("EPSG:4326","EPSG:32633",always_xy=True)
        a=LineString([(14.999,46),(15.001,46)])
        b=LineString([(15,45.999),(15,46.001)])
        rows=[
           self.t("a",(a,[{"connector_id":"a0","at":0.0},
                           {"connector_id":"shared","at":.5},
                           {"connector_id":"a1","at":1.0}])),
           self.t("b",(b,[{"connector_id":"b0","at":0.0},
                           {"connector_id":"shared","at":.5},
                           {"connector_id":"b1","at":1.0}]))
        ]
        g,segments,qa=build_rail_graph(rows,tran)
        self.assertEqual(nx.number_connected_components(g),1)
        self.assertEqual(qa["source_usable_rail_segments"],2)

if __name__=="__main__":unittest.main()
