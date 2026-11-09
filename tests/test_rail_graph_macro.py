"""Unit-only railway graph topology tests: no biological outcomes or real routes."""
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_graph import (
    validate_segment,build_rail_connector_graph,
    shortest_connector_path_m,graph_components,
)


def segment(uid, a, b, length=100.,subtype="rail",class_name="standard_gauge",middle=None):
    c=[{"connector_id":a,"at":0.}]
    if middle:
        c.append({"connector_id":middle,"at":.25})
    c.append({"connector_id":b,"at":1.})
    return {"id":uid,"subtype":subtype,"class":class_name,
            "length_m":length,"connectors":c}


class RailGraphMacroTests(unittest.TestCase):
    def test_shared_connector_not_geometry_joins(self):
        graph,stats=build_rail_connector_graph([
            segment("r1","a","b",100),
            segment("r2","b","c",200),
            segment("r3","d","e",100),
        ])
        self.assertEqual(shortest_connector_path_m(graph,"a","c"),300.)
        self.assertIsNone(shortest_connector_path_m(graph,"a","e"))
        self.assertEqual(graph_components(graph),[3,2])
        self.assertEqual(stats["n_eligible_segments"],3)

    def test_interior_connector_linear_reference(self):
        graph,stats=build_rail_connector_graph([
            segment("r1","a","z",100,middle="b"),
            segment("r2","b","q",50),
        ])
        self.assertEqual(shortest_connector_path_m(graph,"a","b"),25.)
        self.assertEqual(shortest_connector_path_m(graph,"a","q"),75.)
        self.assertEqual(shortest_connector_path_m(graph,"b","z"),75.)
        self.assertEqual(stats["n_connector_links"],3)

    def test_geometrically_crossing_segments_with_distinct_ids_remain_disconnected(self):
        graph,_=build_rail_connector_graph([
            segment("r1","west","east"),segment("r2","north","south")
        ])
        self.assertIsNone(shortest_connector_path_m(graph,"west","south"))

    def test_parallel_shortest_edge_and_path_limit(self):
        graph,_=build_rail_connector_graph([
            segment("r1","a","b",120),
            segment("r2","a","b",100),
            segment("r3","b","c",45),
        ])
        self.assertEqual(shortest_connector_path_m(graph,"a","c"),145.)
        self.assertIsNone(shortest_connector_path_m(graph,"a","c",120))
        self.assertEqual(shortest_connector_path_m(graph,"a","a"),0.)

    def test_disallowed_subway_not_silently_included(self):
        graph,stats=build_rail_connector_graph([
            segment("r1","a","b",100,class_name="standard_gauge"),
            segment("r2","b","c",200,class_name="subway")
        ])
        self.assertIsNone(shortest_connector_path_m(graph,"a","c"))
        self.assertEqual(stats["n_excluded_classes"],1)
        self.assertEqual(stats["n_source_segments"],2)

    def test_incorrect_schema_or_duplicate_source_fails_closed(self):
        for record in (
            segment("x","a","a"),
            segment("x","a","b",0),
            {"id":"bad","subtype":"rail","class":"standard_gauge","length_m":100,
             "connectors":[{"connector_id":"a","at":0},{"connector_id":"b","at":.5}]},
            segment("x","a","b",-1),
        ):
            with self.assertRaises(ValueError):
                build_rail_connector_graph([record])
        with self.assertRaisesRegex(ValueError,"Duplicate Overture"):
            s=segment("same","a","b")
            build_rail_connector_graph([s,s])
        with self.assertRaises(ValueError):
            validate_segment({"id":"road","subtype":"road"},("standard_gauge",))

    def test_zero_or_negative_distance_not_imputed(self):
        graph,_=build_rail_connector_graph([segment("r1","a","b",50)])
        self.assertIsNone(shortest_connector_path_m(graph,"unknown","b"))
        with self.assertRaises(ValueError):
            shortest_connector_path_m(graph,"a","b",-1)


if __name__=="__main__":
    unittest.main()
