#!/usr/bin/env python3
"""v360: detect genuine route=railway/tracks ONLY; don't equate trains or crossings."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_original_OSM_named_rail_route_Slovenia_v360 import (
  PBF_EXPECTED_BYTES,pure_source_eligible_railway_relation,
  named_source_railway_relation,source_osm_way_refs,source_rail_name_coverage,
  site_name_coverage
)
class V360NamedRailWaySourceTests(unittest.TestCase):
    def test_source_2026_static_original_pbf_immutable_size(self):
        self.assertEqual(PBF_EXPECTED_BYTES,313730042)
    def test_distinct_infrastructure_rail_relation_not_passenger_train(self):
        self.assertTrue(pure_source_eligible_railway_relation({"type":"route","route":"railway"}))
        self.assertTrue(pure_source_eligible_railway_relation({"type":"route","route":"tracks"}))
        self.assertFalse(pure_source_eligible_railway_relation({"type":"route","route":"train"}))
        self.assertFalse(pure_source_eligible_railway_relation({"type":"multipolygon","route":"railway"}))
    def test_original_named_route_requires_source_relation_not_gers_feature_name(self):
        self.assertFalse(named_source_railway_relation({"type":"route","route":"railway"}))
        self.assertTrue(named_source_railway_relation({"type":"route","route":"tracks","ref":"R15"}))
        self.assertTrue(named_source_railway_relation({"type":"route","route":"railway","name":"Route X"}))
    def test_original_source_osm_way_id_parsed_with_version_not_any_name(self):
        r={"sources":[{"dataset":"OpenStreetMap","record_id":"w12456@31"},
                      {"dataset":"OpenStreetMap","record_id":"w44@0"},
                      {"dataset":"OpenStreetMap","record_id":"r45@9"},
                      {"dataset":"Some-other","record_id":"w933@2"}]}
        self.assertEqual(source_osm_way_refs(r),{12456,44})
    def test_original_GERS_segment_can_have_multiple_OSM_line_relations(self):
        rail=[
          {"id":"a","subtype":"rail","class":"standard_gauge","sources":[
              {"dataset":"OpenStreetMap","record_id":"w11@2"}]},
          {"id":"b","subtype":"rail","class":"tram","sources":[
              {"dataset":"OpenStreetMap","record_id":"w12@2"}]},
          {"id":"c","subtype":"rail","class":"unknown","sources":[
              {"dataset":"OpenStreetMap","record_id":"w13@1"}]},
          {"id":"d","subtype":"rail","class":"narrow_gauge","sources":[
              {"dataset":"OpenStreetMap","record_id":"w14@1"}]}
        ]
        named={11:{111,112},13:{113}}
        all_routes={11:{111,112},13:{113},14:{114}}
        joined,counts=source_rail_name_coverage(rail,named,all_routes)
        self.assertEqual(len(joined),3)
        self.assertEqual(counts["source_eligible_original_GERS_with_multiple_named_OSM_route_relations"],1)
        self.assertEqual(counts["source_eligible_original_GERS_with_named_OSM_route_relation"],2)
        sites=[
          {"status":"matched_to_real_rail_connector_network","source_matched_segment_id_private":"a"},
          {"status":"matched_to_real_rail_connector_network","source_matched_segment_id_private":"c"},
          {"status":"matched_to_real_rail_connector_network","source_matched_segment_id_private":"d"},
          {"status":"ambiguous_multiple_disconnected_source_tracks"}
        ]
        cover=site_name_coverage(sites,joined)
        self.assertEqual(cover["original_sites_with_multiple_named_OSM_route_relations"],1)
        self.assertEqual(cover["original_sites_with_exactly_one_named_OSM_route_relation"],1)
        self.assertEqual(cover["original_sites_with_zero_named_OSM_route_relations"],1)
        self.assertEqual(cover["original_author_sites_not_unambiguously_matched_to_physical_rail"],1)
if __name__=="__main__":unittest.main()
