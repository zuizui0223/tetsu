#!/usr/bin/env python3
"""v366 independent author original Czech and Slovakia source site-and-rail QA gates."""
import importlib.util,sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_original_CZ_SK_station_GERS_topology_v366 import (
    SOURCE_STUDIES,BBOX_MARGIN_DEGREES,MAX_ORIGINAL_SOURCE_TRANSPORT,
    coordinate_source_window,source_rail_gps_coverage,RELEASE
)

class V366IndependentCzechSlovakRailStationTests(unittest.TestCase):
    def test_original_three_countries_are_real_distinct_field_studies(self):
        self.assertEqual(RELEASE,"2026-11-18.0")
        self.assertEqual(len(SOURCE_STUDIES),2)
        self.assertEqual(SOURCE_STUDIES["CZ"]["n_source_botanical_sites"],35)
        self.assertEqual(SOURCE_STUDIES["CZ"]["n_original_author_rail_stretches"],5)
        self.assertEqual(SOURCE_STUDIES["SK"]["n_source_botanical_sites"],19)
        self.assertEqual(SOURCE_STUDIES["SK"]["n_original_author_rail_stretches"],1)
        self.assertIn("NOT",SOURCE_STUDIES["CZ"]["location_type"])

    def test_original_centroid_vs_DMS_threshold_is_explicit_not_favorably_selected(self):
        self.assertEqual(SOURCE_STUDIES["CZ"]["source_actual_rail_near_max_m"],300.)
        self.assertEqual(SOURCE_STUDIES["SK"]["source_actual_rail_near_max_m"],180.)
        self.assertEqual(SOURCE_STUDIES["CZ"]["source_strict_accuracy_sensitivity_m"],120.)
        self.assertEqual(SOURCE_STUDIES["SK"]["source_strict_accuracy_sensitivity_m"],80.)
        self.assertGreaterEqual(MAX_ORIGINAL_SOURCE_TRANSPORT,20_000_000)

    def test_original_CZ_and_SK_GPS_window_derived_before_plant_similarity(self):
        for sample in [
            [{"lon":14.21,"lat":50.7},{"lon":16.9,"lat":48.8}],
            [{"lon":20.2,"lat":49.13},{"lon":20.8,"lat":49.31}]
        ]:
            w,s,e,n=coordinate_source_window(sample)
            self.assertGreater(e,w)
            self.assertGreater(n,s)
            self.assertGreater(e-max(x["lon"] for x in sample),BBOX_MARGIN_DEGREES-.0001)

    def test_source_physical_station_connectedness_not_fabricated_when_GERS_graph_missing(self):
        src=[{"source_original_group":"A"},{"source_original_group":"A"},
             {"source_original_group":"B"}]
        loc=[{"status":"matched_to_real_rail_connector_network","original_source_track_offset_m":4.},
             {"status":"matched_to_real_rail_connector_network","original_source_track_offset_m":5.},
             {"status":"unmatched_to_global_rail_source","original_source_track_offset_m":None}]
        qa=source_rail_gps_coverage(src,loc,{(0,1):20.,(1,0):20.})
        self.assertEqual(qa["source_potential_original_site_pairs"],3)
        self.assertEqual(qa["source_real_network_connected_original_station_pairs_within_300km"],1)
        self.assertEqual(qa["source_real_network_connected_pairs_SAME_original_author_stretch"],1)
        self.assertEqual(qa["source_real_network_connected_pairs_OTHER_original_author_stretch"],0)
        self.assertEqual(qa["source_original_sites_without_a_qualifying_track_n"],1)

if __name__=="__main__":unittest.main()
