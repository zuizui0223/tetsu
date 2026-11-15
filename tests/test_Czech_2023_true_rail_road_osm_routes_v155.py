"""v15.5 physical graph design tests, using ONLY synthetic network fixtures."""
import hashlib,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_true_OSM_rail_road_paths_v150 import SourceNetwork
from analyze_Czech_2023_true_rail_road_osm_routes_v155 import (
    graph_station_source,stream_sha256,SOURCE_OSM_DATE
)

def fake_sites(n=35):
    return [{"site":f"X{i}","lat":49.+i*.010,"lon":15.+i*.010,
             "section":str(i%5)} for i in range(n)]

def fake_graph():
    rail,road=SourceNetwork(),SourceNetwork()
    rail.append_way([(1000+i,49.+i*.010,15.+i*.010) for i in range(35)])
    road.append_way([(2000+i,49.+i*.010,15.+i*.010) for i in range(35)])
    return rail,road

class CzechOSMSourceGraphTests(unittest.TestCase):
    def test_full_country_source_date_is_2023_near_2020_to_2022_survey(self):
        self.assertEqual(SOURCE_OSM_DATE,"2023-01-01")

    def test_hash_uses_raw_original_binary_bytes(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"source-osm-fake"
            p.write_bytes(b"abc123")
            self.assertEqual(stream_sha256(p),hashlib.sha256(b"abc123").hexdigest())

    def test_real_node_route_gate_counts_35_source_stations_without_invented_pair_absences(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"fake.osm.pbf"
            p.write_bytes(b"test-data")
            with patch("analyze_Czech_2023_true_rail_road_osm_routes_v155.EXPECTED_SOURCE_PBF_BYTES",len(b"test-data")):
                with patch("analyze_Czech_2023_true_rail_road_osm_routes_v155.read_osm_network",side_effect=lambda path:fake_graph()):
                    report,pairs=graph_station_source(fake_sites(),p,return_pairs=True)
            self.assertEqual(report["n_original_35site_pairs"],595)
            self.assertEqual(report["rail_source_graph"]["n_source_station_centroids_connected_to_osm_node"],35)
            self.assertEqual(report["road_source_graph"]["n_source_station_centroids_connected_to_osm_node"],35)
            self.assertEqual(report["rail_source_graph"]["n_osm_physically_connected_station_pairs"],595)
            self.assertEqual(report["n_original_source_rail_sections"],5)
            self.assertLess(report["n_pairs_with_both_real_2023_physical_routes_and_geo_ge_2km"],595)
            self.assertTrue(report["NO_Czech_railway_dispersal_effect_estimated_from_topology_alone"])

    def test_source_file_size_integrity_check_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"fake.osm.pbf";p.write_bytes(b"a")
            with self.assertRaisesRegex(ValueError,"PBF bytes"):
                graph_station_source(fake_sites(),p,return_pairs=False)

if __name__=="__main__":unittest.main()
