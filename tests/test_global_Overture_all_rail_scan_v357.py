#!/usr/bin/env python3
"""Synthetic source integrity invariants before world-scale rail GeoParquet scan."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from extract_entire_world_Overture_rail_parquet_v357 import (
    SHARDS,EXPECTED_SOURCE_PARQUET_OBJECTS,EXPECTED_SOURCE_BYTES,
    EXPECTED_SOURCE_INVENTORY_SHA,assigned_files,RAIL_CLASSES,COLUMNS
)
class V357AllSourceShardsTests(unittest.TestCase):
    def setUp(self):
        self.inv={"all_source_parquet_objects":[
           {"key":f"original_{i:03d}.parquet","bytes":567}
           for i in range(128)]}
    def test_fixed_16_shards_cover_exact128_source_files(self):
        parts=[assigned_files(self.inv,j) for j in range(SHARDS)]
        self.assertEqual(len(parts),16)
        self.assertTrue(all(len(x)==8 for x in parts))
        names=[r["key"] for sub in parts for r in sub]
        self.assertEqual(len(names),EXPECTED_SOURCE_PARQUET_OBJECTS)
        self.assertEqual(len(set(names)),128)
    def test_global_source_not_partial_or_replaced(self):
        self.assertEqual(EXPECTED_SOURCE_PARQUET_OBJECTS,128)
        self.assertEqual(EXPECTED_SOURCE_BYTES,72770647616)
        self.assertEqual(len(EXPECTED_SOURCE_INVENTORY_SHA),64)
        with self.assertRaises(ValueError):assigned_files(self.inv,SHARDS)
    def test_rails_include_subways_trams_in_source_census(self):
        self.assertIn("subway",RAIL_CLASSES)
        self.assertIn("tram",RAIL_CLASSES)
        self.assertIn("standard_gauge",RAIL_CLASSES)
    def test_graph_preserves_connectors_and_source_ids(self):
        self.assertIn("connectors",COLUMNS)
        self.assertIn("id",COLUMNS)
        self.assertIn("sources",COLUMNS)
        self.assertIn("geometry",COLUMNS)
if __name__=="__main__":unittest.main()
