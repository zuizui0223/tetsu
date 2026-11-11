"""Synthetic full-polygon A5 source coverage vs older centroid-only selection checks."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_public_GBIF_FR_A5_polygon_coverage_v67 import (
    polygon_cells,list_all_source_parquets,analyze,BLOCK,BUCKET,ROOT
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid


class FakeQuery:
    def __init__(self,parent,params):self.parent=parent;self.params=params
    def fetchall(self):return [(123,),(456,)]
    def fetchone(self):return (123,)


class FakeDuck:
    def execute(self,query,params=None):
        if "a5_lonlat_to_cell" in query:
            return FakeQuery(self,params)
        if "a5_geometry_to_cells" in query:
            return FakeQuery(self,params)
        raise AssertionError("Unknown A5 source SQL")


class FakeS3:
    def list_objects_v2(self,**args):
        prefix=args["Prefix"]
        if prefix==ROOT:
            return {"CommonPrefixes":[{"Prefix":ROOT+"taxonpartition=Plantae/"}],
                    "IsTruncated":False}
        if "a5_r2=123/" in prefix:
            return {"Contents":[
                {"Key":prefix+"data_"+str(i)+".parquet","Size":5000000+i}
                for i in range(50)
            ],"IsTruncated":False}
        if "a5_r2=456/" in prefix:
            return {"Contents":[
                {"Key":prefix+"data_"+str(i)+".parquet","Size":500000000+i}
                for i in range(4)
            ],"IsTruncated":False}
        return {"IsTruncated":False,"Contents":[]}


def old():
    keys=[ROOT+"taxonpartition=Plantae/a5_r2=123/data_"+str(i)+".parquet"
          for i in range(48)]
    return {"complete":True,"regions":{
        BLOCK:{
            "n_source_objects_actually_scanned":48,
            "source_object_results":[{"source_key":x} for x in keys],
        }
    }}


class GBIFA5PolygonCompletenessTests(unittest.TestCase):
    def test_polygon_overlap_includes_a5_regions_beyond_centroid(self):
        # Original 16 fixed WKT geometries each cover A5 123 and 456;
        # centroid belongs to A5 123. 456 must not be silently omitted.
        ids,source=polygon_cells(FakeDuck())
        self.assertEqual(ids,[123,456])
        self.assertEqual(source["n_original_50km_source_cells"],16)
        self.assertEqual(source["n_A5_centroid_partition_ids"],1)
        self.assertEqual(source["n_A5_polygon_cover_partition_ids"],2)
        self.assertEqual(source["extra_A5_partition_ids_missed_by_centroids"],[456])

    def test_source_Parquet_metadata_integrity_and_oversized_file_counter(self):
        c=FakeS3()
        a,n=list_all_source_parquets(c,ROOT+"taxonpartition=Plantae/a5_r2=123/")
        b,_=list_all_source_parquets(c,ROOT+"taxonpartition=Plantae/a5_r2=456/")
        self.assertEqual(len(a),50)
        self.assertEqual(len(b),4)
        self.assertGreater(b[0]["size_bytes"],400000000)
        self.assertEqual(n,1)
        with self.assertRaises(ValueError):
            list_all_source_parquets(c,"s3://private/source")

    def test_former_48_objects_are_strict_subset_of_full_country_polygon_source(self):
        result=analyze(old(),FakeDuck(),FakeS3())
        self.assertEqual(result["n_total_source_parquet_files_in_all_true_A5_polygon_intersecting_partitions"],54)
        self.assertEqual(result["n_original_lexically_fixed_48_Parquet_source_files"],48)
        self.assertEqual(result["n_additional_source_files_not_covered_by_old_v65_query"],6)
        self.assertEqual(result["n_large_files_omitted_by_old_400MB_filter"],4)
        self.assertTrue(result["complete_source_object_metadata_listing"])
        self.assertTrue(result["not_a_full_GBIF_original_event_population"])

    def test_fail_closed_if_source_contract_48_keys_missing_in_full_cover(self):
        data=old()
        data["regions"][BLOCK]["source_object_results"][5]["source_key"]=ROOT+"taxonpartition=Plantae/a5_r2=999/not-real.parquet"
        with self.assertRaisesRegex(ValueError,"missing from full polygon"):
            analyze(data,FakeDuck(),FakeS3())


if __name__=="__main__":
    unittest.main()
