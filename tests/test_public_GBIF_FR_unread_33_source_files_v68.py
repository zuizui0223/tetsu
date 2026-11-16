"""Source-only synthetic S3 multi-object QA for 33 previously omitted GBIF files."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_public_GBIF_FR_unread_33_source_files_v68 import (
    classify_files,scan_unread,ROOT,BLOCK,
    EXPECTED_ADDITIONAL_FILES,SOURCE_TOTAL_FILES
)


def fixture():
    prefix=ROOT+"taxonpartition=Plantae/a5_r2=123/"
    files=[{"key":prefix+f"data_{i}.parquet","size_bytes":620_000_000 if i>=70 else 50_000_000}
           for i in range(81)]
    previous={"complete":True,"regions":{BLOCK:{
        "n_source_objects_actually_scanned":48,
        "n_actual_interpreted_event_source_rows_in_selected_files_and_bbox":671528,
        "source_object_results":[{"source_key":x["key"]} for x in files[:48]],
    }}}
    geometry={"original_region":BLOCK,
        "n_total_source_parquet_files_in_all_true_A5_polygon_intersecting_partitions":81,
        "source_50km_polygon_plus_centroid_A5_grid_audit":{"polygon_cover_disagreed_with_centroid_lookup":True}}
    return previous,geometry,files


def queried(db,block,item):
    leaf=item["key"].rsplit("/",1)[-1]
    n=int(leaf.removeprefix("data_").removesuffix(".parquet"))
    source_count=100*(n-40) if n%2==0 else 0
    return {
        "n_actual_original_interpreted_rows_in_pilot_bbox":source_count,
        "n_Present_status_source_rows":source_count,
        "n_Absent_status_source_rows":0,
        "n_missing_coordinate_uncertainty_source_rows":source_count//5,
    }


class GBIFUnreadSourceTests(unittest.TestCase):
    def test_33_are_independent_from_original_48_and_source_partition_frozen(self):
        old,geo,files=fixture()
        part=classify_files(old,geo,files)
        self.assertEqual(len(part),EXPECTED_ADDITIONAL_FILES)
        self.assertEqual(len(files),SOURCE_TOTAL_FILES)
        self.assertTrue(all(x["a5_r2"]==123 for x in part))
        self.assertEqual(part[0]["key"],files[48]["key"])

    def test_true_source_positive_zeros_are_all_retained_and_not_species_count(self):
        old,geo,files=fixture()
        result=scan_unread(None,old,geo,files,query=queried)
        self.assertTrue(result["complete"])
        self.assertEqual(result["n_new_additional_source_Parquet_objects_successful"],33)
        self.assertGreater(result["2018to2025_raw_all_Plantae_original_rows_from_new_33_files"],0)
        self.assertLess(result["fraction_of_combined_bounded_81_file_rows_covered_by_old_48"],1)
        self.assertTrue(result["all_source_country_GRIIS_taxon_outcomes_UNOBSERVED"])
        self.assertTrue(result["not_a_globally_COMPLETE_GBIF_census_or_new_DOI"])
        self.assertEqual(len(result["source_additional_object_count_diagnostics"]),33)

    def test_bad_source_object_missing_or_too_large_fails_closed(self):
        old,geo,files=fixture()
        files[60]["size_bytes"]=1_000_000_001
        with self.assertRaisesRegex(ValueError,"1GB"):
            classify_files(old,geo,files)
        old,geo,files=fixture()
        files.pop()
        with self.assertRaises(ValueError):
            classify_files(old,geo,files)

    def test_remote_file_error_yields_explicit_incomplete_not_zero(self):
        old,geo,files=fixture()
        def broken(db,block,item):
            if item["key"].endswith("data_52.parquet"):
                raise RuntimeError("source connection failed")
            return queried(db,block,item)
        result=scan_unread(None,old,geo,files,query=broken)
        self.assertFalse(result["complete"])
        self.assertEqual(result["n_source_files_unread_or_invalid"],1)
        self.assertIsNone(result["2018to2025_combined_81_file_raw_all_Plantae_source_rows_in_same_bbox"])


if __name__=="__main__":
    unittest.main()
