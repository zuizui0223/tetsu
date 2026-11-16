"""Actual public GBIF A5 multi-file source diagnosis, never assume zero after first five."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from probe_public_GBIF_A5_bounded_scan_v65 import (
    MAX_BOUNDED_SOURCE_OBJECTS,MIN_BOUNDED_SOURCE_OBJECTS,
    scan_bounded_source,BUCKET,ROOT
)


def fake_candidates(n=9):
    return [{"key":ROOT+f"taxonpartition=Plantae/a5_r2=123/data_{i}.parquet",
             "size_bytes":100000+i}
            for i in range(n)]


def queried(db,block,candidate):
    import re
    i=int(candidate["key"].split("data_")[-1].split(".")[0])
    count=33 if i==6 else 0
    return {
        "selected_public_parquet_key":candidate["key"],
        "selected_public_parquet_size_bytes":candidate["size_bytes"],
        "n_actual_original_interpreted_rows_in_pilot_bbox":count,
        "n_Present_status_source_rows":count,
        "n_Absent_status_source_rows":0,
        "n_missing_coordinate_uncertainty_source_rows":3 if count else 0,
        "n_distinct_GBIF_dataset_UUIDs_in_sample":2 if count else 0,
    }


class ActualMultiFileGBIFQueryTests(unittest.TestCase):
    def test_source_scan_does_not_stop_after_five_empty_files(self):
        # The historic France v6.4 run gave zero in its first five files.
        # Sixth file is positive in this SYNTHETIC fixture, and all
        # 9 source files are audited, not selected for a good outcome.
        row=scan_bounded_source(None,"fr_burgundy",fake_candidates(9),
                                query=queried,scan_limit=MAX_BOUNDED_SOURCE_OBJECTS)
        self.assertEqual(row["n_source_objects_actually_scanned"],9)
        self.assertEqual(row["n_source_objects_with_zero_actual_source_occurrences"],8)
        self.assertEqual(row["n_source_objects_with_nonzero_actual_source_occurrences"],1)
        self.assertEqual(row["n_actual_interpreted_event_source_rows_in_selected_files_and_bbox"],33)
        self.assertTrue(row["full_centroid_partition_source_object_list_covered"])
        self.assertTrue(row["genuine_GBIF_original_interpreted_rows_verified"])
        self.assertFalse(row["all_200km_original_16_cell_geography_fully_observed"])

    def test_truncation_is_recorded_and_zero_is_not_biological_absence(self):
        rows=scan_bounded_source(None,"fr_burgundy",fake_candidates(55),
                                 query=queried,scan_limit=48)
        self.assertEqual(rows["source_candidate_objects_in_centroid_A5_bound"],55)
        self.assertEqual(rows["n_source_objects_actually_scanned"],48)
        self.assertFalse(rows["full_centroid_partition_source_object_list_covered"])
        self.assertEqual(rows["n_source_objects_with_nonzero_actual_source_occurrences"],1)

    def test_invalid_source_paths_or_too_few_files_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"Insufficient public"):
            scan_bounded_source(None,"fr_burgundy",fake_candidates(5),query=queried)
        c=fake_candidates(10)
        c[5]["key"]="s3://wrong-private-source/data_5.parquet"
        with self.assertRaisesRegex(ValueError,"prefix"):
            scan_bounded_source(None,"fr_burgundy",c,query=queried)
        with self.assertRaisesRegex(ValueError,"original frozen"):
            scan_bounded_source(None,"unselected_block",fake_candidates(10),query=queried)
        with self.assertRaisesRegex(ValueError,"pre-outcome"):
            scan_bounded_source(None,"fr_burgundy",fake_candidates(10),query=queried,scan_limit=5)

    def test_any_IO_failure_is_reported_not_silently_reclassified_as_zero(self):
        def broken(db,block,candidate):
            if candidate["key"].endswith("/data_4.parquet"):
                raise RuntimeError("synthetic source corruption")
            return queried(db,block,candidate)
        row=scan_bounded_source(None,"fr_burgundy",fake_candidates(10),query=broken)
        self.assertEqual(row["n_source_objects_failed_API_or_IO"],1)
        self.assertFalse(row["genuine_GBIF_original_interpreted_rows_verified"])
        self.assertFalse(row["full_centroid_partition_source_object_list_covered"])
        self.assertEqual(row["n_source_objects_actually_scanned"],9)


if __name__=="__main__":
    unittest.main()
