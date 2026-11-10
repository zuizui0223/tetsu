"""Offline QA of anonymous GBIF public Parquet preview and reported missing-file snapshot."""
from __future__ import annotations
import unittest
from datetime import datetime,timezone

from botocore.exceptions import ClientError
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from probe_public_GBIF_partitioned_original_v59 import (
    run,list_public_source_prefix,validate_object_sample,PREVIEW_PATHS,
    SOURCE_PUBLISHED_DWCA_DOI
)


class MockPublicGBIFS3:
    def __init__(self,preview=True,missing_20=True):
        self.preview=preview;self.missing_20=missing_20
        self.calls=[]

    def list_objects_v2(self,**kwargs):
        self.calls.append(("list",kwargs))
        if not self.preview:
            return {"CommonPrefixes":[],"Contents":[],"IsTruncated":False}
        return {
            "CommonPrefixes":[
                {"Prefix":kwargs["Prefix"]+"taxonPartition=Plantae_Tracheophyta/"},
                {"Prefix":kwargs["Prefix"]+"taxonPartition=Animalia_Mammalia/"},
            ],
            "Contents":[],
            "IsTruncated":False,
        }

    def head_object(self,**kwargs):
        self.calls.append(("head",kwargs))
        if kwargs["Key"].endswith("/000020") and self.missing_20:
            raise ClientError(
                {"Error":{"Code":"404","Message":"Not Found"},"ResponseMetadata":{"HTTPStatusCode":404}},
                "HeadObject")
        return {"ContentLength":10000,"ContentType":"application/octet-stream"}


class PublicGBIFOriginalParquetPreflightTests(unittest.TestCase):
    def test_anonymous_source_lists_only_bounded_prefixes(self):
        x=MockPublicGBIFS3()
        v=list_public_source_prefix(x,"gbif-public-data",
                "development/2026-08-01-DWCA/occurrence/p_taxon_a5/")
        self.assertEqual(v["n_partition_prefixes_in_bounded_listing"],2)
        self.assertTrue(v["does_NOT_certify_whole_snapshot_completeness"])
        self.assertEqual(x.calls[0][1]["MaxKeys"],100)
        self.assertEqual(x.calls[0][1]["Delimiter"],"/")

    def test_known_missing_monthly_chunk_is_reported_as_missing_not_zero_species(self):
        x=MockPublicGBIFS3()
        ok=validate_object_sample(x,"gbif-open-data-us-east-1",
                "occurrence/2026-10-01/occurrence.parquet/000019",True)
        missing=validate_object_sample(x,"gbif-open-data-us-east-1",
                "occurrence/2026-10-01/occurrence.parquet/000020",False)
        self.assertTrue(ok["present"])
        self.assertFalse(missing["present"])
        self.assertTrue(missing["matches_expected_presence"])

    def test_real_data_NOT_claimed_merely_from_public_s3_metadata(self):
        mock=MockPublicGBIFS3()
        audit=run(mock,mock)
        self.assertTrue(audit["complete_probe"])
        self.assertTrue(audit["legacy_reproduces_published_missing_chunk_issue"])
        self.assertEqual(len(audit["public_preview_source_prefixes"]),2)
        self.assertTrue(audit["no_original_occurrence_rows_claimed_or_downloaded"])
        self.assertTrue(audit["no_complete_GBIF_event_population_claimed"])
        self.assertEqual(audit["preview_underlying_source_GBIF_DWCA_DOI"],SOURCE_PUBLISHED_DWCA_DOI)

    def test_missing_s3_listing_fails_closed(self):
        mock=MockPublicGBIFS3(preview=False)
        audit=run(mock,mock)
        self.assertFalse(audit["complete_probe"])
        self.assertFalse(audit["preview_list_source_accessible"])
        self.assertIn("preview_"+PREVIEW_PATHS[0],audit["errors"])

    def test_existing_missing_chunk_diagnostic_must_not_fake_missingness(self):
        mock=MockPublicGBIFS3(missing_20=False)
        audit=run(mock,mock)
        self.assertFalse(audit["legacy_reproduces_published_missing_chunk_issue"])
        self.assertTrue(audit["complete_probe"])


if __name__=="__main__":
    unittest.main()
