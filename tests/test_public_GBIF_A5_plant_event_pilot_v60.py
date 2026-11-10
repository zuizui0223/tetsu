"""A5-partitioned public GBIF source is NOT silently reported as complete 50km plant survey."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from probe_public_GBIF_A5_plant_events_v60 import (
    ROOT,BUCKET,FROZEN_BLOCKS,equal_area_source_cells,
    listing,find_small_geographic_plant_files
)


class FakeS3:
    def __init__(self):
        self.calls=[]
    def list_objects_v2(self,**args):
        self.calls.append(args)
        prefix=args["Prefix"]
        if prefix==ROOT:
            return {"CommonPrefixes":[
                {"Prefix":ROOT+"taxonpartition=Animalia/"},
                {"Prefix":ROOT+"taxonpartition=Plantae_Tracheophyta/"},
            ],"IsTruncated":False}
        if "a5_r2=" in prefix and prefix.endswith("/123/"):
            return {"Contents":[{"Key":prefix+"actual.parquet","Size":12345678}],
                    "IsTruncated":False}
        return {"Contents":[],"IsTruncated":False}


class OriginalA5PilotTests(unittest.TestCase):
    def test_preserve_geographic_50km_source_blocks_and_16_cells(self):
        for block in FROZEN_BLOCKS:
            cells=equal_area_source_cells(block)
            self.assertEqual(len(cells),16)
            self.assertTrue(all(x["gbif_geometry_wgs84"] for x in cells))
        with self.assertRaises(ValueError):
            equal_area_source_cells("unfrozen")

    def test_source_only_one_parquet_object_is_not_a_full_archive(self):
        s3=FakeS3()
        parts,files=find_small_geographic_plant_files(s3,[123,456])
        self.assertEqual(len(parts),1)
        self.assertEqual(len(files),1)
        self.assertEqual(files[0]["a5_r2"],123)
        self.assertEqual(files[0]["size_bytes"],12345678)
        self.assertTrue(files[0]["key"].startswith(ROOT+"taxonpartition=Plantae_"))
        self.assertNotIn("Animalia",files[0]["key"])
        self.assertEqual(s3.calls[0]["Bucket"],BUCKET)

    def test_list_metadata_does_not_hide_bounded_truncation(self):
        s3=FakeS3()
        folders,files,truncated=listing(s3,ROOT)
        self.assertEqual(len(folders),2)
        self.assertEqual(len(files),0)
        self.assertFalse(truncated)

    def test_missing_predeclared_plants_fails_closed(self):
        class NoPlants(FakeS3):
            def list_objects_v2(self,**params):
                return {"CommonPrefixes":[{"Prefix":ROOT+"taxonpartition=Animalia/"}],
                        "IsTruncated":False}
        with self.assertRaisesRegex(ValueError,"No anonymous public GBIF Plantae"):
            find_small_geographic_plant_files(NoPlants(),[123])


if __name__=="__main__":
    unittest.main()
