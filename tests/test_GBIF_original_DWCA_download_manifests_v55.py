"""GBIF authenticated DOI download request-spec generator: source-only, no credential use."""
import sys
from pathlib import Path
import json
import tempfile
import unittest
from shapely import wkt

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from prepare_GBIF_original_DWCA_download_manifests_v55 import (
    predicate_for,prepare,DOWNLOAD_URL,PILOT,MAX_WITHIN_SOURCE_POINTS
)


class OfficialOriginalDOIManifestTests(unittest.TestCase):
    def test_all_10_pre_frozen_nonoverlapping_country_block_geometry(self):
        self.assertEqual(len(PILOT),10)
        for block in PILOT:
            req,n=predicate_for(block)
            self.assertEqual(req["format"],"DWCA")
            self.assertEqual(req["sendNotification"],False)
            self.assertLess(n,MAX_WITHIN_SOURCE_POINTS)
            self.assertEqual(req["predicate"]["type"],"and")
            p=req["predicate"]["predicates"]
            geo=[x for x in p if x["type"]=="or"][0]["predicates"]
            self.assertEqual(len(geo),16)
            self.assertTrue(all(x["type"]=="within" for x in geo))
            self.assertTrue(all(wkt.loads(x["geometry"]).exterior.is_ccw for x in geo))

    def test_no_legacy_present_non_detected_status_predicate_or_sample_bias(self):
        req,_=predicate_for("us_southeast")
        p=req["predicate"]["predicates"]
        self.assertEqual({v["key"] for v in p if "key" in v},
            {"COUNTRY","KINGDOM_KEY","TAXON_RANK","HAS_COORDINATE",
             "HAS_GEOSPATIAL_ISSUE","YEAR"})
        self.assertNotIn("occurrenceStatus",json.dumps(req))
        years=[x for x in p if x.get("key")=="YEAR"]
        self.assertEqual({x["type"] for x in years},
                         {"greaterThanOrEquals","lessThanOrEquals"})
        self.assertEqual({x["value"] for x in years},{"2010","2025"})

    def test_auth_credentials_never_needed_or_stored_to_prepare(self):
        with tempfile.TemporaryDirectory() as temp:
            out=prepare(temp)
            self.assertTrue(out["NO_GBIF_download_REQUEST_posted"])
            self.assertTrue(out["NO_DOI_CREATED_by_these_request_manifest_files"])
            self.assertEqual(out["GBIF_official_download_API"],DOWNLOAD_URL)
            self.assertEqual(len(out["source_region_download_specs"]),10)
            files=list(Path(temp).glob("*_official_DWCA_request.json"))
            self.assertEqual(len(files),10)
            self.assertTrue((Path(temp)/"manifest_index.json").exists())
            alltext=" ".join(path.read_text() for path in files)
            self.assertNotIn("username",alltext.lower())
            self.assertNotIn("password",alltext.lower())
            self.assertNotIn("notificationAddresses",alltext)
            self.assertNotIn("occurrenceStatus",alltext)

    def test_unselected_region_is_not_introduced_after_hypothesis_sign(self):
        with self.assertRaisesRegex(ValueError,"original frozen"):
            predicate_for("jp_sado")
        with self.assertRaises(ValueError):
            predicate_for("fr_mediterranean")


if __name__=="__main__":
    unittest.main()
