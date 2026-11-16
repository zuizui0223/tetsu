"""v11.1 independent Dryad DOI and source metadata protection tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Kashmir_2021_rail_Dryad_source_v111 import (
    NLM_CATALOG_DATA_DOI,PAPER_REFERENCE_TYPO_CANDIDATE,
    url_for_dataset,summarize_file_list,audit
)


class KashmirDryadOriginalSourceTests(unittest.TestCase):
    def test_published_catalog_DOI_differs_from_article_reference_typo(self):
        self.assertEqual(NLM_CATALOG_DATA_DOI,"10.5061/dryad.0p2ngf20h")
        self.assertNotEqual(NLM_CATALOG_DATA_DOI,PAPER_REFERENCE_TYPO_CANDIDATE)
        self.assertIn("doi%3A10.5061%2Fdryad.0p2ngf20h",
                      url_for_dataset(NLM_CATALOG_DATA_DOI))

    def test_publisher_fails_closed_for_no_source(self):
        def bad(u):raise OSError("No real public source available")
        r=audit(bad)
        self.assertFalse(r["complete"])
        self.assertFalse(r["source_data_raw_downloaded_or_ecology_fitted"])
        self.assertEqual(r["status"],"original_Dryad_metadata_not_verified_no_biological_result")

    def test_title_DOI_mismatch_is_not_a_verified_railway_dataset(self):
        def wrong(u):return {"title":"Unrelated road ecology","identifier":
                              "doi:10.5061/dryad.0p2ngf20h"}
        r=audit(wrong)
        self.assertFalse(r["complete"])
        self.assertEqual(r["status"],"Dryad_source_title_mismatch_no_biological_result")

    def test_files_are_from_original_Dryad_publisher_metadata_only(self):
        source={
            "https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.0p2ngf20h":{
                "id":10,"title":"Species abundance along the railway of Kashmir Himalaya",
                "identifier":"doi:10.5061/dryad.0p2ngf20h",
                "_links":{"stash:versions":{"href":"/api/v2/datasets/foo/versions"}}},
            "https://datadryad.org/api/v2/datasets/foo/versions":{
                "_embedded":{"stash:versions":[{"id":22,"versionNumber":1,"versionStatus":"published"}]}},
            "https://datadryad.org/api/v2/versions/22/files":{
                "_embedded":{"stash:files":[{"id":3,"path":"plant-plot-survey.csv","size":1024,
                     "digest":"abc","digestType":"md5"}]}}
        }
        r=audit(lambda u:source[u])
        self.assertTrue(r["complete"])
        self.assertEqual(r["n_original_metadata_files"],1)
        self.assertEqual(r["original_publisher_dataset_files"][0]["path"],
                         "plant-plot-survey.csv")
        self.assertFalse(r["source_data_raw_downloaded_or_ecology_fitted"])

    def test_raw_data_are_not_implicitly_present(self):
        x=summarize_file_list({"_embedded":{"stash:files":[{"id":42,"path":"metadata.csv"}]}})
        self.assertEqual(len(x),1)
        self.assertIsNone(x[0]["size"])


if __name__=="__main__":unittest.main()
