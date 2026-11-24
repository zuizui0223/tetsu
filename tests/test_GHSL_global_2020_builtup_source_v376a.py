"""v376a authenticated official JRC global source HEAD and ZIP prefix tests."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_GHSL_global_2020_builtup_source_v376a as m
class OriginalGlobalBuiltupSource(unittest.TestCase):
    def test_exact_2020_original_1km_product(self):
        self.assertIn("GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000",m.SOURCE_URL)
        self.assertTrue(m.SOURCE_URL.endswith(".zip"))
    def test_2020_is_time_before_2025_plants(self):
        self.assertIn("E2020",m.SOURCE_URL)
        self.assertNotIn("E2030",m.SOURCE_URL)
    def test_head_valid_content_length(self):
        x=m.parse_head(200,{"Content-Length":"152527564","ETag":"123"})
        self.assertEqual(x["source_original_content_length_bytes"],152527564)
        self.assertEqual(x["source_original_etag"],"123")
    def test_head_unreadable_length_blocks(self):
        with self.assertRaises(ValueError):
            m.parse_head(200,{"Content-Length":"non-numeric"})
    def test_missing_length_stays_unobserved(self):
        x=m.parse_head(200,{})
        self.assertIsNone(x["source_original_content_length_bytes"])
    def test_actual_zip_signature_required(self):
        self.assertEqual(m.ZIP_MAGIC,b"PK"+bytes([3,4]))
        self.assertFalse(b"<html>".startswith(m.ZIP_MAGIC))
    def test_official_european_commission_source(self):
        self.assertTrue(m.SOURCE_URL.startswith("https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/"))
    def test_dataset_citation(self):
        self.assertEqual(m.DOI,"10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA")
if __name__=="__main__":unittest.main()
