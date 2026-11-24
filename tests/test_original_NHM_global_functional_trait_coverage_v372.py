"""v372 synthetic tests of official NHM 59k source trait availability semantics."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_original_NHM_global_functional_trait_coverage_v372 as m

class NHMGlobalTraitCoverage(unittest.TestCase):
    def test_missing_source_trait_not_assigned_other(self):
        for x in ("","NA","unknown","?","not reported","n/a",None):
            self.assertFalse(m.source_value(x))
        self.assertTrue(m.source_value("annual herb"))
        self.assertTrue(m.source_value("wind"))
    def test_original_plant_zero_numeric_not_missing_by_default(self):
        self.assertTrue(m.source_value(0))
        self.assertTrue(m.source_value("0"))
    def test_exact_original_publisher_source_resource(self):
        self.assertEqual(m.EXPECTED_TOTAL,59369)
        self.assertEqual(m.SOURCE_RESOURCE,"2552d157-9519-476c-995d-906e2433cc19")
    def test_full_source_requires_all_declared_functional_fields(self):
        self.assertIn("dispersal mode",m.COLUMNS)
        self.assertIn("life form",m.COLUMNS)
        self.assertIn("habit",m.COLUMNS)
        self.assertIn("plant max_ height [m]",m.COLUMNS)
        self.assertIn("seed min_ diameter [mm]",m.COLUMNS)
    def test_old_nonstandard_unmapped_categories_not_imputed(self):
        self.assertTrue(m.source_value("climbing shrub"))
        self.assertTrue(m.source_value("anemochorous"))
        self.assertNotEqual(m.source_value("wind"),"anemochorous")
    def test_source_page_limit(self):
        self.assertTrue(0<m.PAGE<m.EXPECTED_TOTAL)
        self.assertEqual(m.DOI,"10.5519/p3dm31kc")

if __name__=="__main__":unittest.main()
