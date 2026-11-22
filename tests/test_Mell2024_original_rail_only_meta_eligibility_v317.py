"""Source Q2 railway-only original table audit tests (no live file in git)."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Mell_2024_railway_meta_original_xlsx_v317 import _col,MELL8,MELL5

class OriginalMetaSourceProtection(unittest.TestCase):
    def test_original_meta_source_files_different_from_species_source(self):
        self.assertEqual(MELL8,"13750_2024_328_MOESM8_ESM.xlsx")
        self.assertNotEqual(MELL5,MELL8)
    def test_true_42_original_columns_xlsx_index(self):
        self.assertEqual(_col("A1"),0)
        self.assertEqual(_col("N45"),13)
        self.assertEqual(_col("AO206"),40)
        self.assertEqual(_col("AP206"),41)
        with self.assertRaises(ValueError):_col("A0")
    def test_no_study_multiplication_in_frozen_result(self):
        import json
        p=Path(__file__).resolve().parents[1]/"data"/"derived"/"tetsu_v317_Mell2024_original_railway_meta_two_studies_summary.json"
        r=json.loads(p.read_text())
        self.assertEqual(r["source_original_Q2_cases"],205)
        self.assertEqual(r["source_original_studies"],47)
        self.assertEqual(sum(r["original_LTI_category_case_counts"].values()),205)
        self.assertEqual(r["original_LTI_category_case_counts"]["Railway"],3)
        self.assertEqual(r["railway_source_distinct_publications"],2)
        self.assertEqual(len(r["railway_cases"]),3)
        self.assertEqual(len(set(x["id"].split("_20")[0] for x in r["railway_cases"])),2)
        self.assertFalse(r["meta_analysis_recommended"])
    def test_independent_response_hard_gate(self):
        import json
        p=Path(__file__).resolve().parents[1]/"data"/"derived"/"tetsu_v317_Mell2024_original_railway_meta_two_studies_summary.json"
        x=json.loads(p.read_text())
        rows=x["railway_cases"]
        self.assertEqual(sum(z["response"]=="Species_richness" for z in rows),2)
        self.assertEqual(sum(z["plant_status"]=="Native" for z in rows),2)
        self.assertEqual(sum(z["country"]=="France" for z in rows),1)
        self.assertEqual(sum(z["country"]=="Australia" for z in rows),2)
        self.assertTrue(any(z.get("continent_typographical_warning") for z in rows))
if __name__=="__main__":unittest.main()
