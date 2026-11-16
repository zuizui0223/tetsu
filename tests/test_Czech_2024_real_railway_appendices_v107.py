"""Original Czech 2024 published XLSX metadata audit source/checksum gates."""
import sys
import importlib.util
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Czech_2024_real_railway_appendices_v107 import (
    audit_workbook, audit_sheet, clean_cell
)


@unittest.skipUnless(importlib.util.find_spec("openpyxl") is not None,
                     "openpyxl is installed in dedicated real-source audit workflow")
class CzechRailwaySourceDiscoveryTests(unittest.TestCase):
    def test_no_invalid_workbook_accepted(self):
        with tempfile.TemporaryDirectory() as td:
            f=Path(td)/"fake.xlsx"
            f.write_text("rail data",encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"OOXML"):
                audit_workbook(f)

    def test_source_workbook_sha_and_sheet_names_from_actual_bytes(self):
        from openpyxl import Workbook
        with tempfile.TemporaryDirectory() as td:
            f=Path(td)/"synthetic.xlsx"
            wb=Workbook()
            s=wb.active
            s.title="recent_historical"
            s.append(["taxon","station 1","track zone"])
            s.append(["Plant A","1","tracks"])
            wb.save(f)
            x=audit_workbook(f)
            self.assertEqual(x["original_source_sheet_count"],1)
            self.assertEqual(x["original_source_sheet_schemas"][0]["reported_max_rows_including_possible_formatting"],2)
            self.assertEqual(x["original_source_sheet_schemas"][0]["first_eight_rows_first_55_cells"][1][0],"Plant A")
            self.assertEqual(len(x["original_source_sha256"]),64)
            self.assertTrue(x["NOT_original_republication_or_new_rail_plant_survey"])

    def test_metadata_cell_sanitization_limits_sample(self):
        self.assertEqual(len(clean_cell("a"*300)),130)
        self.assertIsNone(clean_cell(None))


if __name__=="__main__":
    unittest.main()
