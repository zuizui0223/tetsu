"""v14.7 source-safe OOXML parser tests; only SYNTHETIC cells in this file."""
import tempfile
import sys
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Mell2024_supplementary_railway_meta_v147 import (
    RAIL,examine,xlsx_cell_index,sheet_files,worksheet_rows,get_shared_strings
)

M = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
P = "http://schemas.openxmlformats.org/package/2006/relationships"


def make_fake_original_sheet(path, invalid_index=False):
    sheet=f"""<?xml version="1.0" encoding="UTF-8"?>
<worksheet xmlns="{M}"><sheetData>
  <row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="inlineStr"><is><t>Effect</t></is></c></row>
  <row r="2"><c r="A2" t="s"><v>{99 if invalid_index else 1}</v></c><c r="B2"><v>0.125</v></c></row>
  <row r="3"><c r="A3" t="s"><v>2</v></c><c r="B3"><v>-1.25</v></c></row>
</sheetData></worksheet>"""
    with zipfile.ZipFile(path,"w",zipfile.ZIP_STORED) as z:
        z.writestr("xl/workbook.xml",
            f'<workbook xmlns="{M}" xmlns:r="{R}"><sheets>'
            f'<sheet name="Source meta" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels",
            f'<Relationships xmlns="{P}"><Relationship Id="rId1" '
            f'Type="{R}/worksheet" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml",
            f'<sst xmlns="{M}"><si><t>Transport mode</t></si>'
            f'<si><t>Railways</t></si><si><t>Trail vegetation</t></si></sst>')
        z.writestr("xl/worksheets/sheet1.xml",sheet)
        z.writestr("docProps/custom.xml","a"*3000)


class TestOriginalMetaDiagnostic(unittest.TestCase):
    def test_no_trail_false_positive_and_rail_match(self):
        self.assertTrue(RAIL.search("Railways"))
        self.assertTrue(RAIL.search("rail tracks"))
        self.assertFalse(RAIL.search("Trail vegetation"))
        self.assertFalse(RAIL.search("railwaymanuscript"))  # avoid substrings

    def test_excel_cell_index_bounds(self):
        self.assertEqual(xlsx_cell_index("A2"),0)
        self.assertEqual(xlsx_cell_index("AA31"),26)
        with self.assertRaises(ValueError):
            xlsx_cell_index("2AAA")
        with self.assertRaises(ValueError):
            xlsx_cell_index("ZZ999")

    def test_real_oo_xml_parser_preserves_row_ids_shared_and_numeric(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"fictional.xlsx"
            make_fake_original_sheet(f)
            with zipfile.ZipFile(f) as z:
                sh=sheet_files(z)
                self.assertEqual(sh,[("Source meta","xl/worksheets/sheet1.xml")])
                values=list(worksheet_rows(z,sh[0][1],get_shared_strings(z)))
                self.assertEqual(values[1],(2,{0:"Railways",1:"0.125"}))
                self.assertEqual(values[2],(3,{0:"Trail vegetation",1:"-1.25"}))
            result=examine(f,"effect_cases_8")
            self.assertEqual(result["source_worksheet_count"],1)
            source=result["source_sheets"][0]
            self.assertEqual(source["raw_rail_keyword_row_candidates_shown"],1)
            self.assertEqual(source["candidate_rail_rows_for_manual_column_validation"][0]["source_excel_row"],2)

    def test_malformed_shared_string_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"fake.xlsx"
            make_fake_original_sheet(f,invalid_index=True)
            with self.assertRaisesRegex(ValueError,"shared-string index"):
                examine(f,"study_metadata_5")

    def test_text_file_cannot_be_called_original_supplement(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"no.xlsx"
            f.write_text("not a source XLSX",encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"publisher downloadable"):
                examine(f,"effect_cases_8")


if __name__=="__main__":unittest.main()
