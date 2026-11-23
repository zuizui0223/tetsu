#!/usr/bin/env python3
"""Synthetic source validators: never treat publisher error page as source."""
import io,sys,unittest,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from openpyxl import Workbook
from audit_external_railway_global_thermal_sources_v342 import (
    summarize_xlsx,summarize_worldclim_zip,ANCHORS,SUPPLEMENT_URLS,
)

class ExternalSourceGateTests(unittest.TestCase):
    def test_synthetic_xlsx_schema(self):
        b=Workbook()
        sh=b.active
        sh.title="field"
        sh.append(["Species","Location", "Plot"])
        sh.append(["Example species","Sample1", "RAIL"])
        buf=io.BytesIO();b.save(buf)
        s=summarize_xlsx(buf.getvalue())
        self.assertEqual(s["sheets"][0]["columns"],3)
        self.assertEqual(s["sheets"][0]["row_count_incl_header"],2)
        self.assertEqual(s["sheets"][0]["first_row_field_names"][0],"Species")

    def test_html_error_rejected_even_if_network_status_200(self):
        with self.assertRaises(ValueError):
            summarize_xlsx(b"<html><body>blocked</body></html>")

    def test_worldclim_requires_all_real_tiff_labels(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            for i in range(1,20):
                z.writestr(f"wc2.1_10m_bio_{i}.tif",b"FAKE-TIFF")
        v=summarize_worldclim_zip(buf.getvalue())
        self.assertTrue(v["has_annual_mean_temperature_bio1_tif"])
        self.assertEqual(v["bio_tif_count"],19)

    def test_archive_without_all_bio_files_is_rejected(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("wc2.1_10m_bio_1.tif",b"fake")
        with self.assertRaises(ValueError):
            summarize_worldclim_zip(buf.getvalue())

    def test_anchor_list_fixed_source():
        self.assertEqual(ANCHORS,("Cynodon dactylon","Poa annua"))
        self.assertEqual(len(SUPPLEMENT_URLS),2)


if __name__=="__main__":unittest.main()
