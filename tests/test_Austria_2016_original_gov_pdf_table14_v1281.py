"""v12.8.1 source PDF exact per-site author numeric extraction and no fake seed data."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_Austria_ferreral_osr_actual_sites_v128 import validate_original_table
from audit_Austria_2016_original_gov_pdf_table14_v1281 import (
    isolate_table14,extract_original_pdf_table_pairs,ORIGINAL_SITE_KEYS,
    ORIGINAL_GOV_PDF_SHA256,
)

CSV=(Path(__file__).resolve().parents[1]/"data"/
     "austria_pascher_2016_feral_osr_table14_manual_transcription_v128.csv")


class Original2016GovernmentReportPDFExactPairs(unittest.TestCase):
    def test_all_published_60_site_120_numeric_cells_recoverable_from_table_layout(self):
        rows=validate_original_table(CSV)
        prefix="Austrian official original 2016 report\n"
        body="Table 14: Sample sites and number of sampled feral OSR plants.\n"
        for r in rows:
            body+=(f"{r['source_site_code']} "
                   f"Federal State test "
                   f"{r['sampled_2014']:>4}  {r['sampled_2015']:>4}\n")
        body+="In total, 2,113 feral plant individuals were sampled in 2014 and 2015.\n"
        text=prefix+body
        result=extract_original_pdf_table_pairs(isolate_table14(text))
        self.assertEqual(len(result),60)
        self.assertEqual(set(result),ORIGINAL_SITE_KEYS)
        self.assertEqual(result["BHR1"],(26,47))
        self.assertEqual(result["GBH7"],(1,2))
        self.assertEqual(result["SK11"],(0,0))
        self.assertEqual(result["OEM3"],(165,60))

    def test_manual_table_site_transposition_not_allowed_although_totals_still_match(self):
        rows=validate_original_table(CSV)
        text="Table 14: Sample sites and number of sampled feral OSR plants.\n"
        for i,r in enumerate(rows):
            num1=r["sampled_2014"]
            num2=r["sampled_2015"]
            if r["source_site_code"]=="BHR1":num1+=1
            if r["source_site_code"]=="BHR2":num1-=1
            text+=f"{r['source_site_code']} Source country test {num1} {num2}\n"
        text+="In total, 2,113 feral plant individuals were sampled."
        pdf=extract_original_pdf_table_pairs(isolate_table14(text))
        self.assertEqual(sum(x[0] for x in pdf.values()),1001)
        self.assertNotEqual(pdf["BHR1"],(26,47))
        self.assertNotEqual(pdf["BHR2"],(11,3))

    def test_source_audit_refuses_incomplete_original_pdf_table(self):
        rows=validate_original_table(CSV)
        txt="Table 14: Sample sites and number of sampled feral OSR plants.\n"
        for r in rows[:59]:
            txt+=f"{r['source_site_code']} Original Austrian field report location Federal State original source sampling {r['sampled_2014']} {r['sampled_2015']}\n"
        txt+="In total, 2,113 feral plant individuals were sampled."
        with self.assertRaisesRegex(ValueError,"Cannot faithfully recover all"):
            extract_original_pdf_table_pairs(isolate_table14(txt))

    def test_report_checksum_is_genuine_full_length_sha256(self):
        self.assertEqual(len(ORIGINAL_GOV_PDF_SHA256),64)
        self.assertEqual(ORIGINAL_GOV_PDF_SHA256,
          "52a1910fb24cc9a38da841322ce28e9c97b51758a632146551981edac696125d")


if __name__=="__main__":unittest.main()
