"""Original Slovak2016 19-station PDF preflight synthetic math QA."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
 import fitz
 from audit_Slovakia2016_original_19station_318taxa_pdf_preflight_v310 import preflight,fold
except ModuleNotFoundError as e:
 if e.name!="fitz":raise
 fitz=None
@unittest.skipUnless(fitz is not None,"Source original PDF CI installs PyMuPDF only")
class SourcePDFGate(unittest.TestCase):
 def test_diacritic_fold_and_original_source_semantics(self):
  self.assertEqual(fold("Zoznam zistených taxónov"),"zoznam zistenych taxonov")
 def test_html_not_pdf_fails_closed(self):
  import tempfile
  with tempfile.TemporaryDirectory() as d:
   x=Path(d)/"fake.pdf";x.write_bytes(b"<html>"+b"x"*100000)
   with self.assertRaises(ValueError):preflight(x)
 def test_wrong_number_pages_fails_closed(self):
  import tempfile
  with tempfile.TemporaryDirectory() as d:
   x=Path(d)/"not-original.pdf"
   doc=fitz.open()
   doc.new_page().insert_text((20,20),"Poprad Plavec Zoznam zistenych taxonov")
   doc.save(str(x));doc.close()
   with self.assertRaises(ValueError):preflight(x)
if __name__=="__main__":unittest.main()
