"""v11.2 original Dryad file SHA and spreadsheet provenance test."""
import sys,unittest,tempfile,importlib.util
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Kashmir_2021_original_plot_matrix_v112 import (
  SOURCE_SHA256,SOURCE_N_BYTES,SOURCE_FILE_URL,verify_original_file,
)
class KashmirOriginalMatrixTest(unittest.TestCase):
    def test_locked_original_metadata_is_exact(self):
        self.assertEqual(len(SOURCE_SHA256),64)
        self.assertEqual(SOURCE_N_BYTES,122187)
        self.assertTrue(SOURCE_FILE_URL.endswith("/files/546003/download"))
    def test_not_a_spoofed_published_original_workbook(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"Matrix_plot.xlsx"
            f.write_bytes(b"UNVERIFIED TOY DATA")
            with self.assertRaisesRegex(ValueError,"SHA/bytes differ"):
                verify_original_file(f)
if __name__=="__main__":unittest.main()
