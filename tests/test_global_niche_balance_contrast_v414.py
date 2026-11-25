import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_global_niche_balance_contrast_v414 import contrast
class TestV414(unittest.TestCase):
    def test_unsupported_not_zero(self):
        m={"status":"insufficient common-year support"}
        self.assertEqual(contrast([],[],(),m,m,1)["status"],"insufficient common-year support")
    def test_fail_on_different_cells(self):
        with self.assertRaises(ValueError):
            contrast([{"key":(1,1)}],[{"key":(1,2)}],(),{}, {},1)
if __name__=="__main__":unittest.main()
