import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_niche_temporal_v405 import annual_row,paired,MODES
class TestV405(unittest.TestCase):
    def test_missing_support_is_not_zero(self):
        self.assertIsNone(annual_row({},[49,set(range(20)),set(range(5))],[],MODES[0]))
        self.assertIsNone(annual_row({},[100,set(range(20)),set(range(5))],[],MODES[0]))
    def test_no_model_without_support(self):
        self.assertEqual(paired({y:[] for y in ("2023","2024","2025")},("x",),1)["status"],"insufficient common-year support")
if __name__=="__main__":unittest.main()
