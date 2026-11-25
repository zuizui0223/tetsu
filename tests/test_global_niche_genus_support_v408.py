import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_global_niche_genus_support_v408 import support,compare,FEATURES
class TestV408(unittest.TestCase):
    def test_frozen_gates(self):
        self.assertFalse(support({"pool_cells":2,"movable_fraction":1.}))
        self.assertFalse(support({"pool_cells":3,"movable_fraction":.19}))
        self.assertTrue(support({"pool_cells":3,"movable_fraction":.2}))
    def test_zero_variance_no_invented_smd(self):
        r={k:1. for k in FEATURES}
        self.assertTrue(all(v is None for v in compare([r],[r]).values()))
if __name__=="__main__":unittest.main()
