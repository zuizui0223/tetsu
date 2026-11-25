import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_niche_common_support_v400 import prepare,accepted,MODIFIERS,TARGET
class TestV400(unittest.TestCase):
    def test_symmetric_modifier_and_coverage(self):
        r={k:2. for k in MODIFIERS}
        r.update(key=[0,0],block=[0,0],bio1=5.,trait_coverage=.6,known_niche_species=20)
        r[TARGET]=3.
        a=prepare(r);self.assertEqual(a["rail_cold"],3.)
        self.assertAlmostEqual(a["coverage_cold"],.6)
        self.assertTrue(all(a[k+"_cold"]==2 for k in MODIFIERS))
        self.assertTrue(accepted(a))
        a["trait_coverage"]=.49;self.assertFalse(accepted(a))
if __name__=="__main__":unittest.main()
