import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_turnover_effort_v388 import close_ratio,support
class TestV388(unittest.TestCase):
    def test_support_frozen_and_symmetric(self):
        a=dict(observations=200,observers=20,species=40)
        b=dict(observations=400,observers=40,species=80)
        for mode in ("all","effort_balance","balance_richness","high_effort_balance"):
            self.assertTrue(support(a,b,mode));self.assertTrue(support(b,a,mode))
        c=dict(b,species=81)
        self.assertFalse(support(a,c,"balance_richness"))
        self.assertTrue(support(a,c,"effort_balance"))
        self.assertFalse(support(a,dict(b,observations=401),"effort_balance"))
        self.assertFalse(support(dict(a,observers=19),b,"high_effort_balance"))
    def test_zero_support(self):
        self.assertFalse(close_ratio(0,1,2.))
if __name__=="__main__":unittest.main()
