import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_niche_rail_swap_v404 import pairs,swapped,FEATURES
class TestV404(unittest.TestCase):
    def test_match_only_same_block(self):
        rows=[]
        for i in range(4):
            r={k:1. for k in FEATURES};r.update(continent="Asia",block=(i//2,0))
            rows.append(r)
        self.assertEqual(pairs(rows),[(0,1),(2,3)])
    def test_recompute_interaction(self):
        x=np.array([[1.,2.],[3.,12.]])
        y=swapped(x,[(0,1)],np.random.default_rng(0),0,1,np.array([2.,4.]))
        np.testing.assert_array_equal(np.sort(y[:,0]),[1,3])
        np.testing.assert_array_equal(y[:,1],y[:,0]*[2,4])
        np.testing.assert_array_equal(x,[[1,2],[3,12]])
if __name__=="__main__":unittest.main()
