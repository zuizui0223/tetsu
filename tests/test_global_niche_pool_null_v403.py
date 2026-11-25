import sys,unittest,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_niche_pool_null_v403 import expectation,strata,occupancy_bin
class TestV403(unittest.TestCase):
    def test_exact_variance_against_enumeration(self):
        vals={i:np.array([x,0.,w,0.,0.]) for i,(x,w) in enumerate([(1,3),(2,5),(5,7),(8,2),(10,4)])}
        g={0:"a",1:"a",2:"a",3:"b",4:"b"}
        combos=[a+b for a in itertools.combinations([0,1,2],2) for b in itertools.combinations([3,4],1)]
        v,w,m,mov,ng=expectation([0,1,3],g,vals)
        self.assertAlmostEqual(v,np.mean([np.var([vals[t][0] for t in c]) for c in combos]))
        self.assertAlmostEqual(w,np.mean([np.mean([vals[t][2] for t in c]) for c in combos]))
        self.assertAlmostEqual(m,np.mean([np.mean([vals[t][0] for t in c]) for c in combos]))
        self.assertEqual(mov,1.)
    def test_unknown_taxa_not_pooled(self):
        d=strata([1,2,3,4],{1:1,2:1,3:2,4:3},{3:10,4:10},"genus_occupancy")
        self.assertNotEqual(d[1],d[2]);self.assertEqual(d[3],d[4])
        self.assertEqual(occupancy_bin(8),3)
if __name__=="__main__":unittest.main()
