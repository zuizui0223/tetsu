import unittest
import numpy as np
from analyze_world_railstatus_cross_taxa_population_local import fit,paired,ACTIVE,INACTIVE,POP
class TestJoint(unittest.TestCase):
    def test_direct_double_contrast(self):
        rng=np.random.default_rng(18);a=[];b=[]
        for k in range(30):
            for j in range(6):
                x,z,p=rng.uniform(0,5,3)
                r={"continent":"X","block":k,ACTIVE:x,INACTIVE:z,POP:p}
                a.append(dict(r,log_species=.4*x-.1*z+.2*p+k))
                b.append(dict(r,log_species=.2*x+.05*z+.4*p-k))
        aa,sa=fit(a,(ACTIVE,INACTIVE,POP),True);bb,sb=fit(b,(ACTIVE,INACTIVE,POP),True)
        self.assertEqual(aa["status"],"estimated")
        d=paired(sa,sb)
        self.assertAlmostEqual(d["active"]["estimate"],.2,places=8)
        self.assertAlmostEqual(d["active_minus_inactive"]["estimate"],.35,places=8)
    def test_no_support(self):
        z,s=fit([], (ACTIVE,INACTIVE,POP), True)
        self.assertEqual(z["status"],"insufficient support")
        self.assertIsNone(paired(s,s))
if __name__=="__main__":unittest.main()
