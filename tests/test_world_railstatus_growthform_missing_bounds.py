import unittest,itertools,math
import numpy as np
from analyze_world_railstatus_growthform_missing_bounds import gini_bounds,woodiness_bounds,bound,linear_weights,ACTIVE,INACTIVE
class TestMissingBounds(unittest.TestCase):
    def test_exact_integer_extrema(self):
        counts={"F":3,"C":1};N=7
        vals=[]
        for extra in itertools.product(range(7),repeat=3):
            c=[3,1,0,0,0,0,0]
            for i in extra:c[i]+=1
            vals.append((N*N-sum(x*x for x in c))/(N*(N-1)))
        self.assertEqual(gini_bounds(counts,N),(min(vals),max(vals)))
    def test_complete_and_wood(self):
        l,h=gini_bounds({"F":10,"C":10},20)
        self.assertAlmostEqual(l,200/380);self.assertEqual(l,h)
        self.assertEqual(woodiness_bounds(3,8,10),(.3,.5))
    def test_signed_bound(self):
        self.assertEqual(bound([2,-1],[0,2],[1,3]),[-3.,0.])
    def test_weights_replay(self):
        rng=np.random.default_rng(9);rows=[]
        for b in range(25):
            for i in range(6):
                a,z=rng.uniform(0,3,2)
                rows.append({"continent":"X","block":b,ACTIVE:a,INACTIVE:z,"y":2*a-.3*z+b})
        w=linear_weights(rows,(ACTIVE,INACTIVE))
        y=np.array([r["y"] for r in rows])
        self.assertAlmostEqual(float(w["active_minus_inactive"]@y),2.3,places=8)
if __name__=="__main__":unittest.main()
