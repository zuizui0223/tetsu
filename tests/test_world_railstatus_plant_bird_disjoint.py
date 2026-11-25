import unittest
import numpy as np
from analyze_world_railstatus_plant_bird_disjoint import extract,combine,linear,interval
class TestDisjoint(unittest.TestCase):
    def test_counts(self):
        p={"a":[0,0],"b":[1,0],"full_a":[60,10,20],"full_b":[60,10,20],"split_counts_a":[[30,5,10],[30,5,10]],"split_counts_b":[[30,5,10],[30,5,10]],"split_similarity":{}}
        self.assertEqual(len(extract([p])),2)
        p2=dict(p,full_a=[61,10,20])
        with self.assertRaises(ValueError):extract([p,p2])
    def test_paired_covariance(self):
        f={0:(2.,{("X",0):.1,("X",1):-.1}),1:(3.,{("X",0):.1,("X",1):-.1})}
        z=interval(combine(f,{0:1.,1:-1.}))
        self.assertEqual(z["interval95_approx"],[-1.,-1.])
    def test_regression(self):
        rng=np.random.default_rng(5);rows=[]
        for b in range(25):
            for i in range(6):
                x,z=rng.normal(size=2);rows.append({"continent":"X","block":b,"x":x,"z":z,"y":2*x+.5*z+b})
        self.assertAlmostEqual(linear(rows,("z","x"),"y","x")[0],2.,places=8)
if __name__=="__main__":unittest.main()
