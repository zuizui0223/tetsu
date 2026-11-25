import unittest
import numpy as np
from analyze_world_railstatus_plant_insect_bird_beta import score_fit,fields
class TestBeta(unittest.TestCase):
    def test_block(self):
        rng=np.random.default_rng(9);r=[]
        for b in range(25):
            for i in range(10):
                x,z=rng.normal(size=2);r.append({"continent":"X","block":b,"x":x,"z":z,"y":b+.4*x+.7*z})
        q=score_fit(r,("x","z"),"y",("x","z"),False)
        self.assertAlmostEqual(q["x"][0],.4,places=8);self.assertAlmostEqual(q["z"][0],.7,places=8)
    def test_endpoint(self):
        rng=np.random.default_rng(21);r=[]
        for block in range(25):
            v={(i+block*5,j):rng.normal() for i in range(4) for j in range(4)}
            for a in sorted(v):
                for delta in ((1,0),(0,1)):
                    b=tuple(x+y for x,y in zip(a,delta))
                    if b not in v:continue
                    x,z=rng.normal(size=2)
                    r.append({"continent":"X","block":block,"a":a,"b":b,"x":x,"z":z,"y":v[a]+v[b]+.4*x+.7*z})
        q=score_fit(r,("x","z"),"y",("x","z"),True)
        self.assertAlmostEqual(q["x"][0],.4,places=8);self.assertAlmostEqual(q["z"][0],.7,places=8)
    def test_columns(self):
        for e in (False,True):
            f=fields(("plants","birds","insects","fungi"),e,True);self.assertEqual(len(f),len(set(f)))
if __name__=="__main__":unittest.main()
