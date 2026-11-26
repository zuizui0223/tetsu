import unittest
import numpy as np
import analyze_world_railstatus_selection_weighted_coupling as t
class TestSelection(unittest.TestCase):
    def test_weighted_regression(self):
        rng=np.random.default_rng(4);rows=[]
        for g in range(25):
            for j in range(8):
                x,z=rng.normal(size=2);rows.append({"continent":"X","block":g,"key":(g,j),"x":x,"z":z,"y":2*x+z+g})
        w=np.ones(len(rows))/len(rows)
        q=t.weighted_fit(rows,("z","x"),"y","x",w)
        self.assertAlmostEqual(q[0],2.,places=8)
    def test_native_weight_replay(self):
        rng=np.random.default_rng(9);rows=[]
        for g in range(25):
            for j in range(8):
                x,z,e=rng.normal(size=3);rows.append({"continent":"X","block":g,"key":(g,j),"x":x,"z":z,"y":2*x+z+g+e})
        q=t.weighted_fit(rows,("z","x"),"y","x",t.weights(rows,"native",{},{}))
        z=t.linear(rows,("z","x"),"y","x")
        t.check_effect(t.interval(q),t.interval(z))
    def test_spatial_fold(self):
        r={"continent":"X","block":(0,0),"key":(0,1)}
        self.assertEqual(t.fold(r),t.fold(dict(r,key=(0,2))))
    def test_weights(self):
        rows=[{"continent":"X","block":0,"key":0},{"continent":"X","block":0,"key":1}]
        p={0:.2,1:.4}
        w=t.weights(rows,"ipw",{("X",0):10},p)
        np.testing.assert_allclose(w,[2/3,1/3])
if __name__=="__main__":unittest.main()
