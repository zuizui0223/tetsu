import unittest
import numpy as np
import analyze_world_railstatus_joint_climate_selection as t
class TestJoint(unittest.TestCase):
    def test_nonlinear_fields(self):
        self.assertTrue(set(t.SQUARE_BASE).isdisjoint(t.LOG.values()))
        self.assertNotIn("bio1",t.SQUARE_BASE)
        self.assertNotIn("log_observers",t.SQUARE_BASE)
        self.assertNotIn("log_observations",t.SQUARE_BASE)
    def test_weighted_fit(self):
        rng=np.random.default_rng(9);rows=[]
        for g in range(25):
            for j in range(8):
                x,z,c,e=rng.normal(size=4);rows.append({"continent":"X","block":g,"key":(g,j),"x":x,"z":z,"c":c,"y":2*x+z+c+g+e})
        q=t.weighted_fit(rows,("z","c","x"),"y","x",t.weights(rows,"native",{},{}))
        t.check_effect(t.interval(q),t.interval(t.linear(rows,("z","c","x"),"y","x")))
    def test_effort_status_support(self):
        rng=np.random.default_rng(2);rr=[]
        for g in range(30):
            for j in range(7):
                a,b,z=rng.normal(size=3);rr.append({"continent":"X","block":g,t.b.ACTIVE:a,t.b.INACTIVE:b,"z":z,"both_effort":a+2*b+z+g})
        r,_=t.effort_fit(rr,(t.b.ACTIVE,t.b.INACTIVE,"z"))
        self.assertAlmostEqual(r["effects"]["active_minus_inactive"]["estimate"],-1.,places=8)
if __name__=="__main__":unittest.main()
