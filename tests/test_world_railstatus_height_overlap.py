import unittest
from analyze_world_railstatus_height_overlap import select,outcome,SHARE
class TestOverlap(unittest.TestCase):
    def test_empty(self):
        p,q,_=select([],("x",))
        self.assertEqual(p,[])
        self.assertEqual(outcome(p,q)["status"],"insufficient matched support")
    def test_match_outcome_blind(self):
        a={"key":(0,0),"continent":"X","block":1,"active_normal_km":90.,"inactive_normal_km":10.,SHARE:.9,"x":1.,"trait_value":-99}
        c={"key":(0,1),"continent":"X","block":1,"active_normal_km":10.,"inactive_normal_km":90.,SHARE:.1,"x":1.,"trait_value":99}
        p,q,_=select([a,c],("x",))
        self.assertEqual(q["n_pairs"],1)
        self.assertEqual(q["max_abs_SMD"],0.)
        c["trait_value"]=-999
        self.assertEqual(select([a,c],("x",))[1],q)
    def test_caliper(self):
        a={"key":(0,0),"continent":"X","block":1,"active_normal_km":900.,"inactive_normal_km":100.,SHARE:.9,"x":1.}
        c={"key":(0,1),"continent":"X","block":1,"active_normal_km":10.,"inactive_normal_km":90.,SHARE:.1,"x":1.}
        self.assertEqual(select([a,c],("x",))[1]["n_pairs"],0)
if __name__=="__main__":unittest.main()
