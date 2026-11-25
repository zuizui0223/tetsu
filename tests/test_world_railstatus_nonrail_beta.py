import unittest
from unittest.mock import patch
import analyze_world_railstatus_nonrail_beta as t

class TestNonrail(unittest.TestCase):
    def test_classification(self):
        z={c+"_km":0 for c in t.CATS}
        a=dict(z,active_normal_km=1)
        self.assertEqual(t.classify(z,z),"both_four_status_zero")
        self.assertEqual(t.classify(a,a),"both_active")
        self.assertEqual(t.classify(z,a),"mixed_or_other")
        self.assertEqual(t.classify(z,dict(z,inactive_normal_km=2)),"mixed_or_other")
    def test_shared_support(self):
        def r(b):return {"continent":"X","block":b}
        a=[r(1),r(1),r(2),r(3),r(3)]
        b=[r(1),r(1),r(2),r(2),r(4),r(4)]
        x,y=t.shared_support(a,b)
        self.assertEqual([t.key(r) for r in x],[("X",1)]*2)
        self.assertEqual([t.key(r) for r in y],[("X",1)]*2)
    def test_empty(self):
        self.assertEqual(t.coupling([])[0]["status"],"insufficient support or rank")
    def test_cluster_difference(self):
        with patch.object(t,"coupling",side_effect=[({},(2.,{("X",1):.1,("X",2):-.1})),({},(1.,{("X",1):.1,("X",2):-.1}))]):
            q=t.comparison([],[])
        self.assertEqual(q["rail_minus_zero"]["interval95_approx"],[1.,1.])
if __name__=="__main__":unittest.main()
