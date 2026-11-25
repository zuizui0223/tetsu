import unittest
import numpy as np
from analyze_world_railstatus_woodiness_coverage_audit import split_common,compare,ACTIVE,INACTIVE
class TestCoverageAudit(unittest.TestCase):
    def test_common_blocks(self):
        r=[{"continent":"X","block":0,"coverage":c} for c in (.3,.4,.5,.6)]
        r+=[{"continent":"X","block":1,"coverage":c} for c in (.3,.4,.6)]
        k,l,h=split_common(r);self.assertEqual(k,[("X",0)]);self.assertEqual((len(l),len(h)),(2,2))
    def test_direct_paired_difference(self):
        rng=np.random.default_rng(8);r=[]
        for k in range(30):
            for hi in (False,True):
                for i in range(5):
                    a,b=rng.uniform(0,1,2)
                    expected=.25+.05*a+.02*b
                    excess=(.01 if hi else -.1)*a+.02*b
                    r.append({"continent":"X","block":k,ACTIVE:a,INACTIVE:b,
                       "coverage":.6 if hi else .4,"observations":100,"species":50,
                       "active_normal_km":a,"inactive_normal_km":b,
                       "observed":expected+excess,"expected":expected,"excess":excess})
        q=compare(r,(ACTIVE,INACTIVE))
        self.assertEqual(q["status"],"estimated")
        self.assertAlmostEqual(q["effects"]["active_minus_inactive"]["low_minus_high"]["excess"]["estimate"],-.11,places=9)
if __name__=="__main__":unittest.main()
