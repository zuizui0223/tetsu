import unittest
import numpy as np
from analyze_world_railstatus_niche_pool import fit,ACTIVE,INACTIVE,TAGS
class TestStatusNiche(unittest.TestCase):
    def test_paired_interaction_and_decomposition(self):
        rng=np.random.default_rng(301);rows=[]
        for k in range(35):
            for i in range(5):
                a,b,c=rng.uniform(0,3,3)
                r={"continent":"X","block":k,"bio1":15-10*c,ACTIVE:a,INACTIVE:b,
                   ACTIVE+"_cold":a*c,INACTIVE+"_cold":b*c}
                for t in TAGS:
                    expected=2*a+.1*b+.2*a*c+.3*b*c+k
                    excess=-.3*a*c+.1*b*c
                    r.update({"expected_"+t:expected,"excess_"+t:excess,"observed_"+t:expected+excess})
                rows.append(r)
        z=fit(rows,(ACTIVE,INACTIVE,ACTIVE+"_cold",INACTIVE+"_cold"))
        self.assertEqual(z["status"],"estimated")
        self.assertAlmostEqual(z["effects"]["active_minus_inactive_cold"]["excess_T_variance"]["estimate"],-.4,places=8)
    def test_insufficient(self):
        self.assertEqual(fit([])["status"],"insufficient support")
if __name__=="__main__":unittest.main()
