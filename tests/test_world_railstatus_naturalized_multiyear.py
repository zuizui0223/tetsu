import unittest
import numpy as np
from analyze_world_railstatus_naturalized_multiyear import fit,supported_local,ACTIVE,INACTIVE,POP
class TestNaturalizedStatus(unittest.TestCase):
    def test_count_double_contrast(self):
        rng=np.random.default_rng(123);rows=[]
        for k in range(30):
            for j in range(5):
                a,b,p=rng.uniform(0,5,3)
                d=.4*a-.1*b+.2*p+k;u=.1*a+.05*b+.3*p+2*k
                rows.append({"region":"R","continent":"X","block":k,ACTIVE:a,INACTIVE:b,POP:p,"log1p_documented":d,"log1p_unlisted":u,"log1p_count_contrast":d-u})
        z=fit(rows,(ACTIVE,INACTIVE,POP),True)
        self.assertEqual(z["status"],"estimated")
        self.assertAlmostEqual(z["effects"]["active_minus_inactive"]["log1p_count_contrast"]["estimate"],.45,places=8)
    def test_support(self):
        self.assertEqual(fit([])["status"],"insufficient support")
if __name__=="__main__":unittest.main()
