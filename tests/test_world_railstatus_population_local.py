import unittest
import numpy as np
from analyze_world_railstatus_population_local import fit,join_population,local_support,ACTIVE,INACTIVE,POP
class TestPopulationLocal(unittest.TestCase):
    def test_population_density(self):
        r=[{"key":(1,2)}];p={(1,2):{"source_GHS_POP2020_people":99,"source_GHS_POP2020_approx_valid_pixel_area_m2":1e6}}
        rd={(1,2):{"source_GRIP4_total_land_area_m2":1e6}}
        out,_=join_population(r,p,rd)
        self.assertAlmostEqual(out[0][POP],np.log(100))
    def test_block_fixed_effects(self):
        rng=np.random.default_rng(8);rows=[]
        for k in range(30):
            for j in range(5):
                a,b,p=rng.uniform(0,5,3)
                rows.append({"continent":"X","block":k,ACTIVE:a,INACTIVE:b,POP:p,"log_species":.3*a-.1*b+.2*p+10*k})
        z=fit(rows,(ACTIVE,INACTIVE,POP),True)
        self.assertEqual(z["status"],"estimated")
        self.assertAlmostEqual(z["active_minus_inactive"]["estimate"],.4,places=8)
        self.assertEqual(len(local_support(rows+[{"continent":"X","block":31}])),150)
    def test_insufficient_support(self):
        self.assertEqual(fit([], (ACTIVE,INACTIVE), True)["status"],"insufficient support")
if __name__=="__main__":unittest.main()
