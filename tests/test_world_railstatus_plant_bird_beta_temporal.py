import unittest
import numpy as np
from analyze_world_railstatus_plant_bird_beta_temporal import combine,interval,analyze
class TestTemporal(unittest.TestCase):
    def test_cluster_covariance(self):
        f={"2023":(.2,{("X",0):.02,("X",1):-.02}),"2024":(.4,{("X",0):.02,("X",1):-.02})}
        q=interval(combine(f,{y:.5 for y in f}))
        self.assertAlmostEqual(q["estimate"],.3)
        self.assertAlmostEqual(q["interval95_approx"][1]-.3,1.96*.04)
    def test_empty_support(self):
        z=analyze({y:[] for y in ("2023","2024","2025")})
        self.assertEqual(z["n_common_pairs"],0)
        self.assertIsNone(z["models"]["block_Jaccard"]["earlier_2023_2024"])
if __name__=="__main__":unittest.main()
