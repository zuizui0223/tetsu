import unittest
import numpy as np
from analyze_world_railstatus_bird_temporal_influence import average,deletion
class TestTemporal(unittest.TestCase):
    def test_shared_year_cluster_covariance(self):
        m={y:{"keys":[("X",0),("X",1)],"effects":{"active":{"estimate":2.,"scores":np.array([.1,-.1])}}} for y in ("2023","2024")}
        z=average(m,("2023","2024"),"active")
        self.assertAlmostEqual(z["estimate"],2.)
        self.assertAlmostEqual(z["interval95_approx"][1]-2.,1.96*.2)
    def test_deletion_sufficient_statistics(self):
        m={"keys":[("X",i) for i in range(25)],"grams":np.ones((25,1,1)),"xy":np.ones((25,1))*2,
           "terms":{"active":np.array([1.])}}
        self.assertEqual(deletion(m,{("X",0)})["active"],2.)
        self.assertIsNone(deletion(m,{("X",i) for i in range(10)}))
if __name__=="__main__":unittest.main()
