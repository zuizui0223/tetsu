import unittest,json
import numpy as np
from analyze_world_railstatus_direct_boundary_beta import fit,ACTIVE,INACTIVE
class TestBoundaryJSON(unittest.TestCase):
    def test_numpy_support_serializes(self):
        rng=np.random.default_rng(11);rows=[]
        for b in range(30):
            for i in range(6):
                a=float(i%2);z=float(rng.integers(0,2))
                rows.append({"continent":"X","block":b,"a":(b,i),"b":(b,i+1),ACTIVE:a,INACTIVE:z,"jaccard":.3+.1*a-.04*z})
        r=fit(rows,"jaccard",False,(ACTIVE,INACTIVE))
        self.assertEqual(r["status"],"estimated")
        self.assertAlmostEqual(r["effects"][ACTIVE]["estimate"],.1)
        self.assertEqual(r["n_blocks_residual_active"],30)
        json.dumps(r,allow_nan=False)
if __name__=="__main__":unittest.main()
