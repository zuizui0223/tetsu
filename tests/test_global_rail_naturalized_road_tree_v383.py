import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_rail_naturalized_road_tree_v383 import conditional_coefs,eligible_region_rows
class TestV383(unittest.TestCase):
    def test_region_demean_and_nuisance(self):
        groups=[np.arange(0,5),np.arange(5,10)]
        x=np.tile(np.arange(5.),2);z=np.tile([0.,1.,0.,1.,2.],2)
        y=2*x+3*z+np.repeat([10.,100.],5)
        X=np.column_stack([x,z,y])
        r=conditional_coefs(X,groups,np.ones(10),("x","z"),("x","z"))
        self.assertAlmostEqual(r["x"],2.);self.assertAlmostEqual(r["z"],3.)
        X[:,1]=X[:,0]
        self.assertIsNone(conditional_coefs(X,groups,np.ones(10),("x","z"),("x",))["x"])
    def test_zero_bootstrap_group(self):
        X=np.array([[0.,0.],[1.,2.],[0.,100.],[1.,102.]])
        r=conditional_coefs(X,[np.array([0,1]),np.array([2,3])],np.array([1.,1.,0.,0.]),("x",),("x",))
        self.assertAlmostEqual(r["x"],2.)
    def test_support_is_regional_not_global_count(self):
        rows=[dict(region="A",continent="Asia",block=(i//3,0),log_rail_km=float(i)) for i in range(6)]
        rows+= [dict(region="B",continent="Asia",block=(0,0),log_rail_km=float(i)) for i in range(6)]
        self.assertEqual(len(eligible_region_rows(rows)),6)
if __name__=="__main__":unittest.main()
