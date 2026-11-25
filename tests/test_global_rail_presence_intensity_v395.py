import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_rail_presence_intensity_v395 import add_exposure,mixed,coefficients
class TestV395(unittest.TestCase):
    def test_reference_and_zero(self):
        for n in (0,100):
            r=add_exposure(dict(rail_km=n,log_rail_km=math.log1p(n)))
            self.assertEqual(r["rail_present"],float(n>0))
            self.assertAlmostEqual(r["rail_positive_centered"],0.)
    def test_mixed_support(self):
        rows=[dict(region=k,continent="Asia",block=(0,0),rail_present=p) for k,p in [("a",0),("a",1),("b",1)]]
        self.assertEqual(len(mixed(rows)),2)
    def test_count_contrast(self):
        x=np.tile([0.,1.,1.,0.],2);z=np.tile([0.,1.,0.,1.],2);g=np.repeat([10.,20.],4)
        a=2*x+z+g;b=x+z+2*g
        X=np.column_stack([x,z,a,b,a-b,x])
        c=coefficients(X,[np.arange(4),np.arange(4,8)],np.ones(8),("x","z"),("x",))
        self.assertAlmostEqual(c["x"]["log1p_documented"],2.)
        self.assertAlmostEqual(c["x"]["log1p_count_contrast"],1.)
if __name__=="__main__":unittest.main()
