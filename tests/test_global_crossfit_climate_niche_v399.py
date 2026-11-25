import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_crossfit_climate_niche_v399 import envelope,slopes
class TestV399(unittest.TestCase):
    def test_exclude_neighbor_blocks_and_continent(self):
        vals=np.column_stack((np.arange(30.),np.ones(30)))
        blocks=np.column_stack((np.arange(30)//4,np.zeros(30,int)))
        c=np.array(["Asia"]*4+["Europe"]*26)
        a=envelope(vals,blocks,c,(0,0),"Asia","exclude_3x3_blocks")
        self.assertAlmostEqual(a[0],18.5)
        self.assertAlmostEqual(a[1],9.05)
        self.assertAlmostEqual(envelope(vals,blocks,c,(0,0),"Asia","exclude_continent")[0],16.5)
        self.assertIsNone(envelope(vals[:15],blocks[:15],c[:15],(0,0),"Asia","exclude_continent"))
    def test_weighted_local_slope(self):
        x=np.tile([0.,1.,2.,3.],2);y=2*x+np.repeat([10.,30.],4)
        z=np.column_stack([x,y])
        o=slopes(z,[np.arange(4),np.arange(4,8)],np.ones(8),("rail",),("rail",))
        self.assertAlmostEqual(o["rail"][0],2.)
if __name__=="__main__":unittest.main()
