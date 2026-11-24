import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_naturalized_count_decomposition_v389 import count_outcomes,conditional
class TestV389(unittest.TestCase):
    def test_original_count_partition(self):
        a=count_outcomes(dict(vascular_richness=30,documented_richness=6,species=35,known_fraction=.2))
        self.assertAlmostEqual(a[0],math.log1p(6));self.assertAlmostEqual(a[1],math.log1p(24))
        self.assertAlmostEqual(a[2],a[0]-a[1]);self.assertEqual(a[3],20)
        with self.assertRaises(ValueError):
            count_outcomes(dict(vascular_richness=30,documented_richness=6,species=35,known_fraction=.3))
    def test_region_conditioning_and_additivity(self):
        rail=np.tile([0.,1.,2.,3.,4.],2)
        road=np.tile([0.,1.,0.,1.,2.],2)
        group=np.repeat([10.,50.],5)
        x=2*rail+road+group;y=rail+3*road+2*group
        X=np.column_stack((rail,road,x,y,x-y,rail+road))
        out=conditional(X,[np.arange(5),np.arange(5,10)],np.ones(10),("log_rail_km","log_road_density"))
        self.assertAlmostEqual(out["log_rail_km"]["log1p_documented"],2.)
        self.assertAlmostEqual(out["log_rail_km"]["log1p_unlisted"],1.)
        self.assertAlmostEqual(out["log_rail_km"]["log1p_count_contrast"],1.)
if __name__=="__main__":unittest.main()
