import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from aggregate_original_GRIP4_road_density_v381a import weights,area_aggregate
class TestGRIP(unittest.TestCase):
    def test_exact_overlap_and_missing(self):
        x=weights([0,2,4],[0,1,3,4]);y=weights([2,0],[0,2])
        n,a=area_aggregate(np.array([[10.,30.]]),np.array([[True,False]]),x,y)
        np.testing.assert_allclose(a,[[2.,2.,0.]])
        np.testing.assert_allclose(n,[[20.,20.,0.]])
    def test_valid_zero_retained(self):
        x=weights([0,2],[0,2]);y=weights([0,2],[0,2])
        n,a=area_aggregate(np.array([[0.]]),np.array([[True]]),x,y)
        np.testing.assert_allclose(n,[[0.]])
        np.testing.assert_allclose(a,[[4.]])
if __name__=="__main__":unittest.main()
