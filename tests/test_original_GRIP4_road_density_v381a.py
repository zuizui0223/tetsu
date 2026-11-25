import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from aggregate_original_GRIP4_road_density_v381a import weights,area_aggregate,aggregate_land_density
class TestGRIP(unittest.TestCase):
    def test_land_weighted_density_and_unknown(self):
        x=np.array([0.,2000.,4000.]);y=np.array([2000.,0.])
        xw=weights(x,[0,4000]);yw=weights(y,[0,2000])
        # Each source pixel is 4 km2, but land is 1 and 3 km2.
        land,a,n=aggregate_land_density(np.array([[100.,300.]]),np.array([[True,True]]),
                 np.array([[1.,3.]]),np.array([[True,True]]),x,y,xw,yw)
        np.testing.assert_allclose(land,[[4e6]])
        np.testing.assert_allclose(n/a,[[250.]])
        land,a,n=aggregate_land_density(np.array([[100.,300.]]),np.array([[True,False]]),
                 np.array([[1.,3.]]),np.array([[True,True]]),x,y,xw,yw)
        np.testing.assert_allclose(a,[[1e6]])
        np.testing.assert_allclose(n/a,[[100.]])
        np.testing.assert_allclose(land,[[4e6]])
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
