import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_niche_shared_taxa_v401 import summary,rarefy,metrics,ref_subset,MODES,DRAWS
class TestV401(unittest.TestCase):
    def test_exact_twenty_cell_rarefaction(self):
        x=np.column_stack((np.arange(20.),np.ones(20)))
        v=rarefy(x,123)
        self.assertEqual(v.shape,(DRAWS,5))
        np.testing.assert_allclose(v,np.repeat(summary(x)[None,:],DRAWS,axis=0))
    def test_fixed_size_seed(self):
        x=np.column_stack((np.arange(50.),np.ones(50)))
        np.testing.assert_array_equal(rarefy(x,123),rarefy(x,123))
        self.assertTrue((rarefy(x,123)[:,2]>0).all())
    def test_reference_filters(self):
        a=np.arange(4);b=np.array([[0,0],[1,0],[2,0],[3,0]])
        c=np.array(["Asia","Europe","Asia","Europe"]);r={"block":(0,0),"continent":"Asia"}
        np.testing.assert_array_equal(ref_subset(a,b,c,r,MODES[0]),[2,3])
        np.testing.assert_array_equal(ref_subset(a,b,c,r,MODES[1]),[1,3])
    def test_community_metrics_use_draw_specific_sd(self):
        a=np.array([[[1,0,2,1,1],[3,1,4,2,2]],[[3,0,2,1,1],[1,1,4,2,2]]],float)
        r=metrics(a,{"bio1":0,"log1p_bio12":1})
        self.assertAlmostEqual(r["warm_offset_C"],2)
        self.assertAlmostEqual(r["among_species_temperature_sd_C"],1)
        self.assertAlmostEqual(r["temp_width_C"],3)
if __name__=="__main__":unittest.main()
