import unittest
from analyze_world_railstatus_trait_dispersion import fit,attach_variance
class TestDispersion(unittest.TestCase):
    def test_empty(self):
        self.assertNotEqual(fit([])[0]["status"],"estimated")
    def test_moments(self):
        r={"key":(0,0),"active_normal_km":1.,"inactive_normal_km":1.}
        q={"counts":[60,10,20],"traits":{"HWI":{"denominator":20,"known":20,"mean":30.,"variance":100.,"family_known":0}}}
        src={(0,0):{"target_full":q}}
        z=attach_variance([r],src,"HWI",.5)[0]
        self.assertAlmostEqual(z["trait_value"],20/19)
        self.assertEqual(z["trait_mean"],3.)
if __name__=="__main__":unittest.main()
