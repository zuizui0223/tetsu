import unittest,math
from analyze_world_railstatus_growthform_diversity import diversity,expected_present,attach
class TestSupportDiversity(unittest.TestCase):
    def test_decomposition(self):
        d=diversity({"F":10,"C":5,"E":5})
        self.assertAlmostEqual(d["support_diversity"],250/380)
        self.assertAlmostEqual(d["free_vs_nonfree"],200/380)
        self.assertAlmostEqual(d["nonfree_within"],50/380)
        self.assertEqual(d["rarefied_support_richness"],3)
    def test_exact_rarefaction(self):
        self.assertAlmostEqual(expected_present(30,3,20),1-math.comb(27,20)/math.comb(30,20))
        self.assertEqual(expected_present(30,0,20),0.)
        self.assertEqual(expected_present(30,30,20),1.)
        self.assertEqual(diversity({"F":20})["support_diversity"],0.)
    def test_missing_not_other(self):
        r={"key":(0,0),"species":100,"observations":200,"observers":25}
        t=dict(r,known_support=30,support_counts={"F":25,"C":5})
        self.assertEqual(attach([r],{(0,0):t},.5),[])
        self.assertAlmostEqual(attach([r],{(0,0):t},.25)[0]["free_standing_fraction"],25/30)
if __name__=="__main__":unittest.main()
