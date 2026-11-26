import unittest
from scan_world_railstatus_disturbance_traits import binary,summarize,family_refs
from analyze_world_railstatus_disturbance_traits import attach,estimate
class TestTraits(unittest.TestCase):
    def test_missing_denominator(self):
        t={1:{"annual":1.,"woody":0.,"annual_nonwoody":1.},2:{"annual":0.,"woody":0.,"annual_nonwoody":0.},3:{"annual":0.,"woody":1.}}
        z=summarize({1,2,3,4},t,{},{},"plants")
        self.assertEqual(z["annual"]["denominator"],4)
        self.assertEqual(z["annual"]["known"],3)
        self.assertEqual(z["annual_nonwoody"]["denominator"],2)
        self.assertAlmostEqual(z["annual_nonwoody"]["mean"],.5)
        self.assertIsNone(binary("unknown",{"yes"},{"no"}))
    def test_family_leave_self_out(self):
        t={i:{"HWI":float(i)} for i in range(12)}
        f={i:1 for i in t};r=family_refs(t,f,"birds")
        z=summarize({0},t,f,r,"birds")
        self.assertEqual(z["HWI"]["family_known"],1)
        self.assertAlmostEqual(z["HWI"]["family_excess"],-6.)
    def test_empty(self):
        self.assertNotEqual(estimate([])[0]["status"],"estimated")
if __name__=="__main__":unittest.main()
