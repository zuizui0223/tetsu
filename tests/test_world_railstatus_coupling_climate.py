import unittest
from analyze_world_railstatus_coupling_climate import EXTRA,attach_climate,fit_model
class TestClimate(unittest.TestCase):
    def test_fields(self):
        self.assertEqual(len(EXTRA),6)
        self.assertEqual(len(set(EXTRA)),6)
    def test_empty(self):
        self.assertEqual(attach_climate([],{}),[])
        self.assertNotEqual(fit_model([],("insects",))[0]["status"],"estimated")
if __name__=="__main__":unittest.main()
