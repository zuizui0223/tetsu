import unittest
import scan_world_railstatus_butterfly_edge_annual as s
import analyze_world_railstatus_butterfly_edge_annual as a
class TestAnnual(unittest.TestCase):
    def test_gate(self):
        self.assertTrue(s.full_gate((50,10,20),((25,5,10),(25,5,10))))
        self.assertFalse(s.full_gate((50,10,20),((25,4,10),(25,6,10))))
    def test_validate(self):
        s.validate_close({"x":[.1,2]},{"x":(.1,2)})
        with self.assertRaises(ValueError):s.validate_close({"x":1},{"x":2})
    def test_empty(self):
        z,_=a.compare_years({y:[] for y in s.YEARS})
        self.assertEqual(z["n_common"],0)
        self.assertIsNone(z["earlier_mean"])
if __name__=="__main__":unittest.main()
