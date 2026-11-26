import unittest
from analyze_world_railstatus_annual_triangle import annual_panel,fit_triangle,KEYS
class TestTriangle(unittest.TestCase):
    def test_empty(self):
        z=annual_panel({y:[] for y in ("2023","2024","2025")})
        self.assertEqual(z["n_common_cells"],0)
        self.assertIsNone(z["models"]["flexible_climate"]["earlier_mean"])
    def test_unique_edges(self):
        self.assertEqual(len(set(KEYS)),3)
        self.assertIn("plants__insects",KEYS)
    def test_fit_empty(self):
        self.assertIsNone(fit_triangle([])[1])
if __name__=="__main__":unittest.main()
