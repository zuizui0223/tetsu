import unittest
from analyze_world_railstatus_triangle_within_cell import panel_fields,estimate
class TestTriangleWithin(unittest.TestCase):
    def test_effort_fields(self):
        self.assertEqual(len(panel_fields(False)),12)
        self.assertEqual(len(panel_fields(True)),42)
        self.assertEqual(len(set(panel_fields(True))),42)
    def test_empty(self):
        r=estimate({"2023":[],"2024":[],"2025":[]})
        self.assertEqual(r["status"],"partial or insufficient support")
        self.assertTrue(all(v is None for v in r["pairwise"].values()))
        self.assertIsNone(r["joint"])
if __name__=="__main__":unittest.main()
