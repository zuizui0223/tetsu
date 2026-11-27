import unittest
from analyze_world_railstatus_plant_size_temporal import common,temporal
class TestTemporal(unittest.TestCase):
    def test_empty(self):
        z=temporal({y:[] for y in ("2023","2024","2025")})
        self.assertEqual(z["n_common_cells"],0)
        self.assertIsNone(z["earlier_mean"])
    def test_block_drift(self):
        with self.assertRaises(ValueError):
            common({y:[{"key":(0,0),"continent":"X","block":y}] for y in ("2023","2024","2025")})
if __name__=="__main__":unittest.main()
