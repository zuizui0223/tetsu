import unittest
import numpy as np
from analyze_world_railstatus_local_worldcover import tile_name,tile_keys,fractions,evaluate
class TestWorldCover(unittest.TestCase):
    def test_tile(self):
        self.assertEqual(tile_name(-48,36),"S48E036")
        self.assertEqual(tile_name(0,-3),"N00W003")
        self.assertEqual(len(list(tile_keys((-.1,-.1,.1,.1)))),4)
    def test_fractions(self):
        z=fractions(np.array([10,30,10,50]),np.ones(4),np.ones(4,dtype=bool))
        self.assertEqual(z["fractions"]["tree"],.5)
        self.assertEqual(z["fractions"]["grass"],.25)
        self.assertIsNone(fractions(np.array([10,0]),np.ones(2),np.ones(2,dtype=bool)))
    def test_gate(self):
        self.assertEqual(evaluate([])["status"],"insufficient support")
if __name__=="__main__":unittest.main()
