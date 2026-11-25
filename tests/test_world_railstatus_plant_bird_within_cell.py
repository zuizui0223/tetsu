import unittest
import numpy as np
from analyze_world_railstatus_plant_bird_within_cell import fit_panel,model
class TestWithin(unittest.TestCase):
    def test_panel(self):
        rng=np.random.default_rng(11);by={y:[] for y in ("2023","2024","2025")}
        for b in range(25):
            for k in range(6):
                fixed=rng.normal()
                for j,y in enumerate(by):
                    x,z=rng.normal(size=2)
                    by[y].append({"key":(b,k),"continent":"X","block":b,"x":x,"z":z,"static":fixed,"y":fixed+j+b*j+.7*x-.4*z})
        q=fit_panel(by,("x","z","static"),"y","x","block_year")
        self.assertAlmostEqual(q[0],.7,places=8)
        self.assertIsNone(fit_panel(by,("x","static"),"y","static","year"))
    def test_unbalanced(self):
        with self.assertRaises(ValueError):fit_panel({"a":[{"key":(1,2)}],"b":[]} ,("x",),"y","x")
    def test_empty(self):
        self.assertNotEqual(model({"a":[],"b":[]})["status"],"estimated")
if __name__=="__main__":unittest.main()
