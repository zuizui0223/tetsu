import unittest
import numpy as np
from analyze_world_railstatus_plant_bird_coupling_modifiers import combine,linear,interval,design,MODS
class TestModifiers(unittest.TestCase):
    def test_interaction(self):
        rng=np.random.default_rng(21);rows=[]
        for block in range(30):
            for i in range(12):
                x,a,z=rng.normal(size=3)
                rows.append({"continent":"X","block":block,"x":x,"a":a,"z":z,"xa":x*a,"xz":x*z,"y":block+.4*x+.6*a-.2*z+1.2*x*a-.3*x*z})
        fields=("x","a","z","xa","xz")
        f={k:linear(rows,fields,"y",k) for k in ("xa","xz")}
        self.assertAlmostEqual(interval(combine(f,{"xa":1,"xz":-1}))["estimate"],1.5,places=8)
    def test_distinct_modifiers(self):
        self.assertEqual(len(MODS),len(set(MODS)))
if __name__=="__main__":unittest.main()
