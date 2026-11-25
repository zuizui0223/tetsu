import unittest,math
import numpy as np
from analyze_world_railstatus_plant_insect_bird import extra_attach,fields_for,linear,combine,interval
class TestTriad(unittest.TestCase):
    def test_same_cell_attachment(self):
        r={"key":(1,2)};m={"insects":{(1,2):((60,10,20),((30,5,10),(30,5,10)))}}
        f={"insects":{(1,2):(60,10,20)}}
        z=extra_attach([r],m,f,("insects",))
        self.assertAlmostEqual(z[0]["insects_0_species"],math.log1p(10))
        f["insects"][1,2]=(61,10,20)
        with self.assertRaises(ValueError):extra_attach([r],m,f,("insects",))
    def test_paired_attenuation(self):
        rng=np.random.default_rng(9);rows=[]
        for block in range(30):
            for i in range(10):
                p,e=rng.normal(size=2);a=.8*p+e
                rows.append({"continent":"X","block":block,"p":p,"a":a,"y":.4*p+.6*a+block})
        u=linear(rows,("p",),"y","p");v=linear(rows,("p","a"),"y","p")
        self.assertAlmostEqual(v[0],.4,places=8)
        self.assertGreater(interval(combine({"u":u,"v":v},{"u":1.,"v":-1.}))["estimate"],.2)
    def test_fields(self):
        f=fields_for(("insects","fungi"),True)
        self.assertEqual(len(f),len(set(f)))
if __name__=="__main__":unittest.main()
