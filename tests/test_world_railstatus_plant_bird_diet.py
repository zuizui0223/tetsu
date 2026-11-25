import unittest,math
from scan_world_railstatus_bird_diet_cells import category
from analyze_world_railstatus_plant_bird_diet import diet_attach,one_diet
class TestDiet(unittest.TestCase):
    def test_category(self):
        self.assertEqual(category("Herbivore, terrestrial"),"herbivoreterrestrial")
        self.assertEqual(category("Aquatic predator"),"aquaticpredator")
        self.assertEqual(category("NA"),"na")
        self.assertEqual(category(None),"")
    def test_counts(self):
        h=[[30,5,10],[30,5,10]]
        r={"key":(1,2)}
        for i in (0,1):
            for k,v in zip(("observations","observers","species"),h[i]):r[f"birds_{i}_{k}"]=math.log1p(v)
        q={(1,2):{"halves":h,"guilds":{"fruit_seed":[[15,4,5]]*2,"invertivore":[[15,4,5]]*2}}}
        z=diet_attach([r],q,("fruit_seed","invertivore"))
        self.assertEqual(z[0]["d0_0_species"],math.log1p(5))
        del q[(1,2)]["guilds"]["invertivore"]
        self.assertEqual(diet_attach([r],q,("fruit_seed","invertivore")),[])
    def test_empty(self):
        self.assertNotEqual(one_diet([])[0]["status"],"estimated")
if __name__=="__main__":unittest.main()
