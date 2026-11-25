import unittest
from analyze_world_railstatus_bird_diet_insect import estimate,combine,interval
class TestDietInsect(unittest.TestCase):
    def test_empty(self):
        self.assertNotEqual(estimate([])["status"],"estimated")
    def test_direct_difference(self):
        z={0:(.2,{("X",0):.03,("X",1):-.03}),1:(.5,{("X",0):.03,("X",1):-.03})}
        q=interval(combine(z,{0:-1.,1:1.}))
        self.assertAlmostEqual(q["estimate"],.3)
        self.assertEqual(q["interval95_approx"][0],q["interval95_approx"][1])
if __name__=="__main__":unittest.main()
