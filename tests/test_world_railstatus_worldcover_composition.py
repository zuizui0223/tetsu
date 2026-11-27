import unittest
from analyze_world_railstatus_worldcover_composition import components,decompose,stats,CLASS
class TestComposition(unittest.TestCase):
    def test_identity(self):
        total,amount,allocation,dp,dv=components((.4,.75),(.8,.5))
        self.assertAlmostEqual(total,-.1)
        self.assertAlmostEqual(amount,-.25)
        self.assertAlmostEqual(allocation,.15)
        self.assertAlmostEqual(total,amount+allocation)
    def test_amount_only(self):
        z=components((.4,.5),(.8,.5))
        self.assertAlmostEqual(z[0],-.2)
        self.assertAlmostEqual(z[2],0.)
    def test_empty(self):
        self.assertEqual(stats([],("x",))["status"],"insufficient support")
    def test_denominator_gate(self):
        q={k:0. for k in CLASS.values()};q["built"]=1.
        z={s:{k:{"fractions":q} for k in ("near","background")} for s in ("active","inactive")}
        self.assertIsNone(decompose(z))
if __name__=="__main__":unittest.main()
