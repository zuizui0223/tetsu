import unittest
from analyze_world_railstatus_bird_conditional_hwi import common_stratum_weights,standardized,AXES,HWI,YEARS
class TestConditionedHWI(unittest.TestCase):
    def fixture(self,n,j):
        return {"split_counts_a":[[20,4,n],[20,4,n]],"split_counts_b":[[20,4,n],[20,4,n]],
                "split_similarity":{k:{"jaccard":j,"simpson_similarity":j} for k in ("00","01","10","11")}}
    def test_fixed_strata(self):
        z={y:{"guilds":{}} for y in YEARS}
        for y in YEARS:
            for k,n in (("h1_m1",8),("h2_m1",12)):
                for g,j in zip(AXES[k],(.2,.3)):
                    z[y]["guilds"][g]=self.fixture(n,j)
            for g,j in zip(HWI,(.1,.4)):z[y]["guilds"][g]=self.fixture(20,j)
        z["2023"]["guilds"][AXES["h1_m1"][0]]=self.fixture(6,.2)
        w=common_stratum_weights(z,("h1_m1","h2_m1"))
        self.assertAlmostEqual(w["h1_m1"],1/3)
        q=standardized({},z["2024"],w)
        self.assertAlmostEqual(q["difference_jaccard"],.1)
        self.assertAlmostEqual(q["standardized_minus_raw_jaccard"],-.2)
        del z["2025"]["guilds"][AXES["h2_m1"][1]]
        self.assertEqual(common_stratum_weights(z,("h1_m1","h2_m1")),{"h1_m1":1.})
if __name__=="__main__":unittest.main()
