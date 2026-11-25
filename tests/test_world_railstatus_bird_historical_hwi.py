import unittest
from scan_world_railstatus_bird_historical_hwi import period,PERIOD
from analyze_world_railstatus_bird_historical_hwi import weights,panel,AXES
class TestHistorical(unittest.TestCase):
    def test_period(self):
        for y in range(2018,2023):self.assertEqual(period(f"{y}-01-01"),PERIOD)
        for y in (2017,2023,2024,2026):self.assertIsNone(period(f"{y}-01-01"))
        self.assertEqual(period("2025-01-01"),"2025")
    def test_fixed_weights(self):
        p={"guilds":{}}
        for a,n in (("h1_m1",6),("h2_m1",12)):
            for g in AXES[a]:p["guilds"][g]={"split_counts_a":[[20,4,n]]*2,"split_counts_b":[[20,4,n]]*2}
        self.assertAlmostEqual(weights(p,tuple(AXES))["h1_m1"],1/3)
    def test_empty(self):
        r=panel([],{});self.assertEqual(r["n_raw"],0)
        self.assertNotEqual(r["conditioned"]["habitat_x_migration"]["models"]["block_FE"]["difference_jaccard"]["status"],"estimated")
if __name__=="__main__":unittest.main()
