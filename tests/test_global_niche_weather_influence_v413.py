import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_global_niche_weather_influence_v413 import summary,evaluate,TARGET,COLD,DRY_INT,D
class TestV413(unittest.TestCase):
    def test_fixed_temperature_slope(self):
        f={"coef":{TARGET:{D:2.},COLD:{D:-1.},DRY_INT:{D:.1}}}
        s=summary(f);self.assertEqual(s["dry_rail_15C_zero_decline"],2.)
        self.assertEqual(s["dry_rail_5C_zero_decline"],1.)
        self.assertEqual(s["dry_rail_x_precip_decline"],.1)
    def test_no_model_below_gate(self):
        self.assertEqual(evaluate([],())["status"],"insufficient support")
if __name__=="__main__":unittest.main()
