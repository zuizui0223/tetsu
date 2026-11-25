import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_niche_pool_common_v407 import verify_additive,TAGS
class TestV407(unittest.TestCase):
    def test_additive_null_contrast(self):
        c={}
        for t in TAGS:
            c.update({"observed_"+t:-2,"occupancy_expected_"+t:-1,"occupancy_excess_"+t:-1,
                      "family_occupancy_expected_"+t:-3,"family_occupancy_excess_"+t:1,
                      "family_occupancy_minus_occupancy_"+t:-2})
        verify_additive({"coef":{"rail_cold":c}},("occupancy","family_occupancy"))
        c["family_occupancy_excess_T_variance"]=0
        with self.assertRaises(ValueError):verify_additive({"coef":{"rail_cold":c}},("occupancy","family_occupancy"))
if __name__=="__main__":unittest.main()
