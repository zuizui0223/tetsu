import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from validate_global_rail_local_matched_placebo_v385 import FEATURES,match_pairs,swapped
class AlwaysSwap:
    def integers(self,*a):return 1
class TestV385(unittest.TestCase):
    def test_pairing_no_outcome_or_exposure_and_block(self):
        rows=[]
        for i in range(4):
            r={k:1. for k in FEATURES}
            r.update(continent="Asia",block=(i//2,0),rail_km=float(i),log_rail_km=math.log1p(i),log_species=float(i))
            rows.append(r)
        m,qa=match_pairs(rows)
        self.assertEqual(m,[(0,1),(2,3)])
        for r in rows:r["rail_km"]=100-r["rail_km"];r["log_species"]=999.
        self.assertEqual(match_pairs(rows)[0],m)
        fake,c=swapped(rows,m,AlwaysSwap())
        self.assertEqual(c,4)
        self.assertEqual([r["rail_km"] for r in fake],[99.,100.,97.,98.])
        self.assertEqual([r["rail_km"] for r in rows],[100.,99.,98.,97.])
        self.assertAlmostEqual(fake[0]["log_rail_km"],math.log1p(99.))
    def test_large_covariate_distance_unmatched(self):
        rows=[]
        for i in range(2):
            r={k:float(i) for k in FEATURES};r.update(continent="Asia",block=(0,0),rail_km=0.)
            rows.append(r)
        self.assertEqual(match_pairs(rows)[0],[])
if __name__=="__main__":unittest.main()
