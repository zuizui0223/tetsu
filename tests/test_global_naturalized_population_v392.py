import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from aggregate_original_GHS_POP2020_v391 import ROWS,COLS,output_rows,roll
from analyze_global_naturalized_population_v392 import add_population
class TestPop(unittest.TestCase):
    def test_exact_rollup_and_orientation(self):
        p=np.zeros((ROWS,COLS));a=p.copy()
        p[ROWS-1,0]=42.;a[ROWS-1,0]=1e9
        p[ROWS-2,1]=100.;a[ROWS-2,1]=2e9
        self.assertEqual(roll(p,2).sum(),142.)
        r=output_rows(p,a,2)
        self.assertEqual(len(r),1);self.assertEqual(r[0]["source_equal_area_grid_iy"],0)
        self.assertEqual(r[0]["source_GHS_POP2020_people"],142.)
    def test_missing_is_not_zero_and_land_density(self):
        rows=[dict(key=(0,0)),dict(key=(0,1))]
        road={(0,i):{"source_GRIP4_total_land_area_m2":1e9} for i in range(2)}
        pop={(0,0):{"source_GHS_POP2020_people":10000.,"source_GHS_POP2020_approx_valid_pixel_area_m2":1e9}}
        r,qa=add_population(rows,pop,road)
        self.assertAlmostEqual(r[0]["log_population_density"],math.log1p(10.))
        self.assertEqual(qa["no_population_source"],1)
if __name__=="__main__":unittest.main()
