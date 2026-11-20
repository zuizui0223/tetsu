import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_Slovenia_traffic_urbanization_v305 import fit, model
from extract_Slovenia_worldcover_buffers_v305 import summarize_arrays

class TestUrbanTraffic(unittest.TestCase):
    def test_adjustment_recovers_known_confounding(self):
        rows=[{'x':i%2,'u':i/10+(i%2)*2,'y':3*(i%2)+5*(i/10+(i%2)*2)} for i in range(20)]
        result=fit(rows,'y',['x','u'])
        self.assertAlmostEqual(result['coefficients']['x'],3)
        self.assertAlmostEqual(result['coefficients']['u'],5)
    def test_equal_group_weights_resist_duplication(self):
        rows=[{'group':str(i),'x':i%2,'y':float(i*i)} for i in range(8)]
        duplicate=rows+[dict(rows[0]) for _ in range(30)]
        self.assertAlmostEqual(fit(rows,'y',['x'],'group')['coefficients']['x'],
                               fit(duplicate,'y',['x'],'group')['coefficients']['x'])
    def test_saturated_model_stops(self):
        rows=[{'x':i,'y':i} for i in range(4)]
        self.assertEqual(fit(rows,'y',['x'])['status'],'not_identifiable')
    def test_leave_group_removes_entire_cluster(self):
        rows=[{'section':str(i//3),'x':i,'y':2*i+1} for i in range(18)]
        result=model(rows,'y',['x'],'x')
        self.assertEqual(result['leave_one_group_out']['n_estimable'],6)
        self.assertAlmostEqual(result['leave_one_group_out']['min'],2)
    def test_nodata_not_absence_and_area_weighting(self):
        result=summarize_arrays(np.array([50,10,0]),np.array([200,200,200]),np.array([1,3,4]),500)
        self.assertEqual(result['valid_area_fraction'],.5)
        self.assertEqual(result['built_up_fraction'],.25)
    def test_inner_railway_footprint_exclusion(self):
        result=summarize_arrays(np.array([50,10]),np.array([50,300]),np.array([1,1]),500,100)
        self.assertEqual(result['built_up_fraction'],0)
        self.assertEqual(result['n_pixels'],1)

if __name__=='__main__':unittest.main()
