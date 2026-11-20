import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_region_rail_elevation_transfer_v302 import fit_foreign, predict_score, region_holdout, rail_score

class TransferTests(unittest.TestCase):
    def setUp(self):
        self.x=np.arange(6,dtype=float)/5
        self.arrays={r:[(self.x, 2-self.x)] for r in ('A','B','C')}

    def test_perfect_transfer_and_baseline(self):
        model=fit_foreign(self.arrays)
        self.assertAlmostEqual(model['slope'],-1)
        for anchor in ('lowest','median','none'):
            score=predict_score(self.x,2-self.x,model,anchor)
            self.assertAlmostEqual(score['rmse'],0)
            self.assertGreater(score['gain'],0)

    def test_excluded_region_cannot_change_training(self):
        first=region_holdout(self.arrays)['C']
        self.arrays['C']=[(self.x,np.array([2,8,1,0,3,9.]))]
        second=region_holdout(self.arrays)['C']
        self.assertEqual(first['model'],second['model'])
        self.assertNotEqual(first['scores']['lowest']['gain'],second['scores']['lowest']['gain'])

    def test_region_equal_not_site_equal(self):
        self.arrays['A']=[(self.x,2-self.x)]*10
        self.arrays['B']=[(self.x,2+2*self.x)]
        self.arrays['C']=[(self.x,2+2*self.x)]
        self.assertAlmostEqual(fit_foreign(self.arrays)['slope'],1)

    def test_single_anchor_level_offset(self):
        model=fit_foreign(self.arrays)
        a=predict_score(self.x,2-self.x,model,'lowest')
        b=predict_score(self.x,9-self.x,model,'lowest')
        self.assertAlmostEqual(a['gain'],b['gain'])
        self.assertEqual(a['evaluated_sites'],5)
        cold=predict_score(self.x,9-self.x,model,'none')
        self.assertGreater(cold['rmse'],6)

    def test_transfer_failure_is_signed(self):
        model=fit_foreign(self.arrays)
        score=predict_score(self.x,2+self.x,model,'lowest')
        self.assertLess(score['gain'],0)

    def test_training_count(self):
        with self.assertRaises(ValueError):
            fit_foreign({'A':self.arrays['A']})

    def test_rail_support_not_zero_effect(self):
        rows=[dict(x=x/100,y=1,site=str(i)) for i,x in enumerate(self.x)]
        out=rail_score(rows,self.arrays)
        self.assertEqual(out['status'],'insufficient_sites_or_elevation')
        self.assertNotIn('scores',out)

    def test_target_elevation_extrapolation_reported(self):
        model=fit_foreign(self.arrays)
        score=predict_score(self.x+2,2-self.x,model,'none')
        self.assertEqual(score['outside_training_elevation'],6)

if __name__=='__main__':
    unittest.main()
