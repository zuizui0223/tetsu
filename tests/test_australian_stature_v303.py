import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_australian_stature_v303 import estimate, habitat


class StatureTests(unittest.TestCase):
    def test_known_linear_gradient(self):
        rows=[dict(elevation_m=400+i,log_height=3-.002*i,habitat='plain',block=i//5) for i in range(30)]
        result=estimate(rows)
        self.assertAlmostEqual(result['slope_log_height_per_100m'],-.2,places=10)
        for x in result['leave_block_out_range']:
            self.assertAlmostEqual(x,-.2,places=10)

    def test_habitat_confounding_removed(self):
        rows=[]
        for g in range(2):
            for i in range(20):
                rows.append(dict(elevation_m=400+g*40+i,log_height=float(g),habitat=str(g),block=i//5))
        self.assertGreater(estimate(rows)['slope_log_height_per_100m'],0)
        self.assertAlmostEqual(estimate(rows,True)['slope_log_height_per_100m'],0,places=10)

    def test_habitat_is_metadata_based(self):
        self.assertEqual(habitat('Floodplain - sandplain.'),'drainage_floodplain')
        self.assertEqual(habitat('Sandplain.'),'sandplain')
        self.assertEqual(habitat('Hill top.'),'other')


if __name__=='__main__':unittest.main()
