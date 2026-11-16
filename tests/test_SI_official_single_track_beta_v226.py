import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from audit_Slovenia_official_GR_single_vs_double_taxonomic_beta_v226 import equals_stratum
class TestTrackClassMatching(unittest.TestCase):
    def test_G_R_source_line_equal_single_only(self):
        dat=[{'section':'10','annual':1.0},{'section':'10','annual':3.0},
             {'section':'20','annual':2.0},{'section':'21','annual':0.0},
             {'section':'31','annual':0.5}]
        x=equals_stratum({'10','20','21','31'},dat,'annual')
        self.assertEqual(x['n_locations'],5)
        self.assertEqual(x['n_G_main_lines'],2)
        self.assertEqual(x['n_R_regional_lines'],2)
        self.assertAlmostEqual(x['G_vs_R']['main_minus_regional'],1.75)
    def test_no_actual_R_double_does_not_create_fake_zero_effect(self):
        onlyG=[{'section':'10','height_var':0.5}]
        x=equals_stratum({'10'},onlyG,'height_var')
        self.assertEqual(x['n_R_regional_lines'],0)
        self.assertIsNone(x['G_vs_R']['main_minus_regional'])
    def test_unknown_result_missing_not_fake_negative(self):
        rows=[{'section':'10','height_var':None},{'section':'21','height_var':0.}]
        x=equals_stratum({'10','21'},rows,'height_var')
        self.assertIsNone(x['G_vs_R']['main_minus_regional'])
if __name__=='__main__':unittest.main()
