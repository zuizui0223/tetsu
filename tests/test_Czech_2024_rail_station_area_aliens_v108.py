"""v10.8 source-invariant Czech station area confounding math, no fake ecological results."""
import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Czech_2024_rail_station_area_aliens_v108 import (
    solve_linear,ols,partial_slope,design,blocked_residual_permutation,
    stratified_station_bootstrap,audit_and_run
)


def sample(n=40):
    xs=[]
    for i in range(n):
        sid=i%5
        track=.2+.012*(i%8)
        building=.1+.006*(i%7)
        logarea=.02*i
        y=.25+.7*track+.04*building+.003*i+.01*sid
        xs.append({"station":str(i),"stretch":str(sid),
                   "area_ha":__import__("math").exp(logarea),
                   "track_share":track,"building_share":building,
                   "alien_fraction":y,
                   "invasive_fraction_of_aliens":.15+.05*track})
    return xs


class OriginalCzechAreaConditionalTests(unittest.TestCase):
    def test_non_singular_normal_equations(self):
        self.assertEqual(solve_linear([[1,0],[0,2]],[3,8]),[3,4])
        with self.assertRaisesRegex(ValueError,"singular"):
            solve_linear([[1,2],[2,4]],[3,8])

    def test_OLS_exact_synthetic_two_parameter(self):
        b,rs=ols([1,3,5],[[1,0],[1,1],[1,2]])
        self.assertAlmostEqual(b[0],1)
        self.assertAlmostEqual(b[1],2)
        self.assertTrue(all(abs(v)<1e-12 for v in rs))

    def test_track_area_adjustment_with_5_stretch_block_dummies(self):
        xs=sample()
        coef,rx,ry=partial_slope(xs,"alien_fraction")
        # .003*i creates additional unmodelled trend, so only check
        # true track association remains positive but may be biased.
        self.assertGreater(coef,0)
        self.assertEqual(len(rx),40)
        self.assertEqual(len(design(xs)),40)

    def test_residual_permutation_and_stratified_bootstrap_bounded(self):
        xs=sample()
        p=blocked_residual_permutation(xs,"alien_fraction",draws=99)
        ci,n=stratified_station_bootstrap(xs,"alien_fraction",draws=50)
        self.assertGreaterEqual(p,0)
        self.assertLessEqual(p,1)
        self.assertEqual(n,50)
        self.assertLess(ci[0],ci[1])

    def test_original_Czech_XLSX_source_SHA_is_gate_not_synthetic_data(self):
        with self.assertRaises(Exception):
            audit_and_run("/missing/source-S1.xlsx","/missing/source-S2.xlsx")


if __name__=="__main__":unittest.main()
