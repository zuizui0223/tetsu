"""v316 31-original-transect temperature-affinity geography and exchangeability tests.

Synthetic checks only. Authentic source effects require publisher original SHA gates.
"""
import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
    import numpy as np
    from analyze_Kashmir_paired_thermal_altitude_v316 import (
        slope_model, source_site_permutation, min_support,
        make_pair, build_normalized_plot_pairs,
    )
except ModuleNotFoundError as e:
    if e.name not in ("numpy","analyze_Kashmir_paired_thermal_altitude_v316"):raise
    np=None

@unittest.skipUnless(np is not None,"v316 source model dependencies not installed")
class TestTruePairedThermalAltitude(unittest.TestCase):
    def synthetic(self):
        rows=[]
        for i in range(31):
            stratum="station" if i<16 else "between"
            altitude=1500+13*i
            # Difference = temperature-affinity slope per 100m, not actual Celsius.
            effect=0.2+0.4*(altitude-1500)/100+(0.5 if stratum=="station" else 0.)
            rows.append({"altitude":altitude,"stratum":stratum,"difference":effect,
                         "n_rail_numeric":5,"n_interior_numeric":7})
        return rows
    def test_2014_paired_regression_reproduces_known_slope(self):
        a=slope_model(self.synthetic())
        self.assertTrue(a["estimable"])
        self.assertAlmostEqual(a["altitude_slope_per_100m"],0.4,places=10)
        self.assertAlmostEqual(a["station_effect_at_mean_altitude"],0.5,places=10)
        self.assertEqual(a["n_transects"],31)
        self.assertEqual(a["n_station"],16)
        self.assertEqual(a["n_between"],15)
    def test_low_source_numeric_support_does_not_enter_model(self):
        rows=self.synthetic()
        rows[0]["n_rail_numeric"]=2
        rows[1]["n_interior_numeric"]=0
        kept=min_support(rows)
        self.assertEqual(len(kept),29)
        for row in kept:self.assertGreaterEqual(min(row["n_rail_numeric"],row["n_interior_numeric"]),3)
    def test_drop_station_stratum_blocks_inference(self):
        rows=[r for r in self.synthetic() if r["stratum"]=="station"]
        out=slope_model(rows)
        self.assertFalse(out["estimable"])
        self.assertIsNone(out["altitude_slope_per_100m"])
    def test_constant_altitude_blocked_and_not_falsely_signed(self):
        rows=self.synthetic()
        for r in rows:r["altitude"]=1590
        out=slope_model(rows)
        self.assertFalse(out["estimable"])
        self.assertIn("altitude",out["reasons"][0])
    def test_permutation_only_moves_altitude_within_original_station_stratum(self):
        rows=self.synthetic()
        out=source_site_permutation(rows,seed=316,draws=99)
        self.assertEqual(out["n_original_transects"],31)
        self.assertEqual(out["n_permutations"],99)
        self.assertTrue(0<out["upper_rank"]<=1)
        self.assertTrue(0<out["two_sided_rank"]<=1)
        self.assertEqual(out["whole_site_LOO_n"],31)
        self.assertEqual(out["whole_site_LOO_n_positive"],31)
        self.assertAlmostEqual(out["observed_per_100m"],0.4,places=10)
    def test_all_pairs_share_identical_taxon_ambiguity_mask(self):
        from source_matched_rail_road_v314 import audit_names,unresolved_keys
        def norm(n):
            return {"A species":"A species","B species":"B species",
                    "C species":"C species","C species agg.":None}.get(n)
        # An unresolved source name "C species agg." masks even its
        # valid counterpart, and does so in every habitat and year.
        sets={(2014,"S","Track"):{"A species","C species agg."},
              (2014,"S","Embankment"):{"B species"},
              (2014,"S","Interior Plot 1"):{"B species","C species"},
              (2014,"S","Interior Plot 2"):{"A species"}}
        out=build_normalized_plot_pairs(sets,norm)
        self.assertEqual(out[2014]["S"],({"A species","B species"},{"A species","B species"}))
    def test_numerical_thermometer_is_not_the_trait_index(self):
        one=make_pair({"A species","B species"},{"A species","C species"},
                      {"values":{"A species":2.0,"B species":8.0,"C species":4.0},
                       "statuses":{"A species":"numeric","B species":"numeric","C species":"numeric"},
                       "limits":[1,12]})
        self.assertAlmostEqual(one["difference"],2.0,places=10)
        self.assertEqual(one["n_rail_numeric"],2)
        self.assertEqual(one["n_interior_numeric"],2)
        self.assertFalse(one["estimable"])
if __name__=="__main__":unittest.main()
