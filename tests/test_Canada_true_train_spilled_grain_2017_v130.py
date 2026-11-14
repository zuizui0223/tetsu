"""v13.0 authentic 2017 Canadian train-deposited grain author-summary source tests.

Genuine original 19 site × 3 year railway cargo-kernel traps from
Gangadharan et al. (2017) DOI 10.1111/acv.12336.
No fabricated original 1979 site/visit rows or seed viability.
"""
from __future__ import annotations
import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Canada_true_train_spilled_grain_2017_v130 import (
    run,read_published_author_rail_kernels_summary,
    actual_original_annual_source
)

CSV=Path(__file__).resolve().parents[1]/"data"/"Canada_Gangadharan2017_original_train_grain_seasonal_summary_v130.csv"


class OriginalTrainCargoSeedCollection2017Tests(unittest.TestCase):
    def test_all_original_real_train_author_published_seasonal_summary_values(self):
        rows=read_published_author_rail_kernels_summary(CSV)
        self.assertEqual(len(rows),4)
        r=run(CSV)
        self.assertTrue(r["complete"])
        x=r["original_author_real_train_GRAIN_deposition_source"]
        self.assertEqual(x["n_original_real_railway_deposition_site_units"],19)
        self.assertEqual(x["n_original_track_distance_km"],134)
        self.assertEqual(x["reported_original_field_sample_count"],1979)
        self.assertEqual(sum(x["original_sampling_collector_method_counts"].values()),1979)
        self.assertEqual(x["author_rounded_annual_TOTAL_model_projection_tonnes_134km"],110)
        self.assertAlmostEqual(x[
            "mean_warm_projection_plus_one_observed_winter_projection_tonnes_134km"],109.4)
        self.assertEqual(x["author_rounded_grain_type_share_total_percent"],101)
        self.assertFalse(x["directly_measured_grain_seed_viability_or_germination"])
        self.assertFalse(x["directly_measured_individual_spilled_grain_establishing_in_ballast"])

    def test_original_2015_higher_within_same_sampler_warm_months_but_not_trend_significance(self):
        x=run(CSV)["original_author_real_train_GRAIN_deposition_source"]
        self.assertAlmostEqual(x["source_2015_warm_rate_divided_by_2014_warm_rate"],
                               2.75/.87)
        self.assertAlmostEqual(x[
            "source_2015_warm_masses_divided_by_2014_model_extrapolated_masses"],
            112.8/35.6)
        # Source pooled 1.64 is not a simple mean of 3 year point estimates.
        self.assertAlmostEqual(x[
            "source_unweighted_mean_of_3_reported_yearwise_warm_rates_NOT_POOLED"],
            (1.12+.87+2.75)/3)
        self.assertNotAlmostEqual(x["source_all_warm_season_pooled_source_mean_g_m2_day"]
                                  if "source_all_warm_season_pooled_source_mean_g_m2_day" in x
                                  else x["original_all_warm_season_pooled_source_mean_g_m2_day"],
                                  (1.12+.87+2.75)/3)

    def test_published_reported_train_kernel_value_changed_is_not_accepted(self):
        raw=CSV.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"train.csv"
            p.write_text(raw.replace("April-October,2014,2014,basket,0.87",
                                     "April-October,2014,2014,basket,0.88"))
            with self.assertRaisesRegex(ValueError,"summary changed"):
                read_published_author_rail_kernels_summary(p)

    def test_author_winter_SCREEN_is_not_falsely_harmonized_into_basket(self):
        x=run(CSV)
        years=x["author_original_season_by_year_and_method_results"]
        self.assertEqual(years[3]["original_collector"],"screen")
        self.assertEqual(years[0]["original_collector"],"basket")
        self.assertEqual(years[3]["author_reported_mean_g_per_m2_day"],1.52)
        self.assertEqual(years[3]["author_extrapolated_tonnes_134km_NOT_physically_recovered"],44.6)

    def test_missing_original_source_year_not_silent_zero(self):
        rows=CSV.read_text(encoding="utf-8").splitlines()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"train.csv"
            p.write_text("\n".join(rows[:-1])+"\n",encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Missing published"):
                read_published_author_rail_kernels_summary(p)

if __name__=="__main__":unittest.main()
