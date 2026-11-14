"""v12.8 original 60 observed Austrian feral Brassica napus field sites.

Tests use original published numeric Table14 manually transcribed and
independently reconciled with source 2016 government report's 60/44
site counts and 2014/2015 1001+1112 sampled plant totals. No real
train seed trap counts are represented or tested.
"""
from __future__ import annotations
import csv,math,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_Austria_ferreral_osr_actual_sites_v128 import (
    validate_original_table,source_presence,observed_count_summary,
    conditional_stratified_exact,fisher_exact_two_sided,
    risk_and_source_sensitivity,run,
)

CSV=Path(__file__).resolve().parents[1]/"data"/"austria_pascher_2016_feral_osr_table14_manual_transcription_v128.csv"


class FeralAustrianSourcePlantTests(unittest.TestCase):
    def test_true_author_table14_exact_site_ids_years_and_44_positive(self):
        rows=validate_original_table(CSV)
        self.assertEqual(len(rows),60)
        stats=observed_count_summary(rows)
        self.assertEqual(stats["n_original_author_field_sites"],60)
        self.assertEqual(stats["n_site_sources_with_feral_OSR_sampled_in_either_year"],44)
        self.assertEqual(stats["original_2014_sampled_plant_individuals_NOT_plant_census"],1001)
        self.assertEqual(stats["original_2015_sampled_plant_individuals_NOT_plant_census"],1112)
        self.assertEqual(stats["original_author_total_sampled_plant_individuals_NOT_seeds"],2113)
        self.assertEqual(stats["n_original_sites_with_plant_sampling_in_BOTH_years"],33)
        self.assertEqual(stats["n_original_sites_with_plant_sampled_ONLY_2014"],3)
        self.assertEqual(stats["n_original_sites_with_plant_sampled_ONLY_2015"],8)
        self.assertEqual(stats["n_original_sites_with_no_sampled_feral_OSR_in_EITHER_year"],16)

    def test_actual_cultivation_area_strata_preserve_ordinary_station_vs_road_only(self):
        rows=validate_original_table(CSV)
        self.assertEqual(source_presence(rows,"BHR"),(10,10))
        self.assertEqual(source_presence(rows,"SR"),(8,11))
        self.assertEqual(source_presence(rows,"BHK"),(6,10))
        self.assertEqual(source_presence(rows,"SK"),(4,11))
        result=risk_and_source_sensitivity(rows)
        self.assertAlmostEqual(result["original_station_minus_road_observed_detection_fraction"],
                               16/20-12/22)
        self.assertEqual(result["n_original_non_hotspot_station_source_sites"],20)
        self.assertEqual(result["n_original_non_hotspot_road_source_sections"],22)
        self.assertEqual(result["n_source_station_sites_with_OSR_samples"],16)
        self.assertEqual(result["n_source_road_sections_with_OSR_samples"],12)
        self.assertEqual(result["original_42_site_leave_one_out_n"],42)
        self.assertEqual(result["original_42_site_leave_one_out_opposite_sign_vs_original"],0)

    def test_exact_stratified_hypergeometric_probability_normalized_and_source_p(self):
        stat=conditional_stratified_exact([(10,10,8,11),(6,10,4,11)])
        self.assertAlmostEqual(
            stat["exact_two_sided_centered_stratified_tail_probability"],
            0.08251251704602575,places=12)
        self.assertAlmostEqual(stat["source_station_successes_exchangeability_expected"],
                               10*18/21+10*10/21)
        self.assertTrue(stat["NOT_randomized_railway_seed_carriage_treatment_test"])
        self.assertAlmostEqual(fisher_exact_two_sided(16,4,12,10),
                               0.10804559816334419,places=12)
        self.assertAlmostEqual(fisher_exact_two_sided(10,0,8,3),
                               .2142857142857143,places=12)

    def test_published_table14_source_missing_or_wrong_sample_count_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"bad.csv"
            source=CSV.read_text(encoding="utf-8")
            p.write_text(source.replace("BHR1,26,47","BHR1,26,48"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"original Table14 sample counts"):
                validate_original_table(p)
            p.write_text(source.replace("GBH7,1,2","GBH6,1,2"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"site codes"):
                validate_original_table(p)
            p.write_text(source.replace("BHR2,11,3","BHR1,11,3"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"duplicate"):
                validate_original_table(p)

    def test_reconcile_original_source_study_not_real_seed_transport(self):
        res=run(CSV)
        self.assertEqual(res["version"],"v12.8")
        self.assertEqual(res["original_author_60_site_2yr_sample_validation"][
            "n_site_sources_with_feral_OSR_sampled_in_either_year"],44)
        self.assertFalse(res["integration_with_previous_tetsu_research"][
            "arrival_measured_directly_by_Pascher2016"])
        self.assertTrue(res["original_site_count_not_total_species_flora_or_cargo_seeds"])


if __name__=="__main__":
    unittest.main()
