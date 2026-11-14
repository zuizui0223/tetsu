"""v13.1 source verified Basel Swiss govt 8-year two-site OSR monitoring tests."""
from __future__ import annotations
import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Swiss_Basel_rail_freight_port_GM_rapeseed_2013_2020_v131 import (
    ORIGINAL,YEARS,STATION,PORT,
    parse_original_government_source_csv,summarize_source,verify_government_original_pdf,
)

SRC=(Path(__file__).resolve().parents[1]/"data"/
    "Swiss_Basel_CantonalLab_2020_GM_rapeseed_station_port_2013_2020_v131.csv")


class BaselCantonGovernment2013To2020SourceTests(unittest.TestCase):
    def test_original_16_source_site_year_population_and_independently_tested_plant_counts(self):
        rows=parse_original_government_source_csv(SRC)
        self.assertEqual(len(rows),16)
        out=summarize_source(rows)
        self.assertEqual(out["sum_of_original_number_of_PCR_tested_plants_2013_2020"],7908)
        self.assertEqual(out["government_report_original_GM_positive_plant_tests_2013_2020"],459)
        self.assertAlmostEqual(out["government_report_GM_positive_fraction_of_PCR_tests"],459/7908)
        self.assertEqual(out["n_true_original_monitoring_field_sites"],2)
        self.assertEqual(out["n_original_site_year_observation_units"],16)
        self.assertTrue(out["site_contrast_not_randomized_controlled_rail_cargo_experiment"])

    def test_st_johann_zero_all_plants_2017_but_2018_2020_nonzero_not_gm(self):
        out=summarize_source(parse_original_government_source_csv(SRC))
        s=out["source_site_year_time_series"][STATION]
        self.assertEqual(s["source_start_estimated_feral_population_2013"],639)
        self.assertEqual(s["source_end_estimated_feral_population_2020"],2)
        self.assertEqual(s["source_years_zero_estimated_feral_plants"],[2017])
        self.assertEqual(s["n_original_site_years_any_feral_plants"],7)
        self.assertEqual(s["original_source_genetic_test_count_2013_2020_sum_of_tests_NOT_INDEPENDENT_plants"],198)
        self.assertEqual(out["station_report_observed_GT73_GM_no_longer_identified_since_year"],2016)

    def test_port_nonzero_all_eight_years_while_gm_source_held_separate_from_population(self):
        out=summarize_source(parse_original_government_source_csv(SRC))
        p=out["source_site_year_time_series"][PORT]
        self.assertEqual(p["n_original_site_years_any_feral_plants"],8)
        self.assertEqual(p["original_peak_estimated_feral_plants"],9247)
        self.assertEqual(p["original_year_peak_estimated_plants"],2014)
        self.assertEqual(p["source_end_estimated_feral_population_2020"],752)
        self.assertEqual(p["original_source_genetic_test_count_2013_2020_sum_of_tests_NOT_INDEPENDENT_plants"],7710)
        self.assertTrue(out["port_GM_plant_events_present_through_2020"])

    def test_source_population_vs_tested_count_corruption_refused(self):
        old=SRC.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"official-csv.csv"
            f.write_text(old.replace("St_Johann_freight_station,2013,639,100",
                                     "St_Johann_freight_station,2013,639,101"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"published Table1/2"):
                parse_original_government_source_csv(f)
            f.write_text(old.replace("Kleinhuningen_Rhine_port,2014,9247,1333",
                                     "Kleinhuningen_Rhine_port,2014,9248,1333"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"published Table1/2"):
                parse_original_government_source_csv(f)

    def test_original_official_pdf_source_wont_accept_a_placeholder(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"fake.pdf"
            f.write_text("%PDF-fake header\n")
            with self.assertRaisesRegex(ValueError,"not valid substantive PDF"):
                verify_government_original_pdf(f)


if __name__=="__main__":unittest.main()
