"""v14.4 two-track evidence cohort gate tests, source manuscript units."""
from __future__ import annotations
import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_tetsu_two_track_original_study_eligibility_v144 import (
    inventory,result
)
SOURCE=Path(__file__).resolve().parents[1]/"data"/"literature"/"tetsu_two_track_original_cohorts_v144.csv"


class GenuineStudyOutcomeGate(unittest.TestCase):
    def test_original_source_papers_stages_and_reanalysis_not_false_combined_meta(self):
        raw=inventory(SOURCE)
        r=result(raw)
        self.assertEqual(r["n_unique_published_source_study_documents"],17)
        self.assertEqual(r["n_independent_countries_with_TRUE_ready_B_incidence_matrix"],2)
        self.assertEqual(r["n_effect_sizes_meta_harmonized_for_rail_vs_nonrail_or_rail_vs_road_with_same_trait_method"],0)
        self.assertEqual(r["n_actual_direct_rail_freight_cargo_grain_x_same_site_viable_germination_x_adult_matched_cohorts"],0)
        self.assertEqual(set(r["primary_B_source_true_plot_by_taxon_matrices_obtained"]),
                         {"SI_Sabic_2026","CZ_Kutlvasr_2024"})
        self.assertFalse(r["source_2012_France_real_71site_break_vegetation_raw_matrix_acquired"])

    def test_non_independent_Swiss_molecular_related_cohorts_cannot_be_invented_as_replicates(self):
        r=result(inventory(SOURCE))
        groups=r["SAME_geography_or_overlapping_source_cohorts_MUST_not_count_as_independent_experiments"]
        self.assertIn(["CH_Schoenenberger_2012","CH_Schulze_2014_port","CH_BaselCanton_2020"],groups)
        self.assertIn(["FR_Penone_2012","FR_Blanchet_2015"],groups)

    def test_France_RATP_394_source_original_plots_and_170_distinct_alerts_NOT_18th_study(self):
        studies=inventory(SOURCE)
        source={r["study_key"]:r for r in studies}
        self.assertEqual(len(studies),17)
        r=source["FR_Rakotoson_2025"]
        self.assertEqual(r["track"],"B_D_MANAGEMENT_CANDIDATE")
        self.assertEqual(int(r["original_independent_units"]),394)
        self.assertIn("170 DISTINCT staff alerts",r["original_field_unit"])
        self.assertEqual(r["true_site_taxon_matrix_obtained"],"false")
        self.assertEqual(r["true_site_coordinates_verified"],"false")
        self.assertEqual(r["usable_for_common_species_beta"],"false")
        self.assertNotIn("FR_Rakotoson_2025",
                         result(studies)["primary_B_source_true_plot_by_taxon_matrices_obtained"])

    def test_France_RATP_170_reported_alerts_are_not_proven_mowing_or_observed_species_matrix(self):
        s=SOURCE.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"fake.csv"
            p.write_text(s.replace(
                 "FR_Rakotoson_2025,10.1111/wre.70025,FR,B_D_MANAGEMENT_CANDIDATE",
                 "FR_Rakotoson_2025,10.1111/wre.70025,FR,B"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"RATP original 394"):
                inventory(p)
            s2=s.replace(
                 "200m to alert requests NOT actual mechanical management",
                 "200m to actual randomized mechanical mowing treatments")
            p.write_text(s2,encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"RATP original 394"):
                inventory(p)

    def test_treating_french_published_results_as_retrieved_species_matrix_fails(self):
        s=SOURCE.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"fake.csv"
            p.write_text(s.replace(
                 "FR_Penone_2012,10.1016/j.biocon.2012.01.041,FR,B_CANDIDATE,railway_verge_fixed_area_plot,71,2012,71 original sites 5x1m2 quadrats on two Paris suburban railway lines,false",
                 "FR_Penone_2012,10.1016/j.biocon.2012.01.041,FR,B_CANDIDATE,railway_verge_fixed_area_plot,71,2012,71 original sites 5x1m2 quadrats on two Paris suburban railway lines,true"))
            with self.assertRaisesRegex(ValueError,"2012 original 71-site France"):
                inventory(p)

    def test_bogus_common_railroad_study_list_rejected_if_site_taxon_not_obtained(self):
        s=SOURCE.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"bad.csv"
            p.write_text(s.replace(
                "FR_Penone_2012,10.1016/j.biocon.2012.01.041,FR,B_CANDIDATE,railway_verge_fixed_area_plot",
                "FR_Penone_2012,10.1016/j.biocon.2012.01.041,FR,B,railway_verge_fixed_area_plot"))
            with self.assertRaisesRegex(ValueError,"only two currently verified"):
                inventory(p)

    def test_source_record_cannot_count_czech_unverified_station_GPS(self):
        s=SOURCE.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"bad.csv"
            p.write_text(s.replace(
                "CZ_Kutlvasr_2024,10.23855/preslia.2024.247,CZ,B,railway_station_assemblage,39,2020-2022,39 whole unequal area railway stations nested 5 source stretches,true,false,true",
                "CZ_Kutlvasr_2024,10.23855/preslia.2024.247,CZ,B,railway_station_assemblage,39,2020-2022,39 whole unequal area railway stations nested 5 source stretches,true,true,true"))
            with self.assertRaisesRegex(ValueError,"No verified original Czech"):
                inventory(p)


if __name__=="__main__":unittest.main()
