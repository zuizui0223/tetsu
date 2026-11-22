"""v320 original complete French-Belgian railway vs matrix source support guards.

No external source in the repository. Real 420-page report must be fetched in
dedicated GitHub Actions and checked by the source parser.
"""
import json,unittest
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"data"/"derived"/"tetsu_v320_FranceBelgium_Likoto_original6rail_sites_matched_plant_source_support_summary.json"

class FrenchBelgianOriginalReportScope(unittest.TestCase):
    def setUp(self):
        self.j=json.loads(SOURCE.read_text(encoding="utf8"))
    def test_original_six_railway_candidates_are_one_study_not_six_countries(self):
        r=self.j
        self.assertEqual(r["source_original_pdf_pages"],420)
        self.assertEqual(r["source_original_24_site_infrastructure_inventory_locations"],24)
        self.assertEqual(r["original_railway_location_candidates_from_TOC"],6)
        self.assertEqual(r["original_railway_high_speed_locations"],3)
        self.assertEqual(r["original_railway_conventional_locations"],3)
        self.assertEqual(sorted(x["code"] for x in r["source_original_railway_site_list"]),[17,18,19,20,21,24])
        self.assertEqual(r["n_original_published_studies_added"],0)
    def test_matched_botanical_plots_can_exist_but_original_species_matrix_is_unverified(self):
        r=self.j
        self.assertTrue("where possible" in r["source_botany_sampling_note"])
        self.assertEqual(r["n_source_verified_same_area_railway_vs_local_matrix_species_biological_plot_pairs"],0)
        self.assertEqual(r["n_source_verified_railway_species_by_plot_incidence_cells"],0)
        self.assertTrue(any("NOT zero" in a for a in r["interpretation_limits"]))
        self.assertTrue(r["no_new_railway_climate_or_transport_effect_estimated"])
    def test_source_pages_and_legend_are_traceable(self):
        r=self.j
        self.assertEqual(r["n_original_railway_locations_with_text_extracted_infrastructure_matrix_map_legend"],5)
        for item in r["source_original_railway_site_list"]:
            lo,hi=item["chapter_pages_pdf_1based"]
            self.assertTrue(100<lo<=hi<=345)
            self.assertTrue(item["author_la_flore_narrative_pages"])
            self.assertTrue(all(lo<=p<=hi for p in item["author_la_flore_narrative_pages"]))
            self.assertTrue(all(lo<=p<=hi for p in item["botanist_infrastructure_matrix_map_legend_pages"]))
if __name__=="__main__":unittest.main()
