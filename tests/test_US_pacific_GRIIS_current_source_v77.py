"""v7.7 fixed Pacific country GRIIS source contract tests (synthetic fixtures).

No LIVE GBIF event rows, seed flux, and no invented train/plant effects.
"""
import sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_US_pacific_GRIIS_current_source_v77 import (
    BLOCK,COUNTRY,PERIODS,MAX_STATUS_UNKNOWN_RATE,
    dataset_source_params,dataset_facets,source_cell,
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid

class FixedPacificGRIISSourceContracts(unittest.TestCase):
    def setUp(self):
        self.cell=make_grid(BLOCK)[0]
    def test_not_selected_by_biological_outcome_and_not_legacy_PRESENT_filter(self):
        self.assertEqual(BLOCK,"us_pacific")
        self.assertEqual(COUNTRY,"US")
        self.assertEqual(len(make_grid(BLOCK)),16)
        p=dataset_source_params(self.cell,PERIODS[1])
        self.assertNotIn("occurrenceStatus",p)
        self.assertEqual(p["country"],"US")
        self.assertEqual(p["facet"],["datasetKey","occurrenceStatus"])
    def test_wrong_period_or_country_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"country/year"):
            dataset_source_params(dict(self.cell,country="CA"),PERIODS[0])
        with self.assertRaisesRegex(ValueError,"country/year"):
            dataset_source_params(self.cell,(2019,2025))
    def test_real_query_contract_synthetic_facets_distinguishes_detected(self):
        result={"count":100,"facets":[
          {"field":"DATASET_KEY","counts":[
            {"name":"d1","count":55},{"name":"d2","count":30},{"name":"d3","count":15}]},
          {"field":"OCCURRENCE_STATUS","counts":[
            {"name":"DETECTED","count":99},{"name":"NOT_DETECTED","count":1}]}]}
        with patch("audit_US_pacific_GRIIS_current_source_v77.get_json",return_value=result):
            r=dataset_facets(object(),self.cell,PERIODS[1])
        self.assertEqual(r["indexed"],100)
        self.assertEqual(r["n_dataset_source_categories_at_least"],3)
        self.assertAlmostEqual(r["status_not_detected_or_unknown"],.01)
    def test_detected_source_status_and_taxonomy_are_not_presumed_true_absence(self):
        fac={"indexed":100,"n_dataset_source_categories_at_least":3,
             "top_dataset_fraction":.55,"status_not_detected_or_unknown":.01}
        full={
          "source_indexed_all_Plantae_records":100,
          "source_non_DETECTED_or_unknown_fraction":.01,
          "source_speciesKey_full_pagination":True,
          "n_source_speciesKey_facet_pages":1,
          "speciesKey_count_by_all_taxa":{"1":25,"2":25,"3":20,"4":20,"5":10}}
        with patch("audit_US_pacific_GRIIS_current_source_v77.dataset_facets",return_value=fac):
            with patch("audit_US_pacific_GRIIS_current_source_v77.whole_source_species",return_value=full):
                r=source_cell(object(),self.cell,PERIODS[1],set(full["speciesKey_count_by_all_taxa"]))
        self.assertTrue(r["source_eligible"])
        self.assertEqual(len(r["observed_GRIIS_country_specieskey_records"]),5)
        self.assertTrue(r["NOT_independent_observer_survey_or_confirmed_plant_presence"])
    def test_high_unknown_status_blocks_not_silently_accepted(self):
        fac={"indexed":100,"n_dataset_source_categories_at_least":3,
             "top_dataset_fraction":.55,"status_not_detected_or_unknown":.03}
        full={"source_indexed_all_Plantae_records":100,
              "source_non_DETECTED_or_unknown_fraction":.03,
              "source_speciesKey_full_pagination":True,
              "n_source_speciesKey_facet_pages":1,
              "speciesKey_count_by_all_taxa":{"1":25,"2":25,"3":20,"4":20,"5":10}}
        with patch("audit_US_pacific_GRIIS_current_source_v77.dataset_facets",return_value=fac):
            with patch("audit_US_pacific_GRIIS_current_source_v77.whole_source_species",return_value=full):
                r=source_cell(object(),self.cell,PERIODS[1],set(full["speciesKey_count_by_all_taxa"]))
        self.assertFalse(r["source_eligible"])
        self.assertEqual(r["exclude_reason"],"unknown_or_not_detected_source_gt2pct")
        self.assertEqual(MAX_STATUS_UNKNOWN_RATE,.02)

if __name__=="__main__":unittest.main()
