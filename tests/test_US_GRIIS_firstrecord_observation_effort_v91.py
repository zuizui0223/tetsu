"""v9.1 synthetic safeguards: GBIF 2026 indexed record opportunity
annual rate and source taxon first recorded. NOT new colonization,
freight seed carriage or any railway stress adaptation.
"""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_US_GRIIS_firstrecord_observation_effort_v91 import (
    MIN_STATES,PERIOD_YEARS,source_state_rows,descriptive,spearman_source)

def source_fixture():
    states={}
    for i in range(1,52):
        p=200+i
        q=350+i*3
        gb={
            "before_2022":{"source_eligible":True,
                "n_all_Plantae_indexed_records":p,
                "n_US_country_GRIIS_indexed_records":50+i},
            "after_2022":{"source_eligible":True,
                "n_all_Plantae_indexed_records":q,
                "n_US_country_GRIIS_indexed_records":100+i*3,
                "n_US_country_GRIIS_observed_speciesKey_categories":50+i}
        }
        states[f"{i:02d}"]={
            "source_windows":gb,
            "two_year_source":{
                "two_period_quality":i<=48,
                "n_early_only_source_taxon_keys":5+(i%4),
                "n_late_only_first_recorded_source_taxon_keys":15+(i%6)}}
    return {"complete":True,"source_periods":{
        "before_2022":"2018,2021","after_2022":"2023,2025"},
        "source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows":True,
        "states":states}

class RecordingProcessQA(unittest.TestCase):
    def test_source_early_four_years_vs_late_three_and_no_2022(self):
        self.assertEqual(PERIOD_YEARS,{"before_2022":4,"after_2022":3})
        self.assertEqual(MIN_STATES,40)
    def test_all_48_source_qualified_states_and_three_unknown_not_zero(self):
        rows,missing=source_state_rows(source_fixture())
        self.assertEqual(len(rows),48)
        self.assertEqual(missing,["49","50","51"])
        self.assertAlmostEqual(rows[0]["post_before_all_Plantae_annual_record_ratio"],
                               ((350+3)/3)/((200+1)/4))
        self.assertTrue(rows[0]["post_only_first_recorded_taxon_categories_NOT_immigrations"]>0)
    def test_source_annual_rate_increase_is_descriptive_not_plant_establishment(self):
        rows,missing=source_state_rows(source_fixture())
        x=descriptive(rows,missing)
        self.assertEqual(x["n_qualified_Census_US_states_for_both_year_windows"],48)
        self.assertEqual(x["annual_all_Plantae_source_record_ratio_post_over_pre"][
            "n_states_more_post_records_per_year"],48)
        self.assertEqual(x["n_states_more_post_only_recorded_taxon_keys_than_pre_only"][
            "post_only_greater"],48)
        self.assertTrue(x["NO_invasion_seed_railroad_mowing_wind_adaptation_inferred"])
    def test_2022_year_overlap_rejected(self):
        x=source_fixture()
        x["source_periods"]["before_2022"]="2018,2022"
        with self.assertRaisesRegex(ValueError,"source windows"):
            source_state_rows(x)
    def test_ineligible_but_missing_source_cannot_be_imputed_zero(self):
        x=source_fixture()
        x["states"]["01"]["two_year_source"]["two_period_quality"]=False
        x["states"]["01"]["source_windows"]["before_2022"]={}
        rows,miss=source_state_rows(x)
        self.assertEqual(len(rows),47)
        self.assertIn("01",miss)
        self.assertNotIn("01",[v["Census_FIPS_state_code"] for v in rows])
    def test_late_only_source_category_can_never_exceed_late_recorded_taxa(self):
        x=source_fixture()
        x["states"]["01"]["two_year_source"][
            "n_late_only_first_recorded_source_taxon_keys"]=10000
        with self.assertRaisesRegex(ValueError,"malformed"):
            source_state_rows(x)
    def test_source_spearman_is_undefined_not_significant_if_no_variance(self):
        self.assertIsNone(spearman_source([1]*48,list(range(48))))
        self.assertIsNone(spearman_source([1,2],[2,1]))
if __name__=="__main__":unittest.main()
