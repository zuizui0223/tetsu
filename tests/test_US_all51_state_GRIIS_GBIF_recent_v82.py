"""Synthetic fixed 51-state GBIF GRIIS source contract tests, not plant effects."""
import sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from harvest_US_all51_state_GRIIS_GBIF_recent_v82 import run,PERIOD,MIN_SOURCE_STATES

def census():
    return {f"{n:02d}":{"name":f"State {n}","abbreviation":f"S{n}"}
            for n in range(1,52)}
def baseline():
    return {"complete":True,"FAF_release":"FAF6.0 / 2022 BASE YEAR / STATE OD CSV",
            "data":{"n_distinct_state_code_categories_in_source":51}}
def source_result(*_):
    return {"source_complete":True,"source_eligible":True,
            "source_GRIIS_speciesKey_indexed_counts":{"1001":10},
            "n_source_GRIIS_country_introduced_species_categories":1}

class AllStateSourceBounds(unittest.TestCase):
    def test_exact_51_states_recent_only_and_minimum_source_support(self):
        self.assertEqual(PERIOD,(2018,2025))
        self.assertEqual(MIN_SOURCE_STATES,40)
    def test_complete_51_state_country_status_cannot_be_inferred_as_state_introduced(self):
        with patch("harvest_US_all51_state_GRIIS_GBIF_recent_v82.valid_griis",
                   return_value={"US":set(str(x) for x in range(1000,1200))}), \
             patch("harvest_US_all51_state_GRIIS_GBIF_recent_v82.get_state_period",
                   side_effect=source_result), \
             patch("harvest_US_all51_state_GRIIS_GBIF_recent_v82.time.sleep"):
            result=run({},census(),None,baseline())
        self.assertTrue(result["complete"])
        self.assertEqual(result["summary"]["n_Census_FIPS_51_state_source_completed"],51)
        self.assertEqual(result["summary"]["n_distinct_US_GRIIS_country_list_keys_observed_2018_2025"],1)
        self.assertTrue(result["no_colonization_causality_native_status_or_seed_transport_assigned"])
    def test_source_less_than_40_prevents_ecological_link_readiness(self):
        def results(_,name,period,keys):
            i=int(name.split()[-1])
            return {"source_complete":True,"source_eligible":i<=35}
        with patch("harvest_US_all51_state_GRIIS_GBIF_recent_v82.valid_griis",
                   return_value={"US":set(str(x) for x in range(1000,1200))}), \
             patch("harvest_US_all51_state_GRIIS_GBIF_recent_v82.get_state_period",
                   side_effect=results), \
             patch("harvest_US_all51_state_GRIIS_GBIF_recent_v82.time.sleep"):
            data=run({},census(),None,baseline())
        self.assertTrue(data["complete"])
        self.assertFalse(data["source_quality_minimum_for_freight_network_comparison"])
        self.assertEqual(data["summary"]["n_state_FIPS_with_adequate_2018_2025_GRIIS_plant_records"],35)
    def test_wrong_FA6_country_frame_fails_before_GBIF(self):
        bad=baseline()
        bad["data"]["n_distinct_state_code_categories_in_source"]=50
        with self.assertRaisesRegex(ValueError,"FAF6"):
            run({},census(),None,bad)
    def test_missing_state_not_impersonated_as_absence(self):
        with self.assertRaisesRegex(ValueError,"Census"):
            run({},{"01":{"name":"State 1"}},None,baseline())
if __name__=="__main__":unittest.main()
