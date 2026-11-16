"""Source-only Census-FIPS GBIF USA stateProvince GRIIS taxonomy QA tests."""
import sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from harvest_US_state_GRIIS_GBIF_plant_source_v81 import (
    PILOT_FIPS,SOURCE_PERIODS,parse_census_state_table,
    query_params,get_state_period,MAX_NON_DETECTED_OR_UNKNOWN)

def census_fixture():
    text="STATE|STUSAB|STATE_NAME|STATENS\n"
    codes=["01","02","04","05","06","08","09","10","11","12","13","15",
           "16","17","18","19","20","21","22","23","24","25","26",
           "27","28","29","30","31","32","33","34","35","36","37",
           "38","39","40","41","42","44","45","46","47","48","49",
           "50","51","53","54","55","56"]
    for c in codes:
        text+=f"{c}|S{c}|State{c}|0000000{c}\n"
    text+="60|AS|American Samoa|001\n"
    return text

def gbif_fixture(non=1):
    return {"count":200,"facets":[
        {"field":"speciesKey","counts":[
            {"name":str(i),"count":25} for i in (1001,1002,1003,1004,1005)]},
        {"field":"occurrenceStatus","counts":[
            {"name":"DETECTED","count":200-non},
            {"name":"NOT_DETECTED","count":non}]},
        {"field":"datasetKey","counts":[
            {"name":"uuidA","count":90},{"name":"uuidB","count":60},
            {"name":"uuidC","count":50}]},
    ]}

class FrozenStateSourceTests(unittest.TestCase):
    def test_six_or_more_geographic_pilot_not_selected_from_freight_effect(self):
        self.assertEqual(PILOT_FIPS,("06","12","17","19","36","37","48","53"))
        self.assertEqual(SOURCE_PERIODS,((2010,2017),(2018,2025)))
    def test_official_Census_FIPS_51_state_plus_DC_source_contract(self):
        x=parse_census_state_table(census_fixture())
        self.assertEqual(len(x),51)
        self.assertNotIn("60",x)
        self.assertIn("06",x)
        self.assertEqual(x["06"]["name"],"State06")
    def test_gbif_state_query_no_legacy_broken_present_enum(self):
        p=query_params("California",(2018,2025))
        self.assertEqual(p["stateProvince"],"California")
        self.assertNotIn("occurrenceStatus",p)
        self.assertEqual(p["facet"],["speciesKey","occurrenceStatus","datasetKey"])
    def test_actual_GRIIS_source_taxa_counts_only_not_true_state_introduction(self):
        with patch("harvest_US_state_GRIIS_GBIF_plant_source_v81.source_get",return_value=gbif_fixture(1)):
            x=get_state_period(None,"California",(2018,2025),
                               {"1001","1002","1003","1004","1005"})
        self.assertTrue(x["source_eligible"])
        self.assertEqual(x["n_source_GRIIS_country_introduced_species_categories"],5)
        self.assertEqual(x["n_source_GRIIS_country_introduced_species_records"],125)
        self.assertTrue(x["GRIIS_country_level_introduced_DOES_NOT_imply_state_level_introduced"])
    def test_high_notdetected_rate_source_rejected_not_species_absent(self):
        with patch("harvest_US_state_GRIIS_GBIF_plant_source_v81.source_get",return_value=gbif_fixture(10)):
            x=get_state_period(None,"California",(2018,2025),
                               {"1001","1002","1003","1004","1005"})
        self.assertFalse(x["source_eligible"])
        self.assertIn("too_many_NOT_DETECTED_or_UNKNOWN_status",x["reason_not_qualified"])
        self.assertEqual(MAX_NON_DETECTED_OR_UNKNOWN,.02)
    def test_invalid_period_state_query_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"frozen year"):
            query_params("California",(2019,2025))
        with self.assertRaisesRegex(ValueError,"frozen year"):
            query_params(" ",(2018,2025))
if __name__=="__main__":unittest.main()
