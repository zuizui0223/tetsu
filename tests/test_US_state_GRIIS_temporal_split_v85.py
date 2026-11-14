"""v8.5 synthetic GBIF source-year separation around single 2022 FAF6
crop freight base-year. No state-introduction/seed inference allowed.
"""
import sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from harvest_US_state_GRIIS_temporal_split_v85 import (
    PERIODS,MIN_BOTH_PERIOD_QUALIFIED_STATE_SOURCES,query_state_years,
    get_window,audit_temporal_window,run)

def state_source():
    return {f"{n:02d}":{"name":f"State {n}","abbreviation":f"S{n}"}
            for n in range(1,52)}

def document():
    return {"complete":True,"source_eligible":True,
            "US_country_GRIIS_positive_speciesKey_record_counts":{
                "1001":15,"1002":15,"1003":15,"1004":15,"1005":15},
            "n_all_Plantae_indexed_records":200}

def response(non=0):
    return {"count":200,"facets":[
        {"field":"SPECIES_KEY","counts":[{"name":str(i),"count":25}
                    for i in range(1001,1006)]},
        {"field":"DATASET_KEY","counts":[{"name":k,"count":v}
                    for k,v in (("d1",100),("d2",50),("d3",50))]},
        {"field":"OCCURRENCE_STATUS","counts":[
                    {"name":"DETECTED","count":200-non},
                    {"name":"NOT_DETECTED","count":non}]}
    ]}

def previous(states):
    return ({"complete":True,"FAF_release":"FAF6.0 / 2022 BASE YEAR / STATE OD CSV"},
            {"complete":True,"n_requested_state_fips_codes":51,
             "states":{c:{} for c in states}})

class VerifiedSourceWindowSafety(unittest.TestCase):
    def test_2022_is_excluded_and_no_future_2026_year(self):
        self.assertEqual(PERIODS,{"before_2022":"2018,2021","after_2022":"2023,2025"})
        self.assertEqual(MIN_BOTH_PERIOD_QUALIFIED_STATE_SOURCES,40)
        self.assertEqual(query_state_years("California","after_2022")["year"],"2023,2025")
        self.assertNotIn("occurrenceStatus",query_state_years("California","before_2022"))
        with self.assertRaisesRegex(ValueError,"FROZEN"):
            query_state_years("California","2022")
    def test_indexed_source_facet_is_only_recorded_category(self):
        with patch("harvest_US_state_GRIIS_temporal_split_v85.source_get",
                   return_value=response(non=0)):
            x=get_window(None,"California","after_2022",set(str(i) for i in range(1001,1006)))
        self.assertTrue(x["source_eligible"])
        self.assertEqual(x["n_US_country_GRIIS_observed_speciesKey_categories"],5)
        self.assertTrue(x["NOT_species_absence_nor_establishment_nor_state_introduced_status"])
    def test_bad_undetected_status_not_called_established_presence(self):
        with patch("harvest_US_state_GRIIS_temporal_split_v85.source_get",
                   return_value=response(non=12)):
            x=get_window(None,"California","before_2022",set(str(i) for i in range(1001,1006)))
        self.assertFalse(x["source_eligible"])
        self.assertIn("source_NOT_DETECTED_or_unknown_status_above_2pct",
                      x["source_eligibility_failures"])
    def test_late_only_source_is_not_colonization(self):
        early=document();late=document()
        early["US_country_GRIIS_positive_speciesKey_record_counts"].pop("1005")
        x=audit_temporal_window(early,late)
        self.assertEqual(x["n_late_only_first_recorded_source_taxon_keys"],1)
        self.assertTrue(x["n_late_only_species_recorded_is_NOT_biological_introduction"])
    def test_original_51_state_geography_and_2022_base_year_contract(self):
        states=state_source();fa,fl=previous(states)
        with patch("harvest_US_state_GRIIS_temporal_split_v85.valid_griis",
                   return_value={"US":set(str(i) for i in range(1000,1250))}), \
             patch("harvest_US_state_GRIIS_temporal_split_v85.get_window",
                   return_value=document()), \
             patch("harvest_US_state_GRIIS_temporal_split_v85.time.sleep"):
            out=run({},states,None,fa,fl)
        self.assertTrue(out["complete"])
        self.assertEqual(out["summary"]["n_both_window_GRIIS_US_recording_source_quality_states"],51)
        self.assertTrue(out["source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows"])
        self.assertTrue(out["not_seed_delivery_causation_nor_mowing_wind_adaptation"])
    def test_inadequate_two_window_state_support_is_reported_not_imputed(self):
        states=state_source();fa,fl=previous(states)
        def uneven(_session,name,window,_keys):
            code=int(name.split()[-1])
            d=document()
            if code>29 and window=="before_2022":
                d["source_eligible"]=False
            return d
        with patch("harvest_US_state_GRIIS_temporal_split_v85.valid_griis",
                   return_value={"US":set(str(i) for i in range(1000,1250))}), \
             patch("harvest_US_state_GRIIS_temporal_split_v85.get_window",
                   side_effect=uneven), \
             patch("harvest_US_state_GRIIS_temporal_split_v85.time.sleep"):
            out=run({},states,None,fa,fl)
        self.assertEqual(out["summary"]["n_both_window_GRIIS_US_recording_source_quality_states"],29)
        self.assertFalse(out["n_qualified_source_states_gate_40"])
        self.assertTrue(out["complete"])
if __name__=="__main__":unittest.main()
