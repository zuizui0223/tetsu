"""v9.3 synthetic source QA: iNaturalist research-grade dataset dominates
GBIF indexed plant source, but same dataset != same observer/survey.
NO freight seed delivery, botanical colonization or rail-stress adaptation.
"""
import sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_US48_iNaturalist_GBIF_source_dominance_v93 import (
    INAT_GBIF_DATASET_UUID,SOURCE,WINDOWS,MIN_STATES,
    source_frame,analyse_inat_state,run)

def baseline_state():
    return {"state_name":"State 1",
        "two_year_source":{"two_period_quality":True},
        "source_windows":{
            "before_2022":{"n_all_Plantae_indexed_records":100,
              "n_US_country_GRIIS_indexed_records":45,
              "US_country_GRIIS_positive_speciesKey_record_counts":{
                   "1001":15,"1002":15,"1003":15}},
            "after_2022":{"n_all_Plantae_indexed_records":200,
              "n_US_country_GRIIS_indexed_records":60,
              "US_country_GRIIS_positive_speciesKey_record_counts":{
                   "1001":15,"1002":15,"1004":30}}}}
def provider_stages():
    k=INAT_GBIF_DATASET_UUID
    return {
     "before_2022":{
        "dataset_source":{"top_dataset_UUID_counts":{k:90},
                          "top100_source_count_coverage":.99},
        "INAT_scoped":{"original_indexed_all_Plantae_records":90,
                "n_US_country_GRIIS_indexed_records":30,
                "US_country_GRIIS_positive_speciesKey_record_counts":{
                   "1001":15,"1002":15}}},
     "after_2022":{
        "dataset_source":{"top_dataset_UUID_counts":{k:190},
                          "top100_source_count_coverage":.999},
        "INAT_scoped":{"original_indexed_all_Plantae_records":190,
                "n_US_country_GRIIS_indexed_records":45,
                "US_country_GRIIS_positive_speciesKey_record_counts":{
                   "1001":15,"1004":30}}}
    }

def full_frame():
    states={}
    census={}
    for i in range(1,52):
        code=f"{i:02d}"
        s=baseline_state()
        s["state_name"]=f"State {i}"
        s["two_year_source"]["two_period_quality"]=i<=48
        states[code]=s
        census[code]={"name":f"State {i}","abbreviation":f"ST{i}"}
    out={
        "complete":True,
        "source_periods":{"before_2022":"2018,2021","after_2022":"2023,2025"},
        "source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows":True,
        "states":states
    }
    return out,census

class SourceDominanceTests(unittest.TestCase):
    def test_official_doi_accepted_and_fixed_48states(self):
        self.assertEqual(INAT_GBIF_DATASET_UUID,"50c9509d-22c7-4a22-a47d-8c48425ef4a7")
        self.assertIn("Research-grade",SOURCE)
        self.assertEqual(WINDOWS,("before_2022","after_2022"))
        self.assertEqual(MIN_STATES,40)
    def test_48_state_fips_eligible_selection_never_by_rail_signal(self):
        raw,census=full_frame()
        chosen=source_frame(raw,census,set(str(i) for i in range(1000,1200)))
        self.assertEqual(len(chosen),48)
        self.assertNotIn("49",chosen)
        raw["source_periods"]["after_2022"]="2022,2025"
        with self.assertRaisesRegex(ValueError,"source universe"):
            source_frame(raw,census,set(str(i) for i in range(1000,1200)))
    def test_exact_non_iNat_indexed_record_count_not_species_presence(self):
        m=analyse_inat_state(baseline_state(),provider_stages())
        early=m["two_GBIF_year_source_periods"]["before_2022"]
        late=m["two_GBIF_year_source_periods"]["after_2022"]
        self.assertEqual(early["NON_iNaturalist_ALLPlantae_indexed_source_records"],10)
        self.assertEqual(late["NON_iNaturalist_ALLPlantae_indexed_source_records"],10)
        self.assertEqual(early["NON_iNaturalist_country_GRIIS_indexed_source_records"],15)
        self.assertEqual(late["NON_iNaturalist_country_GRIIS_indexed_source_records"],15)
        self.assertFalse(m["NON_iNaturalist_20_GRIIS_source_records_available_in_BOTH_windows"])
    def test_inat_late_only_counts_source_not_immigration(self):
        m=analyse_inat_state(baseline_state(),provider_stages())
        self.assertEqual(m["all_provider_late_only_country_GRIIS_recorded_taxon_categories"],1)
        self.assertEqual(m["all_provider_late_only_categories_also_recorded_late_in_INAT"],1)
        self.assertEqual(m["fraction_all_provider_late_only_categories_also_in_iNaturalist_late_source"],1.)
        self.assertTrue(m["iNaturalist_dataset_key_same_in_both_periods_NOT_same_observer_or_detection"])
        self.assertTrue(m["NOT_actual_establishment_or_new_alien_species"])
    def test_provider_scoped_taxa_must_be_subset_original_all_gbif(self):
        original=baseline_state()
        p=provider_stages()
        p["after_2022"]["INAT_scoped"]["US_country_GRIIS_positive_speciesKey_record_counts"]["9999"]=10
        with self.assertRaisesRegex(ValueError,"not subset"):
            analyse_inat_state(original,p)
    def test_same_source_by_uuid_not_invented_as_independent_field_survey(self):
        original,census=full_frame()
        inat=INAT_GBIF_DATASET_UUID
        def fake_provider(_sess,name,win,total):
            return {"top_dataset_UUID_counts":{inat:90 if win=="before_2022" else 190},
                    "top100_source_count_coverage":1.0}
        def fake_griis(_sess,name,win,uuid,expected,keys):
            return provider_stages()[win]["INAT_scoped"]
        with patch("audit_US48_iNaturalist_GBIF_source_dominance_v93.provider_facet",
                   side_effect=fake_provider), \
             patch("audit_US48_iNaturalist_GBIF_source_dominance_v93.provider_GRIIS_facets",
                   side_effect=fake_griis):
            result=run(None,original,census,set(str(i) for i in range(1000,1200)))
        self.assertTrue(result["complete"])
        self.assertEqual(result["summary"]["n_original_48state_GBIF_provider_scoped_source_completed"],48)
        self.assertEqual(result["summary"]["n_states_iNaturalist_share_increased"],48)
        self.assertTrue(result["NO_riding_train_seed_causality_mowing_gust_or_adaptation"])
if __name__=="__main__":unittest.main()
