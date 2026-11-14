"""v9.2 synthetic GBIF datasetKey-by-period stability source QA.
Never infer independent same-observer field surveys or new invasions.
"""
import sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_US_GBIF_matched_provider_temporal_v92 import (
    PILOT_FIPS,PERIODS,N_PROVIDER,
    common_provider_params,provider_taxon_params,
    provider_facet,select_prior_common_provider_pairs,
    provider_GRIIS_facets,evaluate_state)

A="11111111-1111-4111-8111-111111111111"
B="22222222-2222-4222-8222-222222222222"
C="33333333-3333-4333-8333-333333333333"
D="44444444-4444-4444-8444-444444444444"

def dataset_response():
    return {"count":600,"facets":[
        {"field":"DATASET_KEY","counts":[
            {"name":A,"count":300},{"name":B,"count":200},{"name":C,"count":100}]},
        {"field":"OCCURRENCE_STATUS","counts":[{"name":"DETECTED","count":600}]}
    ]}

def source_response():
    return {"count":300,"facets":[
        {"field":"SPECIES_KEY","counts":[{"name":str(k),"count":30} for k in range(1001,1011)]},
        {"field":"OCCURRENCE_STATUS","counts":[{"name":"DETECTED","count":300}]}
    ]}

class MatchedGBIFSourceTests(unittest.TestCase):
    def test_original_eight_FIPS_and_both_disjoint_2022_source_bins(self):
        self.assertEqual(PILOT_FIPS,("06","12","17","19","36","37","48","53"))
        self.assertEqual(PERIODS,{"before_2022":"2018,2021","after_2022":"2023,2025"})
        self.assertEqual(N_PROVIDER,3)
    def test_same_provider_source_params_not_false_all_provider_facet(self):
        p=provider_taxon_params("California","after_2022",A,0)
        self.assertEqual(p["datasetKey"],A)
        self.assertEqual(p["year"],"2023,2025")
        self.assertEqual(p["facet"],["speciesKey","occurrenceStatus"])
        q=common_provider_params("California","before_2022")
        self.assertNotIn("datasetKey",q)
        self.assertIn("datasetKey",q["facet"])
        with self.assertRaisesRegex(ValueError,"UUID"):
            provider_taxon_params("California","after_2022","fake-id")
    def test_source_top_provider_facet_count_reconciled(self):
        with patch("audit_US_GBIF_matched_provider_temporal_v92.source_get",return_value=dataset_response()):
            result=provider_facet(None,"California","before_2022",600)
        self.assertEqual(result["top100_source_count_coverage"],1.)
        self.assertEqual(len(result["top_dataset_UUID_counts"]),3)
        self.assertTrue(result["top100_not_the_exhaustive_provider_universe"])
    def test_source_provider_selection_only_min_early_late_allPlantae_count(self):
        a={"top_dataset_UUID_counts":{A:1000,B:90,C:40,D:800}}
        b={"top_dataset_UUID_counts":{A:20,B:900,C:50,D:100}}
        providers,m=select_prior_common_provider_pairs(a,b)
        self.assertEqual(providers,[D,B])
        self.assertTrue(m["provider_selected_independent_of_GRIIS_taxa_or_FAF_rail_effect"])
    def test_true_dataset_scoped_speciesKey_facet_does_not_impute_absences(self):
        with patch("audit_US_GBIF_matched_provider_temporal_v92.source_get",return_value=source_response()):
            r=provider_GRIIS_facets(None,"California","before_2022",A,300,
                                    {str(x) for x in range(1001,1012)})
        self.assertEqual(r["n_US_country_GRIIS_indexed_records"],300)
        self.assertEqual(r["n_US_country_GRIIS_recorded_species_key_categories"],10)
        self.assertTrue(r["same_dataset_UUID_does_not_identify_same_observers_or_plots"])
    def test_source_late_only_categories_matched_provider_not_biological_immigration(self):
        original={
            "source_windows":{
                "before_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{
                    "1001":5,"1002":5,"1005":5},
                    "n_all_Plantae_indexed_records":100},
                "after_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{
                    "1001":5,"1003":5,"1004":5},
                    "n_all_Plantae_indexed_records":200}
            }}
        selected={
            A:{
              "before_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{"1001":5},
                             "n_US_country_GRIIS_indexed_records":30,
                             "original_indexed_all_Plantae_records":20},
              "after_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{"1001":5,"1003":5},
                            "n_US_country_GRIIS_indexed_records":40,
                            "original_indexed_all_Plantae_records":60}
            }
        }
        s=evaluate_state({}, {},original,selected,[A],{})
        self.assertEqual(s["original_late_only_source_taxon_categories_NOT_colonizations"],2)
        self.assertEqual(s["n_original_allsource_late_only_taxa_also_visible_from_retained_providers"],1)
        self.assertEqual(s["matched_provider_only_late_taxon_categories"],1)
        self.assertEqual(s["fraction_original_late_only_supported_by_common_UUID_providers"],.5)
        self.assertTrue(s["recorded_late_only_taxa_NOT_new_biological_introductions"])
    def test_earlier_species_from_different_dataset_not_called_new_for_entire_GBIF(self):
        orig={
         "source_windows":{
            "before_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{"1001":5,"1002":5},
                           "n_all_Plantae_indexed_records":100},
            "after_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{"1001":5,"1002":5,"1003":5},
                          "n_all_Plantae_indexed_records":100}}}
        provider={
            A:{
              "before_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{"1001":5},
                             "n_US_country_GRIIS_indexed_records":30,
                             "original_indexed_all_Plantae_records":40},
              "after_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{"1001":5,"1002":5,"1003":5},
                            "n_US_country_GRIIS_indexed_records":60,
                            "original_indexed_all_Plantae_records":40}
            }}
        value=evaluate_state({}, {},orig,provider,[A],{})
        self.assertEqual(value["n_matched_provider_late_only_already_recorded_early_by_OTHER_datasets"],1)
        self.assertEqual(value["n_original_allsource_late_only_taxa_also_visible_from_retained_providers"],1)
    def test_dataset_scope_count_fails_closed_when_drifting(self):
        with patch("audit_US_GBIF_matched_provider_temporal_v92.source_get",return_value=dataset_response()):
            with self.assertRaisesRegex(ValueError,">2.5%"):
                provider_facet(None,"California","after_2022",1000)
if __name__=="__main__":unittest.main()
