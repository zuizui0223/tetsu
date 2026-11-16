"""Offline GBIF 50km datasetKey/source recording bias QA and current-status schema tests."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_gbif_dataset_composition_6country_v50 import (
    source_query,_facet_rows,full_dataset_profile,run_one_block,BLOCKS_SIX
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid,BLOCKS
from analyze_6country_50km_network_plants_v44 import PERIODS


def good_source():
    return {"count":120,"facets":[
        {"field":"DATASET_KEY","counts":[
            {"name":"dataset_a","count":70},
            {"name":"dataset_b","count":30},
            {"name":"dataset_c","count":20}]},
        {"field":"OCCURRENCE_STATUS","counts":[
            {"name":"DETECTED","count":119},
            {"name":"NOT_DETECTED","count":1}]},
    ]}


def fake_flora(block="us_midwest"):
    country=BLOCKS[block][0]
    cells={}
    for i,r in enumerate(make_grid(block)):
        eligible=i in (0,1)
        cells[r["cell_id"]]={
            "cell_id":r["cell_id"],
            "country":country,
            "source_pre_outcome_both_period_recording_eligible":eligible,
            "periods":{
                period:{
                    "alien_observed_minimum_5_species_and_20_records":eligible,
                    "indexed_all_plantae_source_records":120}
                for period in PERIODS
            }
        }
    return {"complete":True,"source_errors":{},
            "regions":{block:{"country":country,"cells":cells}}}


class DatasetConcentrationSourceTests(unittest.TestCase):
    def test_no_broken_legacy_OCCURRENCE_STATUS_present_query(self):
        p=source_query(make_grid("fr_burgundy")[0],"2018_2025",0)
        self.assertNotIn("occurrenceStatus",p)
        self.assertEqual(p["facet"],["datasetKey","occurrenceStatus"])
        self.assertEqual(p["datasetKey.facetLimit"],5000)
        self.assertEqual(p["limit"],0)
        self.assertEqual(p["country"],"FR")
        self.assertEqual(p["year"],"2018,2025")
        with self.assertRaises(ValueError):
            source_query(make_grid("fr_burgundy")[0],"2018_2025",100)

    def test_true_facet_hhi_and_two_orthogonal_marginals(self):
        p=full_dataset_profile(None,make_grid("us_midwest")[0],"2018_2025",120,
              request=lambda session,params:good_source(),sleep=lambda _:None)
        self.assertTrue(p["source_profile_valid_for_recording_bias_covariates"])
        self.assertEqual(p["status_source_counts"],{"DETECTED":119,"NOT_DETECTED":1})
        self.assertAlmostEqual(p["top_datasetKey_source_share"],70/120,places=5)
        self.assertEqual(p["n_distinct_source_dataset_keys"],3)
        self.assertAlmostEqual(p["dataset_source_HHI_including_unattributed_bucket"],
                               (70**2+30**2+20**2)/120**2,places=5)
        self.assertTrue(p["NOT_GRIIS_aliens_only_source_dataset_distribution"])
        self.assertTrue(p["raw_GBIF_occurrence_ID_recordedBy_coordinates_NOT_stored"])

    def test_missing_source_dataset_categories_are_unknown_not_shared(self):
        data=good_source()
        data["facets"][0]["counts"][0]["count"]=50
        p=full_dataset_profile(None,make_grid("us_midwest")[0],"2018_2025",120,
              request=lambda session,params:data,sleep=lambda _:None)
        self.assertEqual(p["missing_datasetKey_indexed_records"],20)
        self.assertFalse(p["source_profile_valid_for_recording_bias_covariates"])

    def test_non_detected_not_silently_counted_as_introduced_plants(self):
        data=good_source()
        data["facets"][1]["counts"][0]["count"]=60
        data["facets"][1]["counts"][1]["count"]=60
        p=full_dataset_profile(None,make_grid("us_midwest")[0],"2018_2025",120,
              request=lambda session,params:data,sleep=lambda _:None)
        self.assertFalse(p["source_profile_valid_for_recording_bias_covariates"])
        self.assertEqual(p["status_NOT_DETECTED_share"],.5)

    def test_duplicate_facet_source_and_count_drift_fail_closed(self):
        data=good_source()
        data["facets"][0]["counts"][1]["name"]="dataset_a"
        with self.assertRaises(ValueError):
            full_dataset_profile(None,make_grid("us_midwest")[0],"2018_2025",120,
                request=lambda session,params:data,sleep=lambda _:None)
        q=full_dataset_profile(None,make_grid("us_midwest")[0],"2018_2025",1000,
            request=lambda session,params:good_source(),sleep=lambda _:None)
        self.assertFalse(q["source_profile_valid_for_recording_bias_covariates"])

    def test_source_data_quality_gate_uses_preexisting_introduced_cells_not_new_outcome(self):
        calls=[]
        def mocked(session,params):
            calls.append(params)
            return good_source()
        data=run_one_block("us_midwest",fake_flora(),
                           request=mocked,sleep=lambda _:None)
        self.assertTrue(data["complete"])
        self.assertEqual(data["n_original_country_GRIIS_source_eligible_cells"],2)
        self.assertEqual(data["n_source_dataset_facet_complete_both_period_cells"],2)
        self.assertEqual(data["n_source_covariate_valid_both_period_cells"],2)
        self.assertEqual(len(calls),4)
        self.assertTrue(all("occurrenceStatus" not in x for x in calls))
        self.assertEqual(len(data["cells"]),16)
        with self.assertRaises(ValueError):
            run_one_block("not_frozen",fake_flora(),request=mocked,sleep=lambda _:None)


if __name__=="__main__":
    unittest.main()
