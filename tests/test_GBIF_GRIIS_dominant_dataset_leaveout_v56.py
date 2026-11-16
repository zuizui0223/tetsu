"""Synthetic regression of full GBIF speciesKey taxa vs single dominant dataset source subtraction."""
import sys
from pathlib import Path
import unittest
import copy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_gbif_griis_dominant_dataset_leaveout_v56 import (
    query_species,whole_source_species,split_top_and_other,run_one,
    MIN_COUNTRY_SOURCE_ALIEN_RECORDS
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid,BLOCKS
from analyze_6country_50km_network_plants_v44 import PERIODS,BLOCKS_SIX
from analyze_macro_rail_assemblages_v22 import COUNTRIES
from test_gbif_dataset_composition_v50 import fake_flora


def checklists():
    countries=sorted(set(COUNTRIES.values())-{"TW"})
    return {"complete":True,"errors":{},"countries":{
        country:{
            "country":country,
            "full_source_pages_verified":True,
            "gbif_backbone_plant_species_key_set":[str(100+i) for i in range(110)],
            "n_directly_matched_plant_backbone_species_keys":110
        } for country in countries
    }}


def species_response(params):
    dataset=params.get("datasetKey")
    if dataset is None:
        source={"100":60,"101":20,"102":15,"103":10,"104":10,"999":5}
        n=120
    elif dataset=="dataset_a":
        source={"100":45,"101":15,"102":5,"999":5}
        n=70
    else:
        raise AssertionError("Not source-selected top-one dataset")
    return {
        "count":n,
        "facets":[
            {"field":"SPECIES_KEY","counts":[
                {"name":key,"count":v} for key,v in source.items()
            ]},
            {"field":"OCCURRENCE_STATUS","counts":[
                {"name":"DETECTED","count":n-1},
                {"name":"NOT_DETECTED","count":1}
            ]}
        ]
    }


def dataset_source(block="us_midwest"):
    f=fake_flora(block)
    entries={}
    for cid,c in f["regions"][block]["cells"].items():
        include=c["source_pre_outcome_both_period_recording_eligible"]
        entries[cid]={
            "country":BLOCKS[block][0],
            "original_GRIIS_two_period_source_eligible":include,
            "periods":{p:{
                "source_profile_valid_for_recording_bias_covariates":True,
                "top_datasetKey_source_id":"dataset_a",
                "source_indexed_all_Plantae_records":120,
                "dataset_key_record_counts":{"dataset_a":70,"dataset_b":30,"dataset_c":20}
            } for p in PERIODS} if include else {}
        }
    return {"complete":True,"errors":{},"block":block,"country":BLOCKS[block][0],
            "cells":entries}


class ActualDatasetGBIFSpeciesFacetSourceTests(unittest.TestCase):
    def test_source_taxon_facet_without_incompatible_present(self):
        cell=make_grid("fr_burgundy")[0]
        p=query_species(cell,"2018_2025","top_dset",0)
        self.assertNotIn("occurrenceStatus",p)
        self.assertEqual(p["datasetKey"],"top_dset")
        self.assertEqual(p["facet"],["speciesKey","occurrenceStatus"])
        self.assertEqual(p["year"],"2018,2025")
        self.assertEqual(p["limit"],0)
        with self.assertRaises(ValueError):
            query_species(cell,"2018_2025","top_dset",51)

    def test_full_vs_top_dataset_partition_restores_original_record_counts(self):
        cell=make_grid("us_midwest")[0]
        allsrc=whole_source_species(None,cell,"2018_2025",request=lambda c,p:species_response(p),sleep=lambda _:None)
        top=whole_source_species(None,cell,"2018_2025","dataset_a",request=lambda c,p:species_response(p),sleep=lambda _:None)
        old={"indexed_all_plantae_source_records":120,
             "n_observed_GRIIS_introduced_plant_species_keys":5}
        split=split_top_and_other(allsrc,top,set(str(100+i) for i in range(110)),old)
        self.assertEqual(split["n_alien_species_keys_all_indexed_sources"],5)
        self.assertEqual(split["n_species_keys_seen_only_in_top_dataset"],0)
        self.assertEqual(split["n_species_keys_seen_only_in_other_datasets"],2)
        self.assertEqual(split["n_species_keys_seen_in_top_AND_other_datasets"],3)
        self.assertEqual(split["n_indexed_country_GRIIS_alien_records_total"],115)
        self.assertEqual(split["n_indexed_country_GRIIS_alien_records_top_dataset"],65)
        self.assertEqual(split["n_indexed_country_GRIIS_alien_records_other_datasets"],50)
        self.assertTrue(split["source_alien_min5_species_and20records_without_top_dataset"])
        self.assertTrue(split["dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys"])

    def test_bad_filtered_species_exceeds_total_causes_strict_source_failure(self):
        cell=make_grid("us_midwest")[0]
        a=whole_source_species(None,cell,"2018_2025",request=lambda c,p:species_response(p),sleep=lambda _:None)
        b=whole_source_species(None,cell,"2018_2025","dataset_a",request=lambda c,p:species_response(p),sleep=lambda _:None)
        b["speciesKey_count_by_all_taxa"]["100"]=65
        with self.assertRaisesRegex(ValueError,"exceed SAME-query total"):
            split_top_and_other(a,b,set(str(100+i) for i in range(110)),
                               {"indexed_all_plantae_source_records":120})

    def test_historical_count_drift_does_not_allow_source_imputation(self):
        cell=make_grid("ca_ontario")[0]
        a=whole_source_species(None,cell,"2018_2025",request=lambda c,p:species_response(p),sleep=lambda _:None)
        b=whole_source_species(None,cell,"2018_2025","dataset_a",request=lambda c,p:species_response(p),sleep=lambda _:None)
        with self.assertRaisesRegex(ValueError,"drift"):
            split_top_and_other(a,b,set(str(100+i) for i in range(110)),
                               {"indexed_all_plantae_source_records":500})

    def test_real_source_cell_frame_status_and_dual_period_are_preserved(self):
        block="us_midwest";f=fake_flora(block)
        d=dataset_source(block)
        called=[]
        def req(client,params):
            called.append(params)
            return species_response(params)
        out=run_one(block,f,d,checklists(),request=req,sleep=lambda _:None)
        self.assertTrue(out["complete"])
        self.assertEqual(out["source_original_GRIIS_eligible_cells"],2)
        self.assertEqual(out["n_both_period_full_total_top_and_other_source_cells"],2)
        self.assertEqual(out["n_both_period_without_top_dataset_min5species_20records"],2)
        self.assertEqual(len(called),8)
        self.assertTrue(all("occurrenceStatus" not in q for q in called))
        self.assertEqual(len(out["source_cells"]),16)
        self.assertTrue(out["late_period_median_fraction_source_GRIIS_species_keys_in_top_AND_other_datasets"]>0)
        self.assertNotIn("recordedBy",str(out))

    def test_invalid_country_source_fails_closed(self):
        with self.assertRaises(ValueError):
            run_one("not_frozen",fake_flora(),dataset_source(),checklists(),
                    request=lambda c,p:species_response(p),sleep=lambda _:None)


if __name__=="__main__":
    unittest.main()
