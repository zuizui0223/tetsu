"""Synthetic strict current GBIF status FACET + first-300 original row QA."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
from shapely import wkt

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_gbif_detected_original_rows_v48b import (
    params_without_broken_legacy_status,parse_status_facet,run
)
from audit_gbif_original_row_quality_v48 import frozen_cells
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX


def floras():
    return {"complete":True,"source_errors":{},"regions":{
        block:{"country":BLOCKS[block][0],"cells":{
            cell["cell_id"]:{
                "source_pre_outcome_both_period_recording_eligible":True,
                "periods":{
                    period:{"indexed_all_plantae_source_records":3,
                            "alien_observed_minimum_5_species_and_20_records":False}
                    for period in ("2010_2017","2018_2025")
                }
            } for cell in frozen_cells(block)
        }} for block in BLOCKS_SIX
    }}


def row(params):
    p=wkt.loads(params["geometry"]).centroid
    year=int(params["year"].split(",")[0])
    return {"key":42,"datasetKey":"source_dataset","basisOfRecord":"HUMAN_OBSERVATION",
        "recordedBy":"Sensitive observer not stored",
        "occurrenceID":"identifier-not-stored",
        "eventDate":str(year)+"-06-15","year":year,
        "speciesKey":100,"rank":"SPECIES",
        "decimalLongitude":p.x,"decimalLatitude":p.y,
        "coordinateUncertaintyInMeters":20.
    }


class GBIFDetectedOriginalRowsTests(unittest.TestCase):
    def test_current_get_legacy_present_filter_is_explicitly_removed(self):
        c=frozen_cells("us_midwest")[0]
        p=params_without_broken_legacy_status(c,"2018_2025")
        self.assertNotIn("occurrenceStatus",p)
        self.assertEqual(p["limit"],300)
        self.assertEqual(p["country"],"US")

    def test_DETECTED_and_NOT_DETECTED_facet_exhaustive_status(self):
        a=parse_status_facet({"count":101,"facets":[{
            "field":"OCCURRENCE_STATUS",
            "counts":[{"name":"DETECTED","count":100},
                      {"name":"NOT_DETECTED","count":1}]
        }]})
        self.assertEqual(a["indexed_source_count"],101)
        self.assertEqual(a["statuses"]["NOT_DETECTED"],1)
        self.assertEqual(a["unclassified_count"],0)
        with self.assertRaises(ValueError):
            parse_status_facet({"count":1,"facets":[{"field":"OCCURRENCE_STATUS",
                   "counts":[{"name":"DETECTED","count":2}]}]})

    def test_missing_status_not_assumed_detected(self):
        p=parse_status_facet({"count":10,"facets":[{
            "field":"OCCURRENCE_STATUS",
            "counts":[{"name":"DETECTED","count":7}]}]})
        self.assertEqual(p["unclassified_count"],3)
        self.assertEqual(p["statuses"]["DETECTED"],7)
        with self.assertRaises(ValueError):
            parse_status_facet({"count":1,"facets":[]})

    def test_run_exact_24_bins_no_raw_names_saved(self):
        calls=[]
        def fake(_session,params):
            calls.append(params)
            if params.get("facet")=="OCCURRENCE_STATUS":
                return {"count":1,"facets":[{"field":"OCCURRENCE_STATUS","counts":[
                    {"name":"DETECTED","count":1}]}]}
            return {"count":1,"results":[row(params)],"endOfRecords":True}
        out=run(floras(),query=fake,sleep=lambda _:None)
        self.assertTrue(out["complete"])
        self.assertEqual(len(calls),48)
        self.assertEqual(out["summary"]["n_facet_plus_original_row_source_bins_complete"],24)
        self.assertEqual(out["summary"]["n_firstpage_actual_source_occurrences_inspected"],24)
        self.assertEqual(out["summary"]["n_facet_DETECTED_indexed_records_full_24_selected_slots"],24)
        self.assertEqual(out["summary"]["n_facet_NOT_DETECTED_indexed_records_full_24_selected_slots"],0)
        self.assertNotIn("Sensitive observer",str(out))
        self.assertTrue(out["no_raw_personal_names_or_precise_geocoordinates_exported"])

    def test_prior_positive_counts_cannot_reconcile_zero_live_detected(self):
        def fake(_session,params):
            if params.get("facet")=="OCCURRENCE_STATUS":
                return {"count":0,"facets":[]}
            return {"count":0,"results":[],"endOfRecords":True}
        out=run(floras(),query=fake,sleep=lambda _:None)
        self.assertFalse(out["complete"])
        self.assertEqual(len(out["errors"]),24)
        self.assertEqual(out["summary"]["n_facet_plus_original_row_source_bins_complete"],0)


if __name__=="__main__":
    unittest.main()
