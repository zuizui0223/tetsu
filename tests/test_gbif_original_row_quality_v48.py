"""Strict synthetic GBIF original occurrence first-page metadata quality QA.

No real observer names, coordinates or GBIF keys are retained by the output.
"""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
from shapely import wkt

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_gbif_original_row_quality_v48 import (
    frozen_cells,source_params,audit_source_page,run,PAGE_LIMIT
)
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX,PERIODS
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS


def row(cell,period,key=1234,override=None):
    p=wkt.loads(cell["gbif_geometry_wgs84"]).centroid
    yr=int(period.split("_")[0])
    data={
        "key":key,"datasetKey":"dataset-source-uuid",
        "basisOfRecord":"HUMAN_OBSERVATION",
        "recordedBy":"Individual recorder name MUST NOT be exported",
        "occurrenceID":"occ"+str(key),
        "eventID":"event"+str(key),
        "eventDate":f"{yr}-06-15",
        "year":yr,
        "rank":"SPECIES","speciesKey":8012,
        "decimalLongitude":p.x,
        "decimalLatitude":p.y,
        "coordinateUncertaintyInMeters":30,
        "issues":[]
    }
    if override:data.update(override)
    return data


def flora_source():
    return {"complete":True,"source_errors":{},"regions":{
        block:{"country":BLOCKS[block][0],
               "cells":{c["cell_id"]:{
                   "source_pre_outcome_both_period_recording_eligible":True,
                   "periods":{
                       "2010_2017":{
                           "indexed_all_plantae_source_records":2,
                           "alien_observed_minimum_5_species_and_20_records":False},
                       "2018_2025":{
                           "indexed_all_plantae_source_records":2,
                           "alien_observed_minimum_5_species_and_20_records":False}
                   }
               } for c in frozen_cells(block)}
        } for block in BLOCKS_SIX
    }}


class OriginalGBIFRowSourceTests(unittest.TestCase):
    def test_two_fixed_cell_centers_per_country_no_outcome_selected(self):
        self.assertEqual(len(BLOCKS_SIX),6)
        for block in BLOCKS_SIX:
            ids=[r["cell_id"] for r in frozen_cells(block)]
            self.assertEqual(ids,[f"{block}_x1_y1",f"{block}_x2_y2"])

    def test_bounded_api_page_query_matches_original_spatial_geometry(self):
        c=frozen_cells("us_midwest")[0]
        q=source_params(c,"2018_2025")
        self.assertEqual(q["country"],"US")
        self.assertEqual(q["limit"],PAGE_LIMIT)
        self.assertEqual(q["offset"],0)
        self.assertEqual(q["year"],"2018,2025")
        self.assertIn("POLYGON((",q["geometry"])
        self.assertEqual(q["hasGeospatialIssue"],"false")
        with self.assertRaises(ValueError):
            source_params(c,"2018_2025",limit=301)

    def test_actual_records_only_aggregate_recorder_missingness_and_uncertainty(self):
        cell=frozen_cells("fr_burgundy")[0]
        a=row(cell,"2010_2017",1234)
        b=row(cell,"2010_2017",1235,override={
            "recordedBy":None,"coordinateUncertaintyInMeters":None,
            "eventDate":a["eventDate"],"speciesKey":a["speciesKey"],
            "decimalLatitude":a["decimalLatitude"],
            "decimalLongitude":a["decimalLongitude"],
            "occurrenceID":None
        })
        item=audit_source_page(cell,"2010_2017",{
            "count":2,"results":[a,b],"endOfRecords":True
        })
        self.assertTrue(item["records_exhaustively_scanned"])
        self.assertEqual(item["n_missing_recordedBy_first_page"],1)
        self.assertEqual(item["n_missing_coordinateUncertaintyInMeters_first_page"],1)
        self.assertEqual(item["n_same_rounded_point_species_eventdate_collisions_first_page"],1)
        self.assertEqual(item["n_duplicate_dataset_occurrenceID_first_page"],0)
        self.assertEqual(item["n_coords_outside_original_query_polygon_first_page"],0)
        self.assertTrue(item["raw_row_ID_coordinates_and_human_recorder_names_NOT_saved"])
        self.assertNotIn(a["recordedBy"],str(item))

    def test_bounded_first_page_does_not_invent_full_dataset(self):
        cell=frozen_cells("ca_ontario")[1]
        rows=[row(cell,"2018_2025",1000+i) for i in range(PAGE_LIMIT)]
        payload={"count":1000,"results":rows,"endOfRecords":False}
        x=audit_source_page(cell,"2018_2025",payload)
        self.assertFalse(x["records_exhaustively_scanned"])
        self.assertTrue(x["first_page_is_nonrandom_if_more_records"])
        self.assertEqual(x["actual_original_record_rows_examined"],300)
        self.assertEqual(x["indexed_source_search_total"],1000)
        self.assertEqual(x["fraction_of_indexed_records_returned"],.3)

    def test_invalid_nonterminal_missing_first_page_never_claims_full(self):
        c=frozen_cells("gb_yorkshire")[0]
        with self.assertRaisesRegex(ValueError,"Unexpected incomplete"):
            audit_source_page(c,"2010_2017",{"count":300,"results":[row(c,"2010_2017")],
                                         "endOfRecords":False})
        with self.assertRaises(ValueError):
            audit_source_page(c,"2010_2017",{"count":False,"results":[],"endOfRecords":True})

    def test_human_source_row_coordinates_are_checked_not_saved(self):
        c=frozen_cells("de_saxony")[0]
        bad=row(c,"2018_2025",2345,override={
            "decimalLatitude":0.,"decimalLongitude":0.,
            "coordinateUncertaintyInMeters":100000.
        })
        x=audit_source_page(c,"2018_2025",{"count":1,"results":[bad],"endOfRecords":True})
        self.assertEqual(x["n_coords_outside_original_query_polygon_first_page"],1)
        self.assertEqual(x["n_coordinate_uncertainty_gt50000m_first_page"],1)

    def test_six_country_24_nonrandom_source_slots_can_be_queried_without_retaining_raw(self):
        flora=flora_source()
        called=[]
        def fake(_session,params):
            called.append(params)
            from shapely import wkt
            center=wkt.loads(params["geometry"]).centroid
            yr=int(params["year"].split(",")[0])
            return {"count":1,"results":[{
                "key":12345,"datasetKey":"source_uuid","basisOfRecord":"PRESERVED_SPECIMEN",
                "recordedBy":"Privacy-sensitive human recorder",
                "occurrenceID":"id","eventID":"event",
                "eventDate":f"{yr}-06-15","year":yr,
                "rank":"SPECIES","speciesKey":100,
                "decimalLongitude":center.x,"decimalLatitude":center.y,
                "coordinateUncertaintyInMeters":100.
            }],"endOfRecords":True}
        result=run(flora,query=fake,sleep=lambda _:None)
        self.assertTrue(result["complete"])
        self.assertEqual(len(called),24)
        self.assertEqual(result["summary"]["n_original_record_source_bins"],24)
        self.assertEqual(result["summary"]["n_bins_successful"],24)
        self.assertEqual(result["summary"]["n_original_row_records_inspected"],24)
        self.assertNotIn("Privacy-sensitive human recorder",str(result))
        self.assertTrue(result["not_DOI_citable_bulk_occurrence_download"])


if __name__=="__main__":
    unittest.main()
