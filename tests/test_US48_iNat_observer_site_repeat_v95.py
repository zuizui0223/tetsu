"""Synthetic GBIF original-occurrence iNaturalist source QA tests.

Observer labels (recordedByID/recordedBy) are processed in memory only:
no names, user IDs, geocoordinates or occurrence GBIF IDs in artifacts.
Source pages are NONRANDOM samples, not repeated vegetation surveys.
"""
from __future__ import annotations
import json,sys,unittest
from pathlib import Path
from unittest.mock import patch
from pyproj import Transformer
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_US48_iNat_observer_site_repeat_v95 import (
    PAGES,PAGE_SIZE,CELL_KM,EXPECTED_BINS,
    normalize_observer,source_query,fetch_original_rows,
    cell_of_record,summarize_source_window,temporal_overlap,
    state_record_summary,run,INAT_GBIF_DATASET_UUID)


def record(n,year,observer,lon,lat,species="1001"):
    return {"gbifID":str(n),"year":year,"recordedByID":[observer] if observer else [],
            "recordedBy":"fallbackuser" if observer is None else observer,
            "decimalLongitude":lon,"decimalLatitude":lat,
            "speciesKey":species,"datasetKey":INAT_GBIF_DATASET_UUID}


class HistoricalINatObserverSiteSample(unittest.TestCase):
    def test_frozen_two_periods_source_pages_and_large_source_state_frame(self):
        self.assertEqual(PAGES,(0,1000))
        self.assertEqual(PAGE_SIZE,100)
        self.assertEqual(CELL_KM,(10,50))
        self.assertEqual(EXPECTED_BINS,{"before_2022":"2018,2021",
                                        "after_2022":"2023,2025"})
        self.assertEqual(source_query("Iowa","before_2022",0)["datasetKey"],
                         INAT_GBIF_DATASET_UUID)
        self.assertEqual(source_query("Iowa","after_2022",1000)["year"],
                         "2023,2025")
        with self.assertRaisesRegex(ValueError,"frozen"):
            source_query("Iowa","2022",0)
    def test_observer_prefer_recordedByID_and_no_generic_identity(self):
        self.assertEqual(normalize_observer({"recordedByID":[" USER_17 "],
                                             "recordedBy":"Someone Else"}),"user_17")
        self.assertEqual(normalize_observer({"recordedByID":[],
                                             "recordedBy":" JoHN  DOE "}),"john doe")
        self.assertIsNone(normalize_observer({"recordedBy":"unknown"}))
        self.assertIsNone(normalize_observer({}))
    def test_coarse_10km_and_50km_cells_shift_only_on_original_coords(self):
        t=Transformer.from_crs("EPSG:4326","EPSG:5070",always_xy=True)
        a=record(1,2019,"observerA",-93.617,41.599)
        same=record(2,2023,"observerA",-93.618,41.600)
        self.assertEqual(cell_of_record(a,t,10),cell_of_record(same,t,10))
        self.assertEqual(cell_of_record(a,t,50),cell_of_record(same,t,50))
        self.assertIsNone(cell_of_record(record(3,2023,"X",None,None),t,10))
    def test_same_observer_and_source_grid_recurs_without_leaking_person(self):
        early=[record(1,2019,"sensitive_observer_a",-93.617,41.599),
               record(2,2019,"early_only",-100.0,40.0),
               record(3,2019,"anonymous",-99.0,38.0)]
        late=[record(4,2023,"sensitive_observer_a",-93.618,41.60),
              record(5,2023,"different",-90.0,40.0),
              record(6,2024,"observer_z",-93.618,41.600)]
        a=summarize_source_window(early,{"1001"})
        b=summarize_source_window(late,{"1001"})
        r=temporal_overlap(a,b)
        saved=state_record_summary(a,b,r)
        self.assertEqual(r["n_observer_labels_in_BOTH_bounded_samples"],1)
        self.assertEqual(r["n_late_source_rows_from_observer_seen_in_early_sample"],1)
        self.assertGreaterEqual(r["coarse_site_overlap"]["10km"][
                             "n_geographical_cells_repeated_in_both_sample_year_bins"],1)
        self.assertEqual(r["coarse_site_overlap"]["10km"][
                         "n_late_records_with_same_observer_same_coarse_site_as_early"],1)
        text=json.dumps(saved)
        self.assertNotIn("sensitive_observer_a",text)
        self.assertNotIn("recordedBy",text)
        self.assertNotIn("-93.618",text)
        self.assertNotIn('"gbifID"',text)
    def test_original_GBIF_search_offset_pages_count_and_scope_conservation(self):
        calls=[]
        def fake(_session,query):
            offset=query["offset"];calls.append(offset)
            return {"count":30000,"results":[
                record(offset+i,2019,"A",-93.6,41.5) for i in range(1,101)]}
        with patch("audit_US48_iNat_observer_site_repeat_v95.source_get",
                   side_effect=fake):
            rows,quality=fetch_original_rows(None,"Iowa","before_2022",30000)
        self.assertEqual(calls,[0,1000])
        self.assertEqual(len(rows),200)
        self.assertTrue(quality["nonrandom_ordered_search_API_probe_not_randomized_sampling"])
        self.assertEqual(quality["n_original_source_sampled_records"],200)
        with patch("audit_US48_iNat_observer_site_repeat_v95.source_get",
                   return_value={"count":10000,"results":[
                       record(i,2019,"A",-93,41) for i in range(1,101)]}):
            with self.assertRaisesRegex(ValueError,"drift"):
                fetch_original_rows(None,"Iowa","before_2022",30000)
    def test_full_48_state_source_outputs_aggregates_only(self):
        census={f"{i:02d}":{"name":f"State {i}"}
                for i in range(1,52)}
        inat={
            "complete":True,"source_original_GBIF_qualified_states":48,
            "iNaturalist_original_GBIF_dataset_UUID":INAT_GBIF_DATASET_UUID,
            "states":{f"{i:02d}":{
                "complete":True,"name":f"State {i}",
                "source_result":{"two_GBIF_year_source_periods":{
                    win:{"iNaturalist_original_all_Plantae_record_source_count":10000}
                    for win in EXPECTED_BINS}}}
                for i in range(1,49)}}
        def fake(_session,name,win,denominator):
            obs="stable_mystery_user"
            return [record(1 if win=="before_2022" else 2,
                           2019 if win=="before_2022" else 2023,
                           obs,-93.6,41.59)],{
                "source_indexed_total_records_from_v93":denominator,
                "n_original_source_sampled_records":1,
                "source_GBIF_offset_pages":[{"requested_original_GBIF_offset":0,
                    "n_actual_original_rows":1}],
                "nonrandom_ordered_search_API_probe_not_randomized_sampling":True}
        with patch("audit_US48_iNat_observer_site_repeat_v95.fetch_original_rows",
                   side_effect=fake):
            output=run({"complete":True,"source_periods":dict(EXPECTED_BINS)},inat,set(str(n) for n in range(1000,1200)),
                       census,None)
        self.assertTrue(output["complete"])
        self.assertEqual(output["summary"]["n_Census_state_sources_completed"],48)
        self.assertEqual(output["summary"][
            "n_states_with_at_least_one_same_observer_label_in_both_bounded_samples"],48)
        self.assertNotIn("stable_mystery_user",json.dumps(output))
        self.assertTrue(output["no_raw_observer_names_ids_record_ids_or_latlon_saved"])

if __name__=="__main__":unittest.main()
