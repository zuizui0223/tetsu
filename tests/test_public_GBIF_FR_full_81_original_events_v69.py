"""v6.9 synthetic qa for 81 original GBIF France Parquet source keys and event inclusion.

Tests do NOT query or represent actual GBIF events; actual source is
read only by the v6.9 GitHub Actions job from public GBIF S3.
"""
from __future__ import annotations
import copy,sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_public_GBIF_FR_full_81_original_events_v69 import (
    candidate_keys,compare_previous,OLD_BBOX_ROWS,
    EXTRA_BBOX_ROWS,COMBINED_BBOX_ROWS,COMPLETE_CANDIDATES,BLOCK
)
from test_public_GBIF_original_FR_GRIIS_COL_v66 import manifest as _old_source
from test_public_GBIF_FR_unread_33_source_files_v68 import fixture as _all_files


def metadata():
    old=_old_source()
    row=old["regions"][BLOCK]
    row["n_actual_interpreted_event_source_rows_in_selected_files_and_bbox"]=OLD_BBOX_ROWS
    for i,v in enumerate(row["source_object_results"]):
        v["original_interpreted_source_rows_in_200km_bounding_envelope"]=OLD_BBOX_ROWS if i==0 else 0
    _pre,_geo,allfiles=_all_files()
    assert {v["source_key"] for v in row["source_object_results"]}=={v["key"] for v in allfiles[:48]}
    geo={
        "original_region":BLOCK,
        "n_total_source_parquet_files_in_all_true_A5_polygon_intersecting_partitions":81,
        "n_additional_source_files_not_covered_by_old_v65_query":33,
        "source_50km_polygon_plus_centroid_A5_grid_audit":{
            "polygon_cover_disagreed_with_centroid_lookup":True,
        }
    }
    rest=[{
        "source_Parquet_key":f["key"],
        "source_compressed_bytes":f["size_bytes"],
        "n_original_Plantae_rows_in_2018to2025_FR_200km_bbox":EXTRA_BBOX_ROWS if i==0 else 0
    } for i,f in enumerate(allfiles[48:])]
    new={
        "block":BLOCK,"source_country":"FR","complete":True,
        "source_query_errors":[],"n_source_files_unread_or_invalid":0,
        "n_new_additional_source_Parquet_objects_successful":33,
        "2018to2025_raw_all_Plantae_original_rows_from_new_33_files":EXTRA_BBOX_ROWS,
        "2018to2025_combined_81_file_raw_all_Plantae_source_rows_in_same_bbox":COMBINED_BBOX_ROWS,
        "source_additional_object_count_diagnostics":rest
    }
    return old,geo,new


def event_source(new_extra=0,precise_extra=0):
    cells={f"{BLOCK}_x{i}_y{j}":{
        "COL_strict_country_GRIIS_speciesKey_event_counts":{},
        "COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts":{}
    } for i in range(4) for j in range(4)}
    first=f"{BLOCK}_x0_y0"
    cells[first]["COL_strict_country_GRIIS_speciesKey_event_counts"]={"COL_1":22435+new_extra}
    cells[first]["COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts"]={"COL_1":17528+precise_extra}
    return {
        "n_strict_country_GRIIS_distinct_event_species_counts":22435+new_extra,
        "n_strict_country_GRIIS_event_species_known_uncertainty_le5000m":17528+precise_extra,
        "cells":cells
    }


def old_event():
    return {"geographic_block":BLOCK,
        "n_pre_outcome_verified_source_Parquet_objects":48,
        "source_only_period_analyzed":"2018_2025",
        "periods":{"2018_2025":event_source()}
    }


class GBIFAll81FRSourceTests(unittest.TestCase):
    def test_81_keys_are_source_verified_from_original_48_plus_unread_33(self):
        a,b,c=metadata()
        keys,status=candidate_keys(a,b,c)
        self.assertEqual(len(keys),COMPLETE_CANDIDATES)
        self.assertEqual(len(set(keys)),81)
        self.assertEqual(status["n_all_81_bbox_rows_EXPECTED"],697577)
        self.assertTrue(status["candidate_A5_polygon_vs_centroid_disagreement"])
        self.assertTrue(status["NOT_A_FULL_GBIF_RECORD_POPULATION_CENSUS"])
        self.assertEqual(len(status["all_81_source_keys_sha256"]),64)

    def test_missing_source_file_or_681k_rows_disagreement_not_called_absence(self):
        a,b,c=metadata()
        c["source_additional_object_count_diagnostics"].pop()
        with self.assertRaisesRegex(ValueError,"Missing complete"):
            candidate_keys(a,b,c)
        a,b,c=metadata()
        c["2018to2025_combined_81_file_raw_all_Plantae_source_rows_in_same_bbox"]-=1
        with self.assertRaisesRegex(ValueError,"denominator drift"):
            candidate_keys(a,b,c)

    def test_double_counted_source_uri_or_untrusted_S3_cannot_enter_regional_plant_events(self):
        a,b,c=metadata()
        c["source_additional_object_count_diagnostics"][0]["source_Parquet_key"]=(
            a["regions"][BLOCK]["source_object_results"][0]["source_key"])
        with self.assertRaisesRegex(ValueError,"Duplicated GBIF source file"):
            candidate_keys(a,b,c)
        a,b,c=metadata()
        c["source_additional_object_count_diagnostics"][0]["source_Parquet_key"]="' OR 1=1 --"
        with self.assertRaisesRegex(ValueError,"Unapproved additional"):
            candidate_keys(a,b,c)

    def test_new_event_species_counts_and_precise_event_counts_must_contain_old_source(self):
        old=old_event()
        new={"2018_2025":event_source(new_extra=12,precise_extra=9)}
        q=compare_previous(old,new)
        self.assertEqual(q["net_new_country_GRIIS_strict_COL_unique_event_species"],12)
        self.assertEqual(q["net_new_known5km_GRIIS_events"],9)
        self.assertTrue(q["not_proof_of_new_biological_plant_colonizations"])
        self.assertEqual(q["n_original_species_by_cell_categories_with_increased_event_count"],1)
        self.assertEqual(q["n_new_strict_COL_species_by_cell_occurrence_categories"],0)

    def test_source_count_decrease_and_missing_earlier_taxon_are_not_zero_filled(self):
        old=old_event()
        with self.assertRaisesRegex(ValueError,"cannot lose"):
            compare_previous(old,{"2018_2025":event_source(new_extra=-1,precise_extra=0)})
        new=event_source(2,1)
        new["cells"][f"{BLOCK}_x0_y0"]["COL_strict_country_GRIIS_speciesKey_event_counts"]={
            "COL_2":22437
        }
        with self.assertRaisesRegex(ValueError,"decreased a previously observed species"):
            compare_previous(old,{"2018_2025":new})

    def test_never_reclassify_nonmapped_taxa_into_native_absence(self):
        a,b,c=metadata()
        c["complete"]=False
        with self.assertRaisesRegex(ValueError,"missing genuine"):
            candidate_keys(a,b,c)
        old=old_event()
        old["source_only_period_analyzed"]="2010_2017"
        with self.assertRaisesRegex(ValueError,"source provenance missing"):
            compare_previous(old,{"2018_2025":event_source()})


if __name__=="__main__":
    unittest.main()
