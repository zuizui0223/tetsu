"""France Burgundy public GBIF original event source crosswalk/precision QA.

All tests are SYNTHETIC. Real Parquet files are queried in v6.6
GitHub Actions only after original bounded v6.5 source provenance.
"""
from __future__ import annotations

import sys
from pathlib import Path
import unittest
from pyproj import Transformer
from shapely import wkt

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_public_GBIF_original_FR_GRIIS_COL_v66 import (
    audit_rows,source_manifest,FIELDS,BLOCK,PREVIEW_ROOT,
    MAX_SOURCE_PARQUET_OBJECTS
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid


def manifest():
    keys=[PREVIEW_ROOT+"taxonpartition=Plantae/a5_r2=123/data_"+str(i)+".parquet"
          for i in range(MAX_SOURCE_PARQUET_OBJECTS)]
    return {"complete":True,"regions":{
        BLOCK:{
            "genuine_GBIF_original_interpreted_rows_verified":True,
            "source_candidate_objects_in_centroid_A5_bound":71,
            "n_source_objects_actually_scanned":48,
            "source_errors":[],
            "n_actual_interpreted_event_source_rows_in_selected_files_and_bbox":8,
            "source_object_results":[{
                "source_key":k,
                "original_interpreted_source_rows_in_200km_bounding_envelope":8 if i==0 else 0,
            } for i,k in enumerate(keys)],
            "n_source_objects_with_nonzero_actual_source_occurrences":1,
        }
    }}


def synthetic_rows():
    c=next(x for x in make_grid(BLOCK) if x["cell_id"]==f"{BLOCK}_x1_y1")
    p=wkt.loads(c["gbif_geometry_wgs84"]).centroid
    columns=list(FIELDS)
    def row(id,sp="COL-A",occ=None,evt="E1",status="PRESENT",unc=10000,
            year=2018,rank="SPECIES",dataset="D1"):
        record={
            "gbifid":str(id),"datasetkey":dataset,"countrycode":"FR",
            "year":year,"specieskey":sp,"occurrencestatus":status,
            "decimallongitude":p.x,"decimallatitude":p.y,
            "coordinateuncertaintyinmeters":unc,
            "occurrenceid":occ or "O"+str(id),
            "eventid":evt,"eventdate":str(year)+"-07-14",
            "taxonrank":rank
        }
        return tuple(record.get(k) for k in columns)
    rows=[
        row(1,occ="O1",unc=10000),
        row(2,occ="O2",unc=1000),
        row(3,sp="COL-B",occ="O3",evt="E2",unc=4000,dataset="D2"),
        row(4,status="ABSENT",occ="O4",evt="E3"),
        row(5,sp="COL-UNKNOWN",occ="O5",evt="E4"),
        row(6,occ="O6",evt="E5",year=2013,unc=200),
        row(1,occ="O1",evt="E1",unc=10000),
        row(8,occ="O1",evt="E1",unc=8000),
        row(9,occ="O9",evt="E10",rank="GENUS"),
    ]
    return rows,columns


class AuthenticFRPublicGBIFEventSourceTests(unittest.TestCase):
    def test_source_manifest_requires_exact_preoutcome_48_public_objects(self):
        keys,src=source_manifest(manifest())
        self.assertEqual(len(keys),48)
        self.assertEqual(src["source_candidate_objects_in_centroid_A5_bound"],71)
        bad=manifest()
        bad["regions"][BLOCK]["source_object_results"][1]["source_key"]="s3://private/secret.parquet"
        with self.assertRaisesRegex(ValueError,"Unapproved external"):
            source_manifest(bad)
        bad=manifest()
        bad["regions"][BLOCK]["source_object_results"].pop()
        with self.assertRaises(ValueError):
            source_manifest(bad)

    def test_event_group_dedup_is_order_independent_for_coordinate_precision(self):
        rows,columns=synthetic_rows()
        late=sum(row[columns.index("year")]==2018 for row in rows)
        a,counts=audit_rows(rows,columns,{"COL-A"},{"COL-A","COL-B"},
                            expected_late_bbox_count=late)
        b,counts_reverse=audit_rows(list(reversed(rows)),columns,
                            {"COL-A"},{"COL-A","COL-B"},
                            expected_late_bbox_count=late)
        cell=f"{BLOCK}_x1_y1"
        late_a=a["2018_2025"]["cells"][cell]
        late_b=b["2018_2025"]["cells"][cell]
        self.assertEqual(late_a["n_GRIIS_strict_COL_event_species"],1)
        self.assertEqual(late_a["n_GRIIS_strict_COL_event_species_known_uncertainty_le5000m"],1)
        self.assertEqual(late_b["n_GRIIS_strict_COL_event_species_known_uncertainty_le5000m"],1)
        self.assertEqual(late_a["n_GRIIS_strict_plus_synonym_COL_event_species"],2)
        self.assertEqual(a["2010_2017"]["n_strict_country_GRIIS_distinct_event_species_counts"],1)
        self.assertEqual(counts["n_source_duplicate_same_event_species_group"],2)
        self.assertEqual(counts_reverse["n_source_duplicate_same_event_species_group"],2)
        self.assertEqual(counts["n_source_duplicate_GBIF_unique_ID"],1)
        self.assertEqual(counts["n_source_duplicate_same_dataset_occurrenceID"],1)
        self.assertEqual(counts["n_source_non_PRESENT_or_unknown_status"],1)
        self.assertEqual(counts["n_source_speciesKey_not_in_country_GRIIS_accepted_unknown"],1)

    def test_conflicting_same_dataset_occurrenceID_event_is_not_chosen_by_input_order(self):
        rows,columns=synthetic_rows()
        conflict=list(rows[-2])
        conflict[columns.index("eventid")]="CONFLICT_OTHER_EVENT"
        altered=list(rows)
        altered[-2]=tuple(conflict)
        with self.assertRaisesRegex(ValueError,"Conflicting event/taxon/cell"):
            audit_rows(altered,columns,{"COL-A"},{"COL-A","COL-B"})

    def test_missing_source_status_does_not_mean_plant_absence_or_zero(self):
        rows,columns=synthetic_rows()
        bad=[list(x) for x in rows]
        for x in bad:
            x[columns.index("occurrencestatus")]="NOT_DETECTED"
        byperiod,stats=audit_rows([tuple(x) for x in bad],columns,
                                   {"COL-A"},{"COL-A","COL-B"})
        self.assertEqual(byperiod["2018_2025"]["n_strict_country_GRIIS_distinct_event_species_counts"],0)
        self.assertGreater(stats["n_source_non_PRESENT_or_unknown_status"],0)

    def test_same_source_archive_period_denominator_mismatch_fails_closed(self):
        rows,columns=synthetic_rows()
        with self.assertRaisesRegex(ValueError,"disagrees with v6.5"):
            audit_rows(rows,columns,{"COL-A"},{"COL-A","COL-B"},
                       expected_late_bbox_count=671528)

    def test_incompatible_source_COL_namespace_cannot_be_treated_as_native(self):
        rows,columns=synthetic_rows()
        with self.assertRaises(ValueError):
            audit_rows(rows,columns,set(),{"COL-A"})
        with self.assertRaises(ValueError):
            audit_rows(rows,columns,{"legacy_old_key"},{"COL-A"})
        with self.assertRaisesRegex(ValueError,"budget exceeded"):
            audit_rows(rows,columns,{"COL-A"},{"COL-A","COL-B"},
                       max_rows=2)


if __name__=="__main__":
    unittest.main()
