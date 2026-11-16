"""Actual public GBIF A5 source original-event 50km country GRIIS QA, synthetic row fixture."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
from shapely import wkt
from pyproj import Transformer
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_public_GBIF_original_GRIIS_rows_v61 import (
    BLOCK,SOURCE_FILE,audit_rows,source_query
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid
from ingest_GBIF_DWCA_original_events_v58 import source_cell_frame,_cell_of_lonlat

COLS=("gbifid","datasetkey","countrycode","year","specieskey",
      "occurrencestatus","decimallongitude","decimallatitude",
      "coordinateuncertaintyinmeters","occurrenceid","eventid",
      "eventdate","taxonrank")


def rec(p,uid,sp="100",yr=2018,status="PRESENT",unc=1000,
        event="eventA",occ=None,ds="source_dataset"):
    return (uid,ds,"US",yr,sp,status,p.x,p.y,unc,occ or uid,event,
            f"{yr}-05-15","SPECIES")


class FakeDB:
    def execute(self,sql):
        self.sql=sql
        return self
    def fetchall(self):
        return [(x,) for x in COLS]


class PublicGBIFRealRowSourceQATests(unittest.TestCase):
    def test_original_2500km2_source_cell_location(self):
        block=BLOCK
        c=next(x for x in make_grid(block) if x["cell_id"]==f"{block}_x1_y1")
        p=wkt.loads(c["gbif_geometry_wgs84"]).centroid
        frame=source_cell_frame(block)
        project=Transformer.from_crs(4326,"EPSG:6933",always_xy=True)
        self.assertEqual(_cell_of_lonlat(block,p.x,p.y,project,frame),
                         f"{block}_x1_y1")

    def test_same_source_file_schema_and_original_source_bbox_are_guarded(self):
        fake=FakeDB()
        query,fields,frame,meta=source_query(fake)
        self.assertIn("s3://gbif-public-data/development/2026-08-01-DWCA",query)
        self.assertEqual(set(fields),set(COLS))
        self.assertIn('countrycode',query.lower())
        self.assertIn('2010 AND 2025',query)
        self.assertEqual(len(frame["wkt"]),16)
        self.assertTrue(meta["bbox_west"]<meta["bbox_east"])

    def test_row_grain_presence_dedup_and_precision_excludes_uncertain(self):
        block=BLOCK
        c=next(x for x in make_grid(block) if x["cell_id"]==f"{block}_x1_y1")
        p=wkt.loads(c["gbif_geometry_wgs84"]).centroid
        rows=[
            rec(p,"id1","100",unc=1000,event="eventA"),
            rec(p,"id2","100",unc=1000,event="eventA"),
            rec(p,"id3","101",unc=12000,event="eventB"),
            rec(p,"id4","100",unc=None,event="eventC"),
            rec(p,"id5","999",unc=500,event="eventD"),
            rec(p,"id6","100",status="ABSENT",event="eventE"),
            rec(p,"id1","100",event="eventF"),
            rec(p,"id7","100",yr=2014,event="eventG"),
        ]
        summary,counts=audit_rows(rows,COLS,{"100","101"})
        late=summary["2018_2025"][f"{block}_x1_y1"]
        early=summary["2010_2017"][f"{block}_x1_y1"]
        self.assertEqual(late["n_removed_duplicate_dataset_event_species"],1)
        self.assertEqual(late["n_source_GRIIS_keys"],2)
        self.assertEqual(late["country_GRIIS_specieskeys_original_row_counts"],{"100":2,"101":1})
        self.assertEqual(late["country_GRIIS_specieskeys_precise_5000m_counts"],{"100":1})
        self.assertEqual(late["n_missing_coordinate_uncertainty"],1)
        self.assertEqual(late["n_known_uncertainty_gt5000m"],1)
        self.assertEqual(early["country_GRIIS_specieskeys_original_row_counts"],{"100":1})
        self.assertEqual(counts["speciesKey_not_in_country_GRIIS_status_UNKNOWN"],1)
        self.assertEqual(counts["non_PRESENCE_or_unknown_status_rows_not_used_as_plant_presence"],1)
        self.assertEqual(counts["missing_or_duplicate_GBIF_unique_source_key"],1)
        self.assertNotIn("id1",str(summary))

    def test_nonfrozen_source_region_not_accepted(self):
        fake=FakeDB()
        with self.assertRaises(ValueError):
            source_query(fake,"nz_north")
        with self.assertRaises(ValueError):
            audit_rows([],COLS,{"100"},"gb_yorkshire")

    def test_source_accepted_taxa_cannot_be_assumed_native_if_absent(self):
        block=BLOCK
        c=next(x for x in make_grid(block) if x["cell_id"]==f"{block}_x0_y0")
        p=wkt.loads(c["gbif_geometry_wgs84"]).centroid
        cells,counts=audit_rows([rec(p,"idA","555")],COLS,{"100"})
        self.assertEqual(counts["speciesKey_not_in_country_GRIIS_status_UNKNOWN"],1)
        self.assertEqual(cells["2018_2025"][f"{block}_x0_y0"]["n_source_GRIIS_keys"],0)


if __name__=="__main__":
    unittest.main()
