"""Offline GBIF DWCA original interpreted ZIP QA; synthetic only, no actual DOI claim."""
from __future__ import annotations
from pathlib import Path
import csv,io,json,sys,tempfile,zipfile
import unittest
from shapely import wkt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from ingest_GBIF_DWCA_original_events_v58 import (
    analyze_archive,metadata_contract,source_cell_frame,_cell_of_lonlat
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid
from analyze_macro_rail_assemblages_v22 import COUNTRIES
from pyproj import Transformer

FAKE_META={"key":"0001005-130906152512535","status":"SUCCEEDED",
           "doi":"10.15468/dl.syntheticTEST",
           "downloadLink":"https://api.gbif.org/occurrence/download/request/0001005-130906152512535.zip"}


def fake_griis():
    countries=set(COUNTRIES.values())-{"TW"}
    return {"complete":True,"errors":{},"countries":{
        country:{"country":country,"full_source_pages_verified":True,
          "n_directly_matched_plant_backbone_species_keys":110,
          "gbif_backbone_plant_species_key_set":[str(100+i) for i in range(110)]
        } for country in countries
    }}


def fake_dwca(path,block="us_midwest"):
    c=next(x for x in make_grid(block) if x["cell_id"]==f"{block}_x1_y1")
    center=wkt.loads(c["gbif_geometry_wgs84"]).centroid
    lon=center.x;lat=center.y
    header=[
        "gbifID","datasetKey","occurrenceID","eventID","eventDate","year",
        "countryCode","speciesKey","taxonRank","kingdom","occurrenceStatus",
        "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
        "recordedBy"
    ]
    def item(uid="1",dataset="A",occ="oa",event="e1",yr="2018",sp="100",
             status="PRESENT",unc="1000",country="US",recorder="Sensitive name"):
        return dict(zip(header,[uid,dataset,occ,event,yr+"-07-01",yr,country,
                   sp,"SPECIES","Plantae",status,str(lon),str(lat),unc,recorder]))
    rows=[
        item(),
        item(uid="2",occ="ob",event="e1"),
        item(uid="3",dataset="B",occ="oc",event="e2",sp="101",unc="10000"),
        item(uid="4",dataset="B",occ="od",event="e3",sp="102",status="ABSENT"),
        item(uid="5",dataset="B",occ="oe",event="e4",sp="999"),
        item(),
        item(uid="7",occ="oa",event="e5"),
        item(uid="8",occ="of",event="",unc="",sp="100"),
        item(uid="9",occ="og",event="e6",country="CA"),
        item(uid="10",occ="oh",event="e7",yr="2014",sp="104",unc="200"),
    ]
    buf=io.StringIO(newline="")
    w=csv.DictWriter(buf,fieldnames=header,delimiter="\t",quoting=csv.QUOTE_MINIMAL)
    w.writeheader()
    w.writerows(rows)
    with zipfile.ZipFile(path,"w",zipfile.ZIP_DEFLATED) as z:
        z.writestr("occurrence.txt",buf.getvalue())
        z.writestr("meta.xml",'<archive xmlns="http://rs.tdwg.org/dwc/text/"/>')
        z.writestr("citations.txt","SYNTHETIC ONLY\n")
    return rows


class OriginalDWCAParserTests(unittest.TestCase):
    def test_metadata_identity_requires_real_source_shaped_key_and_doi(self):
        result=metadata_contract(FAKE_META)
        self.assertEqual(result["original_download_status"],"SUCCEEDED")
        self.assertTrue(result["user_supplied_metadata_not_verified_against_live_GBIF_API"])
        with self.assertRaises(ValueError):
            metadata_contract({**FAKE_META,"status":"PREPARING"})
        with self.assertRaises(ValueError):
            metadata_contract({**FAKE_META,"doi":""})
        with self.assertRaises(ValueError):
            metadata_contract({**FAKE_META,"downloadLink":"file:///tmp/fake.zip"})

    def test_frozen_original_50km_cell_identity(self):
        block="us_midwest"
        origin=source_cell_frame(block)
        self.assertEqual(len(origin["wkt"]),16)
        c=next(x for x in make_grid(block) if x["cell_id"]==f"{block}_x1_y1")
        p=wkt.loads(c["gbif_geometry_wgs84"]).centroid
        cell=_cell_of_lonlat(block,p.x,p.y,Transformer.from_crs(4326,"EPSG:6933",always_xy=True),origin)
        self.assertEqual(cell,"us_midwest_x1_y1")

    def test_real_zip_parsing_on_SYNTHETIC_record_fixture_preserves_privacy_and_temporal_QA(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/"synthetic.zip"
            fake_dwca(f)
            result=analyze_archive(f,"us_midwest",fake_griis(),FAKE_META)
            self.assertTrue(result["actual_original_archive_file_read"])
            self.assertFalse(result["source_metadata_matches_live_GBIF_API"])
            self.assertTrue(result["DOI_archived_event_data_actually_read"])
            self.assertTrue(result["no_train_railway_plant_causality_identified"])
            a=result["periods"]["2018_2025"]["us_midwest_x1_y1"]
            b=result["periods"]["2010_2017"]["us_midwest_x1_y1"]
            self.assertEqual(result["n_original_GBIF_DWCA_rows_processed"],10)
            self.assertEqual(a["n_records_rejected_duplicate_GBIF_ids"],1)
            self.assertEqual(a["n_records_rejected_duplicate_dataset_occurrenceID"],1)
            self.assertEqual(a["n_event_species_raw_composites_deduplicated"],1)
            self.assertEqual(a["country_GRIIS_introduced_event_speciesKey_count_all_valid_geometry"],
                {"100":2,"101":1})
            self.assertEqual(a["country_GRIIS_introduced_event_speciesKey_count_precise_5000m"],
                {"100":1})
            self.assertEqual(a["country_GRIIS_introduced_event_speciesKey_count_sensitivity_25000m"],
                {"100":1,"101":1})
            self.assertEqual(b["country_GRIIS_introduced_event_speciesKey_count_precise_5000m"],
                {"104":1})
            self.assertEqual(result["source_quality_counters"]["n_nonpresence_or_unknown_source_status"],1)
            self.assertEqual(result["source_quality_counters"]["n_country_GRIIS_nonmember_unknown_status"],1)
            self.assertEqual(result["source_quality_counters"]["n_nonmatching_source_country"],1)
            txt=json.dumps(result)
            self.assertNotIn("Sensitive name",txt)
            self.assertNotIn("oa",str(a))
            self.assertNotIn("ob",str(a))

    def test_invalid_archive_not_misreported_as_real_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/"corrupt.zip";f.write_bytes(b"not a DWCA")
            with self.assertRaisesRegex(ValueError,"ZIP not found/valid"):
                analyze_archive(f,"us_midwest",fake_griis(),FAKE_META)
            with self.assertRaises(ValueError):
                analyze_archive(f,"outside",fake_griis(),FAKE_META)


if __name__=="__main__":
    unittest.main()
