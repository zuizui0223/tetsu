"""Offline GRIIS-to-GBIF-backbone country checklist mapping checks."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_griis_country_backbone_keys_v24 import (
    DATASETS,source_params,verify_dataset,parse_species_page,audit_country
)


class GRIISBackboneTests(unittest.TestCase):
    def test_country_lists_are_explicit_and_taiwan_not_falsely_native(self):
        self.assertEqual(len(DATASETS),9)
        self.assertNotIn("TW",DATASETS)
        self.assertEqual(len({x[0] for x in DATASETS.values()}),9)

    def test_source_search_is_checklist_only_not_global_plants(self):
        uuid,_=DATASETS["JP"]
        x=source_params(uuid,300)
        self.assertEqual(x["datasetKey"],uuid)
        self.assertEqual(x["origin"],"SOURCE")
        self.assertEqual(x["rank"],"SPECIES")
        self.assertEqual(x["offset"],300)
        self.assertEqual(x["limit"],300)
        with self.assertRaises(ValueError):
            source_params(uuid,1)

    def test_registry_provenance_is_mandatory(self):
        uuid,_=DATASETS["JP"]
        a=verify_dataset({"key":uuid,"title":"Global Register of Introduced and Invasive Species - Japan",
                         "type":"CHECKLIST"},uuid,"Japan")
        self.assertIn("Japan",a["title"])
        with self.assertRaisesRegex(ValueError,"Wrong GBIF"):
            verify_dataset({"key":"fake","title":a["title"]},uuid,"Japan")
        with self.assertRaisesRegex(ValueError,"Country"):
            verify_dataset({"key":uuid,"title":"Global Register of Introduced and Invasive Species - Germany"},uuid,"Japan")

    def test_species_page_and_bad_count_fail_closed(self):
        x={"count":2,"results":[{"key":101},{"key":102}],"endOfRecords":True}
        n,rows,end=parse_species_page(x,0)
        self.assertEqual(n,2);self.assertEqual(len(rows),2);self.assertTrue(end)
        with self.assertRaisesRegex(ValueError,"count"):
            parse_species_page({"count":False,"results":[]},0)
        with self.assertRaisesRegex(ValueError,"Empty"):
            parse_species_page({"count":2,"results":[]},0)

    def test_plant_species_only_direct_source_backbone_match(self):
        uuid,_=DATASETS["JP"]
        def fake(_client,url,params=None):
            if url.endswith("/dataset/"+uuid):
                return {"key":uuid,
                        "title":"Global Register of Introduced and Invasive Species - Japan",
                        "type":"CHECKLIST"}
            return {"count":4,"endOfRecords":True,"results":[
                {"key":100,"nubKey":111,"kingdom":"Plantae","rank":"SPECIES"},
                {"key":101,"nubKey":111,"kingdom":"Plantae","rank":"SPECIES"},
                {"key":102,"nubKey":None,"kingdom":"Plantae","rank":"SPECIES"},
                {"key":103,"nubKey":555,"kingdom":"Animalia","rank":"SPECIES"}]}
        out=audit_country(None,"JP",get_json=fake,sleep=lambda _x:None)
        self.assertTrue(out["full_source_pages_verified"])
        self.assertEqual(out["n_indexed_source_species_rank_usage_rows"],4)
        self.assertEqual(out["gbif_backbone_plant_species_key_set"],["111"])
        self.assertEqual(out["n_unmatched_plant_name_usages"],1)
        self.assertEqual(out["n_nonplantae"],1)
        self.assertTrue(out["not_in_list_means_unknown_not_native"])

    def test_duplicate_source_rows_are_not_independent_taxa(self):
        uuid,_=DATASETS["JP"]
        def fake(_client,url,params=None):
            if "/dataset/" in url:
                return {"key":uuid,
                        "title":"Global Register of Introduced and Invasive Species - Japan"}
            return {"count":2,"endOfRecords":True,"results":[
                {"key":101,"nubKey":111,"kingdom":"Plantae","rank":"SPECIES"},
                {"key":101,"nubKey":222,"kingdom":"Plantae","rank":"SPECIES"}]}
        with self.assertRaisesRegex(ValueError,"Duplicated"):
            audit_country(None,"JP",get_json=fake,sleep=lambda _x:None)


if __name__=="__main__":
    unittest.main()
