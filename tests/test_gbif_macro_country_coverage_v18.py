"""Synthetic tests for GBIF COUNTRY counts, never biological railway effects."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_gbif_macro_country_coverage_v18 import (
    COUNTRIES, PERIODS, BASIS, BASE, TAXON_SOURCE,
    validate_frame, query_params, parse_count, audit,
)


class MacroGBIFCoverageTests(unittest.TestCase):
    def test_frozen_geographic_temporal_scope(self):
        validate_frame()
        self.assertEqual(len(COUNTRIES),10)
        self.assertEqual(len(set(COUNTRIES)),10)
        self.assertEqual(PERIODS,((2010,2017),(2018,2025)))
        self.assertEqual(BASIS,(None,"HUMAN_OBSERVATION","PRESERVED_SPECIMEN"))
        with self.assertRaisesRegex(ValueError,"contiguous"):
            validate_frame(COUNTRIES,((2010,2017),(2019,2025)))
        with self.assertRaisesRegex(ValueError,"country"):
            validate_frame(("JP","JP","US","FR","DE"),PERIODS)

    def test_count_query_not_capped_by_search_page_size(self):
        p=query_params("JP",(2018,2025),None)
        self.assertEqual(p["limit"],0)
        self.assertEqual(p["year"],"2018,2025")
        self.assertEqual(p["kingdomKey"],6)
        self.assertEqual(p["taxonRank"],"SPECIES")
        self.assertEqual(p["hasCoordinate"],"true")
        self.assertEqual(p["hasGeospatialIssue"],"false")
        self.assertNotIn("basisOfRecord",p)
        self.assertEqual(query_params("JP",PERIODS[0],"HUMAN_OBSERVATION")["basisOfRecord"],"HUMAN_OBSERVATION")
        with self.assertRaisesRegex(ValueError,"frozen"):
            query_params("CN",PERIODS[0],None)

    def test_counts_are_full_server_counts_not_100k_paginated_results(self):
        self.assertEqual(parse_count({"count":12345678,"offset":0,"limit":0}),12345678)
        for bad in ({"count":False},{"count":-1},{"count":"123"},{"count":None},{"limit":30,"count":5},{"count":4,"offset":20}):
            with self.assertRaises(ValueError):
                parse_count(bad)

    def test_synthetic_mock_runs_all_20_source_cells(self):
        queries=[]
        def fake_request(_session,url,params=None):
            queries.append((url,params))
            if url==TAXON_SOURCE:
                return {"canonicalName":"Plantae"}
            self.assertEqual(url,BASE)
            basis=(params or {}).get("basisOfRecord")
            counts={None:5000000,"HUMAN_OBSERVATION":3500000,"PRESERVED_SPECIMEN":1300000}
            return {"count":counts[basis],"limit":0,"offset":0}
        out=audit(None,request_func=fake_request,sleep_func=lambda _:None)
        self.assertTrue(out["full_coverage"])
        self.assertEqual(out["n_successful_country_periods"],20)
        self.assertEqual(len(queries),61)
        self.assertEqual(out["errors"],{})
        self.assertEqual(out["counts"]["JP"][1]["by_basis"]["ALL_BASIS"],5000000)
        self.assertTrue(out["species_level_taxon_status_not_audited"])

    def test_backbone_must_be_verified_from_source(self):
        def bad(_session,url,params=None):
            return {"canonicalName":"Animalia"}
        with self.assertRaisesRegex(ValueError,"Plantae"):
            audit(None,request_func=bad,sleep_func=lambda _:None)

    def test_bad_basis_count_does_not_get_reported_as_success(self):
        def fake(_session,url,params=None):
            if url==TAXON_SOURCE:return {"canonicalName":"Plantae"}
            basis=params.get("basisOfRecord")
            return {"count":100 if basis is None else 80,"limit":0,"offset":0}
        out=audit(None,request_func=fake,sleep_func=lambda _:None)
        self.assertFalse(out["full_coverage"])
        self.assertEqual(out["n_successful_country_periods"],0)
        self.assertEqual(len(out["errors"]),20)


if __name__=="__main__":
    unittest.main()
