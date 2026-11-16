"""Fully synthetic GBIF region x cell taxon-facet QA tests; never actual ecological evidence."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from harvest_gbif_cell_species_facets_v20 import (
    grid_cells,polygon,query_params,parse_page,full_species_facet,
    compute_cell_period,comparison,PERIODS,FACET_PAGE
)


class GBIFCompleteFacetTests(unittest.TestCase):
    def test_nine_fixed_nonempty_geography_cells(self):
        cells=grid_cells("jp_sendai")
        self.assertEqual(len(cells),9)
        self.assertEqual(len({x["cell_id"] for x in cells}),9)
        self.assertEqual({x["country"] for x in cells},{"JP"})
        self.assertEqual(cells[0]["cell_id"],"jp_sendai_x0_y0")
        with self.assertRaisesRegex(ValueError,"frozen"):
            grid_cells("jp_sendai",2)

    def test_source_geometry_and_period_are_fixed(self):
        c=grid_cells("de_leipzig")[0]
        q=query_params(c,PERIODS[1],"speciesKey",200)
        self.assertEqual(q["facet"],"speciesKey")
        self.assertEqual(q["facetOffset"],200)
        self.assertEqual(q["facetLimit"],5000)
        self.assertEqual(q["speciesKey.facetLimit"],5000)
        self.assertEqual(q["speciesKey.facetOffset"],200)
        self.assertEqual(q["country"],"DE")
        self.assertEqual(q["year"],"2018,2025")
        self.assertEqual(q["limit"],0)
        self.assertEqual(q["occurrenceStatus"],"PRESENT")
        self.assertIn("POLYGON((",q["geometry"])
        with self.assertRaisesRegex(ValueError,"Unbounded"):
            query_params(c,PERIODS[1],"speciesKey",0,5001)

    def test_uppercase_gbif_facet_response(self):
        n,items=parse_page({"count":1000,"facets":[{
            "field":"SPECIES_KEY",
            "counts":[{"name":"123","count":300},{"name":"456","count":100}]
        }]},"speciesKey")
        self.assertEqual(n,1000)
        self.assertEqual(items,[("123",300),("456",100)])
        self.assertEqual(parse_page({"count":0,"facets":[]},"speciesKey"),(0,[]))

    def test_page_audits_completeness_on_synthetic_5001_categories(self):
        c=grid_cells("jp_sendai")[0]
        allitems=[{"name":str(i+1000),"count":1} for i in range(5001)]
        queried=[]
        def fake_request(_client,params):
            self.assertEqual(params["facet"],"speciesKey")
            queried.append(params["facetOffset"])
            off=params["facetOffset"]
            return {"count":6000,"facets":[{"field":"SPECIES_KEY",
                      "counts":allitems[off:off+FACET_PAGE]}]}
        x=full_species_facet(None,c,PERIODS[1],request=fake_request,sleep=lambda _:None)
        self.assertEqual(queried,[0,5000])
        self.assertEqual(x["indexed_occurrence_count"],6000)
        self.assertEqual(x["n_species_key_categories"],5001)
        self.assertTrue(x["full_facet_pagination_completed"])
        self.assertTrue(x["not_original_occurrence_rows"])

    def test_repeated_facet_page_fails_instead_of_truncating(self):
        c=grid_cells("jp_sendai")[0]
        row=[{"name":str(i),"count":1} for i in range(5000)]
        def bad(_client,params):
            return {"count":12000,"facets":[{"field":"SPECIES_KEY","counts":row}]}
        with self.assertRaisesRegex(ValueError,"Repeated facet"):
            full_species_facet(None,c,PERIODS[0],request=bad,sleep=lambda _:None)

    def test_taxon_facets_never_invent_absence(self):
        c=grid_cells("jp_sendai")[0]
        def mock(_client,params):
            if params["facet"]=="speciesKey":
                return {"count":301,"facets":[{"field":"SPECIES_KEY","counts":[
                    {"name":"111","count":200},{"name":"222","count":100}]}]}
            return {"count":301,"facets":[{"field":"DATASET_KEY","counts":[
                    {"name":str(i),"count":25} for i in range(4)]}]}
        x=compute_cell_period(None,c,PERIODS[1],request=mock,sleep=lambda _:None)
        self.assertTrue(x["source_eligible_for_recording_effort_screen"])
        self.assertEqual(x["n_source_datasets_at_least"],4)
        self.assertTrue(x["species_keys_are_observed_only_not_biological_absence"])

    def test_comparison_requires_source_eligibility_in_both_periods(self):
        a={"source_eligible_for_recording_effort_screen":True,
           "species_key_occurrence_counts":{"1":5,"2":8}}
        b={"source_eligible_for_recording_effort_screen":True,
           "species_key_occurrence_counts":{"2":8,"3":7}}
        p=comparison({PERIODS[0]:a,PERIODS[1]:b})
        self.assertEqual(p["indexed_taxon_jaccard"],.333333)
        self.assertEqual(p["n_shared_observed_species_keys"],1)
        b["source_eligible_for_recording_effort_screen"]=False
        q=comparison({PERIODS[0]:a,PERIODS[1]:b})
        self.assertIsNone(q["indexed_taxon_jaccard"])
        self.assertFalse(q["eligible"])

    def test_count_and_category_integrity(self):
        bad=[
            {"count":10,"facets":[{"field":"SPECIES_KEY","counts":[{"name":"12","count":-1}]}]},
            {"count":10,"facets":[{"field":"SPECIES_KEY","counts":[{"name":"12","count":True}]}]},
            {"count":10,"facets":[]},
            {"count":False,"facets":[]},
            {"count":10,"facets":[{"field":"SPECIES_KEY","counts":[{"name":"12","count":2},{"name":"12","count":3}]}]},
        ]
        for x in bad:
            with self.assertRaises(ValueError):
                parse_page(x,"speciesKey")


if __name__=="__main__":
    unittest.main()
