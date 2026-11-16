"""Offline country-introduced speciesKeys on pre-fixed equal-area macro cells."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from harvest_griis_macro_50km_facets_v31 import (
    species_params,source_species_paginated,run
)
from audit_gbif_equalarea_macroblocks_v29 import make_grid,PERIODS
from analyze_macro_rail_assemblages_v22 import COUNTRIES


def fake_griis():
    c=sorted(set(COUNTRIES.values())-{"TW"})
    return {
        "complete":True,"errors":{},"countries":{
            country:{
                "country":country,"full_source_pages_verified":True,
                "n_directly_matched_plant_backbone_species_keys":110,
                "gbif_backbone_plant_species_key_set":[str(100+i) for i in range(110)]
            } for country in c
        }
    }


def source_coverage():
    old={}
    for cell in make_grid("jp_tohoku"):
        old[cell["cell_id"]]={
            "country":"JP",
            "period_coverage":{
                f"{p[0]}_{p[1]}":{"source_eligible":True,
                              "indexed_present_plant_record_count":105}
                for p in PERIODS
            }
        }
    return {"complete":True,"source_errors":{},
            "blocks":{"jp_tohoku":{"cells":old}}}


def fake_fetch(_,params):
    if params["facet"]!="speciesKey":
        raise AssertionError("Wrong source taxon facet request")
    return {"count":105,"facets":[{"field":"SPECIES_KEY","counts":[
        {"name":str(100+i),"count":20} for i in range(5)
    ]}]}


class MacroAlienFacetTests(unittest.TestCase):
    def test_country_geo_scope_and_5000_species_facet(self):
        cell=make_grid("jp_tohoku")[0]
        p=species_params(cell,PERIODS[0],0)
        self.assertEqual(p["country"],"JP")
        self.assertEqual(p["year"],"2010,2017")
        self.assertEqual(p["facet"],"speciesKey")
        self.assertEqual(p["limit"],0)
        self.assertEqual(p["speciesKey.facetLimit"],5000)
        self.assertIn("POLYGON((",p["geometry"])
        with self.assertRaises(ValueError):
            species_params(cell,(2019,2025),0)

    def test_5_griis_introduced_keys_vs_unknown_nonlisted(self):
        cell=make_grid("jp_tohoku")[0]
        x=source_species_paginated(None,cell,PERIODS[0],
                                  set(str(100+i) for i in range(110)),
                                  request=fake_fetch,sleep=lambda _:None)
        self.assertEqual(x["n_observed_GRIIS_introduced_plant_species_keys"],5)
        self.assertEqual(x["n_indexed_GRIIS_introduced_plant_records"],100)
        self.assertTrue(x["alien_observed_minimum_5_species_and_20_records"])
        self.assertTrue(x["nonlisted_species_status_unknown_not_native"])

    def test_all_32_cell_periods_exactly_completable_with_source_eligibility(self):
        calls=[]
        def source(client,params):
            calls.append(params)
            return fake_fetch(client,params)
        out=run(["jp_tohoku"],source_coverage(),fake_griis(),client=object(),
                request=source,sleep=lambda _:None)
        self.assertTrue(out["complete"])
        self.assertEqual(len(calls),32)
        self.assertEqual(out["regions"]["jp_tohoku"]["n_complete_introduced_species_facets_two_period_cells"],16)
        self.assertEqual(out["regions"]["jp_tohoku"]["n_introduced_source_minimum_both_period_cells"],16)

    def test_incomplete_checklist_or_bad_prior_coverage_never_becomes_native(self):
        old=source_coverage();old["complete"]=False
        with self.assertRaisesRegex(ValueError,"incomplete"):
            run(["jp_tohoku"],old,fake_griis(),client=object(),
                request=fake_fetch,sleep=lambda _:None)
        old=source_coverage();bad=fake_griis();bad["complete"]=False
        with self.assertRaisesRegex(ValueError,"Incomplete GRIIS"):
            run(["jp_tohoku"],old,bad,client=object(),request=fake_fetch,sleep=lambda _:None)

    def test_unlisted_taxa_are_not_imputed_native_or_present(self):
        cell=make_grid("jp_tohoku")[0]
        x=source_species_paginated(None,cell,PERIODS[0],{"999","100"},
                                  request=fake_fetch,sleep=lambda _:None)
        self.assertEqual(x["n_observed_GRIIS_introduced_plant_species_keys"],1)
        self.assertFalse(x["alien_observed_minimum_5_species_and_20_records"])
        self.assertNotIn("101",x["observed_alien_specieskey_counts"])


if __name__=="__main__":
    unittest.main()
