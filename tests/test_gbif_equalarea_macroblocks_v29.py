"""Strict offline audit of 50 km equal-area cross-country macro source sampling."""
import sys
from pathlib import Path
import unittest

from shapely import wkt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_gbif_equalarea_macroblocks_v29 import (
    BLOCKS,PILOT,PERIODS,validate_frame,make_grid,
    query_params,parse_coverage,run
)


class MacroEqualAreaTests(unittest.TestCase):
    def test_24_geographically_nonoverlapping_predeclared_blocks(self):
        areas=validate_frame()
        self.assertEqual(len(areas),24)
        self.assertEqual(len({v[0] for v in BLOCKS.values()}),8)
        self.assertEqual(len(PILOT),10)
        self.assertNotIn("tw_hsinchu",BLOCKS)
        self.assertNotIn("kr_daejeon",BLOCKS)

    def test_16_cells_per_block_nominal_area_2500km2(self):
        cells=make_grid("jp_tohoku")
        self.assertEqual(len(cells),16)
        self.assertEqual(len({v["cell_id"] for v in cells}),16)
        self.assertTrue(all(x["epsg6933_equal_area_km2"]==2500 for x in cells))
        for row in cells:
            geo=wkt.loads(row["gbif_geometry_wgs84"])
            self.assertTrue(geo.is_valid)
            self.assertGreater(geo.area,0)
        with self.assertRaises(ValueError):
            make_grid("not_frozen")

    def test_country_region_and_period_are_not_inferred_from_observations(self):
        c=make_grid("fr_burgundy")[0]
        req=query_params(c,PERIODS[0])
        self.assertEqual(req["country"],"FR")
        self.assertEqual(req["year"],"2010,2017")
        self.assertEqual(req["taxonRank"],"SPECIES")
        self.assertEqual(req["occurrenceStatus"],"PRESENT")
        self.assertEqual(req["limit"],0)
        self.assertEqual(req["facet"],"datasetKey")
        self.assertEqual(req["datasetKey.facetLimit"],4)
        self.assertIn("POLYGON((",req["geometry"])
        with self.assertRaises(ValueError):
            query_params(c,(2019,2025))

    def test_source_has_at_least_three_datasets_not_single_dominant_dataset(self):
        payload={"count":125,"facets":[{"field":"DATASET_KEY","counts":[
            {"name":"dataset_a","count":110},{"name":"dataset_b","count":10},
            {"name":"dataset_c","count":5}]}]}
        r=parse_coverage(payload)
        self.assertTrue(r["source_eligible"])
        self.assertEqual(r["n_dataset_categories_at_least"],3)
        payload["facets"][0]["counts"]=payload["facets"][0]["counts"][:1]
        self.assertFalse(parse_coverage(payload)["source_eligible"])

    def test_no_source_records_does_not_become_plant_absence(self):
        r=parse_coverage({"count":0,"facets":[]})
        self.assertFalse(r["source_eligible"])
        self.assertEqual(r["indexed_present_plant_record_count"],0)
        for payload in ({"count":None,"facets":[]},{"count":True,"facets":[]},
                        {"count":100,"facets":[]},
                        {"count":100,"facets":[{"field":"DATASET_KEY","counts":[{"name":"x","count":False}]}]}):
            with self.assertRaises(ValueError):
                parse_coverage(payload)

    def test_source_query_audit_uses_actual_32_cell_period_queries(self):
        called=[]
        def fake(_session,url,params):
            called.append((url,params))
            return {"count":200,"facets":[{"field":"DATASET_KEY","counts":[
                {"name":"d1","count":50},{"name":"d2","count":50},{"name":"d3","count":100}]}]}
        source=run(["jp_tohoku"],source_api=fake,sleep=lambda _:None)
        self.assertEqual(len(called),32)
        self.assertTrue(source["complete"])
        self.assertEqual(source["summary"]["n_equal_area_source_cells"],16)
        self.assertEqual(source["summary"]["n_cells_with_both_period_coverage_gate"],16)
        self.assertTrue(source["interpretation_gates"])

    def test_empty_block_list_and_unfrozen_country_fail(self):
        with self.assertRaises(ValueError):
            run([],source_api=lambda *a:None,sleep=lambda _:None)
        with self.assertRaises(ValueError):
            run(["fake_region"],source_api=lambda *a:None,sleep=lambda _:None)


if __name__=="__main__":
    unittest.main()
