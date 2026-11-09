"""Synthetic-source tests for first Overture × GBIF JOIN feasibility, no biological effects."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_overture_gbif_overlap_v19 import (
    REGIONS, RAIL_CLASSES, ROAD_CLASSES, RELEASE,
    validate_regions, sql_bbox, gbif_query_params, parse_gbif_facet
)


class OvertureGBIFSourceTest(unittest.TestCase):
    def test_country_and_replicated_geography_frame(self):
        validate_regions()
        self.assertEqual(len(REGIONS),10)
        self.assertEqual(len({x[0] for x in REGIONS.values()}),10)
        self.assertIn("jp_sendai",REGIONS)
        self.assertIn("de_leipzig",REGIONS)

    def test_actual_source_parquet_geometry_query_scope(self):
        query=sql_bbox(None,REGIONS["jp_sendai"])
        self.assertIn("read_parquet(",query)
        self.assertIn("2026-10-21.0",query)
        self.assertIn("bbox.xmin <= 141.01",query)
        self.assertIn("bbox.ymax >= 38.16",query)
        self.assertIn("subtype IN ('rail','road')",query)
        self.assertNotIn("LIMIT 1",query)

    def test_rail_vs_road_classes_are_distinct(self):
        self.assertIn("standard_gauge",RAIL_CLASSES)
        self.assertNotIn("subway",RAIL_CLASSES)
        self.assertIn("residential",ROAD_CLASSES)
        self.assertNotIn("footway",ROAD_CLASSES)

    def test_gbif_region_query_is_geographic_and_bounded(self):
        p=gbif_query_params(REGIONS["de_leipzig"])
        self.assertIn("POLYGON((",p["geometry"])
        self.assertEqual(p["country"],"DE")
        self.assertEqual(p["limit"],0)
        self.assertEqual(p["year"],"2018,2025")
        self.assertEqual(p["facet"],"speciesKey")
        self.assertEqual(p["speciesKey.facetLimit"],50)
        self.assertEqual(p["hasCoordinate"],"true")
        self.assertEqual(p["hasGeospatialIssue"],"false")

    def test_facets_never_counted_as_complete_species_list(self):
        s={"count":5634,"facets":[{"field":"speciesKey","counts":[
            {"name":"1345","count":120},{"name":"999","count":12}]}]}
        parsed=parse_gbif_facet(s)
        self.assertEqual(parsed["source_record_count"],5634)
        self.assertEqual(parsed["n_top_categories_returned"],2)
        self.assertFalse(parsed["exhaustive_species_inventory"])
        self.assertEqual(parsed["top50_species_key_categories"][0]["n_indexed_occurrences"],120)

    def test_unavailable_or_malformed_facets_fail_closed(self):
        for payload in (
            {"count":12},{"count":12,"facets":[]},
            {"count":12,"facets":[{"field":"familyKey","counts":[]}]},
            {"count":True,"facets":[{"field":"speciesKey","counts":[]}]},
            {"count":-1,"facets":[{"field":"speciesKey","counts":[]}]},
            {"count":12,"facets":[{"field":"speciesKey","counts":[{"name":"abc","count":-2}]}]}
        ):
            with self.assertRaises(ValueError):
                parse_gbif_facet(payload)


if __name__=="__main__":
    unittest.main()
