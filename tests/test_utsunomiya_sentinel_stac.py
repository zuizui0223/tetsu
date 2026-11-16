"""Synthetic STAC catalogue fixtures only; no satellite vegetation effect."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_utsunomiya_sentinel_stac import summarize_items,stac_query_params


def item(year=2024):
    return {"type":"Feature","id":"S2-synthetic-001",
            "properties":{"datetime":f"{year}-07-15T02:10:00Z",
                          "eo:cloud_cover":12},
            "assets":{
                "red":{"href":"https://example.invalid/red.tif"},
                "nir":{"href":"https://example.invalid/nir.tif"},
                "scl":{"href":"https://example.invalid/scl.tif"}
            }}


class SentinelStac(unittest.TestCase):
    def test_fixed_source_time_windows(self):
        x=stac_query_params(2017)
        self.assertEqual(x["datetime"],
                         "2017-06-01T00:00:00Z/2017-09-30T23:59:59Z")
        self.assertEqual(x["query"]["eo:cloud_cover"]["lt"],20)
        with self.assertRaises(ValueError):
            stac_query_params(2016)

    def test_asset_keys_not_asset_dict_values(self):
        result=summarize_items(2024,{
            "type":"FeatureCollection","features":[item()],
            "context":{"matched":1},"links":[]
        })
        self.assertEqual(result["returned_stac_items"],1)
        self.assertEqual(result["number_of_distinct_acquisition_dates"],1)
        self.assertEqual(result["n_items_with_red_nir_and_cloud_assets"],1)
        self.assertEqual(result["stac_asset_keys_observed"]["red"],1)

    def test_missing_cloud_mask_not_analyzable(self):
        x=item()
        x["assets"].pop("scl")
        summary=summarize_items(2024,{"type":"FeatureCollection",
                                       "features":[x],"links":[]})
        self.assertEqual(summary["n_items_with_red_nir_and_cloud_assets"],0)

    def test_unexpected_year_or_cloud_pct_rejected(self):
        with self.assertRaises(ValueError):
            summarize_items(2025,{"type":"FeatureCollection",
                                  "features":[item(year=2024)]})
        x=item()
        x["properties"]["eo:cloud_cover"]=32
        with self.assertRaises(ValueError):
            summarize_items(2024,{"type":"FeatureCollection","features":[x]})

    def test_truncation_reported_not_silently_complete(self):
        summary=summarize_items(2024,{
            "type":"FeatureCollection","features":[item()],
            "links":[{"rel":"next","href":"https://example.invalid/?page=2"}]
        })
        self.assertTrue(summary["is_truncated_or_has_next"])


if __name__=="__main__":
    unittest.main()
