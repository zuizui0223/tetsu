"""Real source acquisition is performed in GitHub Actions, NOT these fixtures."""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_japan_ordinary_rail_2025_s2_vegetation import (
    validate_stac_tile,choose_source_dates,N_ACQUISITION_DATES
)


def fixture(day,cloud=18.,mask=True,offset=-.1):
    bands=[{"scale":.0001,"offset":offset}]
    a={
        "red":{"href":"https://sentinel-cogs.s3.amazonaws.com/red.tif",
               "raster:bands":bands},
        "nir":{"href":"https://sentinel-cogs.s3.amazonaws.com/nir.tif",
               "raster:bands":bands},
    }
    if mask:
        a["scl"]={"href":"https://sentinel-cogs.s3.amazonaws.com/scl.tif"}
    return {
        "type":"Feature",
        "id":f"synthetic-S2-{day}",
        "properties":{"datetime":f"{day}T01:20:00Z",
                      "eo:cloud_cover":cloud},
        "assets":a
    }


class OrdinaryRailCurrentS2Tests(unittest.TestCase):
    def test_valid_reprocessed_reflectance_requires_boa_metadata(self):
        self.assertTrue(validate_stac_tile(fixture("2025-08-01"),2025))
        self.assertFalse(validate_stac_tile(fixture("2025-08-01",mask=False),2025))
        with self.assertRaises(ValueError):
            validate_stac_tile(fixture("2025-08-01",offset=0),2025)

    def test_2025_source_date_selected_without_ndvi_or_vegetation_result(self):
        items=[
            fixture("2025-08-21",5),
            fixture("2025-07-20",1),
            fixture("2025-08-01",77),
            fixture("2025-07-31",77),
            fixture("2025-07-30",1),
        ]
        d,skips=choose_source_dates(items)
        self.assertEqual(len(d),N_ACQUISITION_DATES)
        self.assertEqual(next(iter(d)),"2025-08-01")
        self.assertNotIn("2025-08-21",d)

    def test_duplicate_date_tiles_kept_and_missing_cloudmask_rejected(self):
        items=[fixture("2025-08-03",cloud=3),
               fixture("2025-08-03",cloud=8),
               fixture("2025-08-03",mask=False)]
        # Fixtures with same ID still occupy separate possible tiles in the
        # metadata-only test; no raw imagery required.
        d,skip=choose_source_dates(items)
        self.assertEqual(len(d["2025-08-03"]),2)
        self.assertEqual(skip["not_all_red_nir_scl_assets"],1)

    def test_unexpected_external_cog_host_rejected(self):
        f=fixture("2025-08-01")
        f["assets"]["nir"]["href"]="https://untrusted.example/test.tif"
        with self.assertRaises(ValueError):
            validate_stac_tile(f,2025)


if __name__=="__main__":
    unittest.main()
