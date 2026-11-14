"""Synthetic S2 COG/NDVI QA, no remote imagery or ecology invented."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"scripts"))
from tetsu.sentinel_ndvi import (
    calculate_ndvi, metadata_ambiguous, ndvi_branch_values,
    paired_annual_ndvi,prepost_paired_change
)
from measure_utsunomiya_s2_ndvi_pixels import (
    build_fixed_pairs,choose_year_pixel_items
)
from shapely.geometry import LineString


def item(*, baseline="05.11", boa_applied=True, offset=-0.1):
    bands=[{"nodata":0,"scale":0.0001,"offset":offset}]
    return {
        "properties":{
            "s2:processing_baseline":baseline,
            "earthsearch:boa_offset_applied":boa_applied
        },
        "assets":{"red":{"raster:bands":bands},
                  "nir":{"raster:bands":bands},"scl":{}}
    }


class S2ActualPixelSafetyTests(unittest.TestCase):
    def test_radiometry_metadata_disagreement_is_surfaced(self):
        self.assertTrue(metadata_ambiguous(item()))
        self.assertFalse(metadata_ambiguous(
            item(baseline="00.01",boa_applied=False,offset=0)
        ))

    def test_raw_and_subtract1000_offsets_are_not_same_ndvi(self):
        # Synthetic sample: red 1500 DN and NIR 4500 DN.
        raw=calculate_ndvi(1500,4500,0)
        off=calculate_ndvi(1500,4500,1000)
        self.assertAlmostEqual(raw,.5)
        self.assertAlmostEqual(off,.75)
        self.assertNotEqual(raw,off)
        branches=ndvi_branch_values(1500,4500,item())
        self.assertAlmostEqual(branches["provider_flag_branch"],.5)
        self.assertAlmostEqual(branches["raster_offset_branch"],.75)
        self.assertTrue(branches["calibration_disagreement"])

    def test_nodata_and_bad_denominator_are_not_zero_ndvi(self):
        self.assertIsNone(calculate_ndvi(0,4500))
        self.assertIsNone(calculate_ndvi(1000,1000,1000))
        self.assertIsNone(calculate_ndvi(500,4500,1000))

    def test_frozen_fixed_pixel_pairs_are_not_independent_corridors(self):
        # Synthetic 14.5km straight line in a projected local plane.
        lines=[LineString([(0,0),(14500,0)])]
        sites,pairs,metadata=build_fixed_pairs(lines,n_anchors=9)
        self.assertEqual(len(pairs),18)
        self.assertEqual(len(sites),36)
        self.assertTrue(metadata["not_independent_railway_treatments"])

    def test_prepost_delta_only_on_exact_same_pixel_pairs(self):
        pairs=[("n1","c1"),("n2","c2")]
        data={
            2017:{
                "n1":{"provider_flag_branch":.20},
                "c1":{"provider_flag_branch":.25}
            },
            2024:{
                "n1":{"provider_flag_branch":.45},
                "c1":{"provider_flag_branch":.30},
                "n2":{"provider_flag_branch":.50}
            }
        }
        result=prepost_paired_change(pairs,data,2024,"provider_flag_branch",2017)
        self.assertEqual(result["paired_near_minus_control_delta_NDVI"]["n"],1)
        self.assertAlmostEqual(
            result["paired_near_minus_control_delta_NDVI"]["mean"],.20
        )

    def test_provider_branch_survives_when_stale_raster_offset_invalid(self):
        coords={"p1":(10.,20.)}
        fake={
            "red":{"p1":600},
            "nir":{"p1":2000},
            "scl":{"p1":4}
        }
        with patch(
            "measure_utsunomiya_s2_ndvi_pixels.read_scene_pixels",
            return_value=fake
        ):
            x,quality=choose_year_pixel_items(2024,coords,[item()])
        self.assertEqual(len(x),1)
        self.assertIsNotNone(x["p1"]["provider_flag_branch"])
        self.assertIsNone(x["p1"]["raster_offset_branch"])
        self.assertEqual(x["p1"]["scl_class"],4)
        self.assertEqual(quality["provider_calibration_usable"],1)
        self.assertNotIn("raster_offset_calibration_usable",quality)

    def test_scl_clouds_are_not_clear_even_with_high_ndvi(self):
        coords={"p1":(10.,20.),"p2":(11.,21.)}
        fake={
            "red":{"p1":1500,"p2":1600},
            "nir":{"p1":4500,"p2":4500},
            "scl":{"p1":4,"p2":9}
        }
        with patch(
            "measure_utsunomiya_s2_ndvi_pixels.read_scene_pixels",
            return_value=fake
        ):
            x,quality=choose_year_pixel_items(2024,coords,[item()])
        self.assertEqual(len(x),1)
        self.assertAlmostEqual(x["p1"]["provider_flag_branch"],.5)
        self.assertIn("scl_rejected_9",quality)
        self.assertNotIn("p2",x)


if __name__=="__main__":
    unittest.main()
