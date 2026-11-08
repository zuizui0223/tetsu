"""Synthetic calibration invariants; real pixel values are in CI source audit."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"src"))
from verify_utsunomiya_s2_c1_pixels import (
    collection1_calibrated_pixels,same_scene_calibration_differences
)
from tetsu.sentinel_ndvi import calculate_ndvi


class CrossCollectionPhysicalCalibrationTests(unittest.TestCase):
    def test_c1_metadata_minus1000_matches_legacy_already_corrected_dn(self):
        # Same physical signal, two different storage conventions.
        c1=calculate_ndvi(1600,4500,subtract_dn=1000)
        legacy=calculate_ndvi(600,3500,subtract_dn=0)
        self.assertAlmostEqual(c1,legacy)

    def test_cloud_shadow_and_water_excluded_from_reprocessed_C1(self):
        points={"a":(1,2),"b":(3,4),"c":(5,6)}
        mocked={
            "red":{"a":1600,"b":1600,"c":1600},
            "nir":{"a":4500,"b":4500,"c":4500},
            "scl":{"a":4,"b":9,"c":6}
        }
        with patch(
            "verify_utsunomiya_s2_c1_pixels.read_scene_pixels",
            return_value=mocked
        ):
            values,qc=collection1_calibrated_pixels(points,[{}])
        self.assertEqual(set(values),{"a"})
        self.assertIn("scl_9_excluded",qc)
        self.assertIn("scl_6_excluded",qc)

    def test_cross_collection_concordance_requires_four_actual_pairs(self):
        c1={f"p{x}":{"harmonized":.3} for x in range(4)}
        old={f"p{x}":{"provider_flag_branch":.31} for x in range(4)}
        r=same_scene_calibration_differences(old,c1)
        self.assertEqual(r["same_coordinate_same_date_comparable_valid_points"],4)
        self.assertEqual(r["calibration_cross_source_consistency_gate"],"pass")
        self.assertAlmostEqual(r["absolute_delta_ndvi"]["median"],.01)
        missing=same_scene_calibration_differences(
            {"p0":{"provider_flag_branch":.31}},c1
        )
        self.assertEqual(missing["calibration_cross_source_consistency_gate"],
                         "fail_insufficient_or_disagreeing_source")

    def test_cross_collection_disagreement_does_not_become_vegetation_effect(self):
        c1={f"p{x}":{"harmonized":.6} for x in range(4)}
        old={f"p{x}":{"provider_flag_branch":.1} for x in range(4)}
        r=same_scene_calibration_differences(old,c1)
        self.assertEqual(r["calibration_cross_source_consistency_gate"],
                         "fail_insufficient_or_disagreeing_source")


if __name__=="__main__":
    unittest.main()
