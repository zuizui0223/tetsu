"""Synthetic safety gates for near-track pixel eligibility (no biological data)."""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from tetsu.rail_verge_pixel_eligibility import (
    near_points_from_140m,preconstruction_comparability,analyze_pairs
)


class RailCenterlinePixelEligibility(unittest.TestCase):
    def test_perpendicular_interpolation_returns_exact_20m_35m(self):
        pair=[("T00_L_near","T00_L_control")]
        original={"T00_L_near":(140.,0.),
                  "T00_L_control":(1600.,0.)}
        for offset in (20,35):
            pts,pairs=near_points_from_140m(original,pair,offset)
            self.assertAlmostEqual(pts[pairs[0][0]][0],offset,places=7)
            self.assertAlmostEqual(pts[pairs[0][1]][0],1600,places=7)
        with self.assertRaises(ValueError):
            near_points_from_140m(original,pair,0)
        with self.assertRaises(ValueError):
            near_points_from_140m(original,pair,140)

    def test_current_road_near_does_not_imply_rail_vegetation(self):
        flags=preconstruction_comparability(
            6.,5.,4,4,.32,.36
        )
        self.assertTrue(flags["2017_both_vegetation"])
        self.assertTrue(flags["2024_road_distances_matched"])
        self.assertTrue(flags["2024_near_road_lt25m"])
        self.assertFalse(flags["strict_veg_road_isolated_and_road_matched"])

    def test_strict_road_and_preenvironment_gate(self):
        flags=preconstruction_comparability(
            40.,35.,4,4,.33,.37
        )
        self.assertTrue(flags["strict_veg_road_isolated_and_road_matched"])
        far_bad=preconstruction_comparability(40.,135.,4,4,.33,.37)
        self.assertFalse(far_bad["2024_road_distances_matched"])
        ndvi_bad=preconstruction_comparability(40.,35.,4,4,.33,.75)
        self.assertFalse(ndvi_bad["strict_veg_road_isolated_and_road_matched"])

    def test_non_vegetated_preconstruction_cannot_qualify(self):
        flags=preconstruction_comparability(40.,35.,5,5,.03,.05)
        self.assertTrue(flags["2017_scl_both_clear"])
        self.assertFalse(flags["2017_both_vegetation"])
        self.assertFalse(flags["strict_veg_road_isolated_and_road_matched"])

    def test_source_missing_zero_not_absence_and_no_did(self):
        pairs=[("T00_L_near20m","T00_L_control")]
        road={"T00_L_near20m":50.,"T00_L_control":60.}
        result=analyze_pairs(pairs,{2024:{"T00_L_near20m":.30}},
                             {"T00_L_near20m":.25},
                             {"T00_L_near20m":4},road)
        self.assertEqual(result["n_geographic_pairs"],1)
        self.assertEqual(result["pixel_NDVI_contrasts"]["2024"]["all"]["n_complete"],0)
        self.assertIsNone(result["pixel_NDVI_contrasts"]["2024"]["all"]["median"])

    def test_complete_pair_true_spectral_but_not_independent_treatment(self):
        pairs=[("T00_L_near20m","T00_L_control"),
               ("T15_R_near20m","T15_R_control")]
        source={
            "T00_L_near20m":.40,"T00_L_control":.45,
            "T15_R_near20m":.30,"T15_R_control":.35
        }
        post={
            2024:{
                "T00_L_near20m":.30,"T00_L_control":.47,
                "T15_R_near20m":.40,"T15_R_control":.40
            }
        }
        road={k:35. for k in source}
        scl={k:4 for k in source}
        out=analyze_pairs(pairs,post,source,scl,road)
        x=out["pixel_NDVI_contrasts"]["2024"]["strict_veg_road_isolated_and_road_matched"]
        self.assertEqual(x["n_complete"],2)
        self.assertEqual(x["n_contiguous_along_rail_blocks"],2)
        self.assertTrue(x["one_railway_corridor_not_independent_experiments"])


if __name__=="__main__":
    unittest.main()
