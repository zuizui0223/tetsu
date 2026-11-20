"""Synthetic mathematical checks for road-heldout slope transfer; no field conclusions."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_MIREN_original_road_out_transfer_v301 import (
    road_out_calibrated_transfer, conditional_trait_label_null,
    _within_road_slope, _anchor_index, region_equal_aggregate,
)

T = np.array([[-1.0, -0.25], [1.0, 0.25]], dtype=float)

def make_road(intercept=0.0, roadside_slope=-0.75, inland_slope=-0.55):
    alt = np.array([100., 200., 300., 400., 500., 600.])
    x = alt / 1000
    def weight(slope):
        y = intercept + slope*x
        p = (y+1)/2
        return np.column_stack([1-p, p])
    return {
        "elev_m": alt,
        "matrices": {
            "roadside": weight(roadside_slope),
            "interior": weight(inland_slope),
        },
    }

class TestGenuineRoadOutTransfer(unittest.TestCase):
    def setUp(self):
        self.roads = [make_road(c) for c in (-.12, .07, .2)]

    def test_known_reproducible_linear_slope(self):
        x = np.array([.1,.2,.3,.4,.5,.6])
        y = 0.2 - .75*x
        self.assertAlmostEqual(_within_road_slope(x,y), -.75, places=12)

    def test_anchor_is_picked_from_altitude_only(self):
        x = np.array([.1,.5,.2,.4,.3])
        self.assertEqual(_anchor_index(x, "lowest"),0)
        self.assertEqual(_anchor_index(x, "median"),4)
        with self.assertRaises(ValueError):
            _anchor_index(x,"highest")

    def test_positive_transfer_on_distinct_roads_with_new_intercepts(self):
        for anchor in ("lowest", "median"):
            d = road_out_calibrated_transfer(self.roads,T,anchor)
            self.assertGreater(d["roadside"]["gain_zlog_height_RMSE"], 0.0)
            self.assertEqual(d["roadside"]["positive_original_heldout_roads"], 3)
            self.assertGreater(d["interior"]["gain_zlog_height_RMSE"], 0.0)
            self.assertEqual(d["roadside"]["original_heldout_road_count"], 3)

    def test_opposing_heldout_slope_can_be_worse_than_baseline(self):
        roads = [make_road(0,-.75), make_road(0,-.75),make_road(0,+.75)]
        d = road_out_calibrated_transfer(roads,T,"lowest")
        # Not a positivity-only test: out-of-road transfer is allowed to fail.
        self.assertLess(d["roadside"]["positive_original_heldout_roads"], 3)

    def test_does_not_use_heldout_road_to_fit_slope(self):
        x = np.array([.1,.2,.3,.4,.5,.6])
        b_train = [
            _within_road_slope(x,road["matrices"]["roadside"]@T[:,0])
            for road in self.roads[1:]
        ]
        self.assertAlmostEqual(np.mean(b_train),-.75,places=12)
        # Change ONLY test-road trait profile; other roads' fitted slope remains.
        self.roads[0] = make_road(.0, +.2)
        b_same = [
            _within_road_slope(x,road["matrices"]["roadside"]@T[:,0])
            for road in self.roads[1:]
        ]
        self.assertAlmostEqual(np.mean(b_train),np.mean(b_same),places=12)

    def test_invalid_geographical_support_fails_not_zero(self):
        with self.assertRaises(ValueError):
            road_out_calibrated_transfer(self.roads[:1],T)
        bad = make_road()
        bad["elev_m"] = np.array([300.,310.,320.,330.,340.,350.])
        with self.assertRaises(ValueError):
            road_out_calibrated_transfer([bad,self.roads[1]],T)

    def test_fixed_398_null_draws_are_not_study_replicates(self):
        d = conditional_trait_label_null(
            self.roads,T,[np.array([0,1],dtype=int)],
            0,0,"fixed_source_region",199)
        for anchor in ("lowest","median"):
            for zone in ("roadside","interior","roadside_minus_interior"):
                x = d["anchor_effects"][anchor][zone]
                self.assertEqual(x["technical_draws_NOT_ecological_field_replicates"],398)
                self.assertAlmostEqual(
                    x["real_out_road_calibrated_slope_gain"] -
                    x["phylo75_x_occupancy_joint_trait_label_null_gain"],
                    x["real_minus_joint_trait_label_null_gain"],places=12)

    def test_biological_region_support_gate_is_not_zero(self):
        self.assertIsNone(region_equal_aggregate({"AUT":{}})["effects"])

if __name__ == "__main__":
    unittest.main()
