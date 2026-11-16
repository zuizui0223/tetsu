"""v10.5 source-independent mathematical tests (NOT a replacement for original run)."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_Slovenia_rail_v105_cluster_missingness import (
    percentile, segment_cluster_diagnostic,
    paired_fraction_extreme_bounds, pair_missingness_bounds,
    original_pair_wind_count, section_and_missingness_report,
    verify_against_locked_v104, WIND_STRA,
)


def sourcepair(i, side, section, wind_delta=-.05, known=10, missing=2, nwind_bs=2, nwind_bt=3):
    return {
        "location_id": "L" + str(i),
        "track_side": side,
        "section_number": section,
        "n_original_taxa_known_strategy_BS": known,
        "n_original_taxa_known_strategy_BT": known,
        "n_original_taxa_missing_strategy_BS": missing,
        "n_original_taxa_missing_strategy_BT": missing,
        "wind_and_diaspore_group_shares_BS_minus_BT": {WIND_STRA: wind_delta},
        "_wind_original_counts": {"BS": nwind_bs, "BT": nwind_bt},
    }


class SegmentAndMissingTraitTests(unittest.TestCase):
    def test_quantile_interpolation(self):
        self.assertEqual(percentile([1, 3, 5], .5), 3)
        self.assertEqual(percentile([1, 3], .25), 1.5)
        self.assertIsNone(percentile([], .5))

    def test_railway_section_is_cluster_not_100_sites(self):
        rows = ([{"section_number": 1, "d": 1.0}] * 100 +
                [{"section_number": 2, "d": -1.0}] * 2)
        x = segment_cluster_diagnostic(rows, "d", draws=1000)
        self.assertEqual(x["n_source_segments"], 2)
        self.assertEqual(x["n_original_field_locations"], 102)
        self.assertAlmostEqual(x["equal_field_location_mean_BS_minus_BT"], 98/102, places=5)
        self.assertEqual(x["equal_rail_section_mean_BS_minus_BT"], 0)
        self.assertEqual(x["n_leave_one_section_out_positive_equal_location"], 1)
        self.assertEqual(x["n_leave_one_section_out_negative_equal_location"], 1)

    def test_unknown_trait_bounds_can_reverse_naive_wind_result(self):
        # BS observed 1/3 known wind, 1 unknown, BT 1/4 known wind.
        lo, hi = paired_fraction_extreme_bounds(1, 3, 1, 1, 4, 0)
        self.assertEqual(lo, 0)
        self.assertEqual(hi, 0.25)
        # No all-taxa directional conclusion can be made in this case.
        self.assertTrue(lo <= 0 <= hi)
        with self.assertRaises(ValueError):
            paired_fraction_extreme_bounds(5, 3, 1, 1, 4, 0)

    def test_source_missingness_is_not_coded_as_nonwind(self):
        p = sourcepair(1, "L", 1, known=4, missing=6, nwind_bs=1, nwind_bt=2)
        self.assertAlmostEqual(pair_missingness_bounds(p)["BS"]["coverage"], .4)
        self.assertEqual(original_pair_wind_count(p, "BS"), 1)
        p["_wind_original_counts"]["BS"] = 5
        with self.assertRaises(ValueError):
            original_pair_wind_count(p, "BS")

    def test_first_average_actual_track_sides_within_each_site(self):
        # Synthetic exactly 89 sites, 173 site-sides and 15 source sections.
        sites, pairs = [], []
        for i in range(89):
            seg = i % 15 + 1
            sites.append({"location_id": "L" + str(i), "section_number": seg})
            pairs.append(sourcepair(i, "L", seg, wind_delta=-.1))
            if i < 84:
                pairs.append(sourcepair(i, "R", seg, wind_delta=.02))
        r = section_and_missingness_report(pairs, sites)
        self.assertEqual(r["strategy_complete_case_sensitivity_original_sites"]["0.0"][
            "n_original_locations_with_qualifying_sides"], 89)
        self.assertEqual(r["strategy_complete_case_sensitivity_original_sites"]["0.0"][
            "number_of_original_side_pairs_qualifying"], 173)
        self.assertEqual(r["strategy_complete_case_sensitivity_original_sites"]["1.0"][
            "n_original_locations_with_qualifying_sides"], 0)
        self.assertEqual(r["mean_sign_is_identified_under_every_unknown_trait_assignment"], "not_identified")

    def test_source_receipt_is_mandatory(self):
        with self.assertRaisesRegex(ValueError, "verified original-source receipt"):
            verify_against_locked_v104({}, {}, {"complete": False})


if __name__ == "__main__":
    unittest.main()
