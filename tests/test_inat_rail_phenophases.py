"""Unit tests for evidence-preserving observation quality rules.

Tests are synthetic fixtures; no source API or train measurements used.
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.audit_inat_rail_phenophases import acceptable_observation, near_rail_bin


def example():
    return {
        "id": 42,
        "geojson": {"type": "Point", "coordinates": [139.7, 35.7]},
        "obscured": False, "geoprivacy": None, "captive": False,
        "positional_accuracy": 15, "observed_on": "2024-05-20",
        "photos": [{"id": 1, "license_code": "cc-by"}]
    }


class TestInatPhenophaseQC(unittest.TestCase):
    def test_high_quality(self):
        self.assertEqual(acceptable_observation(example()), (True, "pass"))

    def test_obscured_rejected_even_when_coords_present(self):
        item = example()
        item["obscured"] = True
        self.assertEqual(acceptable_observation(item)[1], "obscured_or_private")

    def test_cultivation_not_mistaken_for_wild_population(self):
        item = example()
        item["captive"] = True
        self.assertEqual(acceptable_observation(item)[1], "captive_or_cultivated")

    def test_accuracy_missing_or_large_rejected(self):
        item = example()
        item["positional_accuracy"] = None
        self.assertEqual(acceptable_observation(item)[1], "accuracy_not_within_100m")
        item["positional_accuracy"] = 500
        self.assertEqual(acceptable_observation(item)[1], "accuracy_not_within_100m")

    def test_photos_mandatory(self):
        item = example()
        item["photos"] = []
        self.assertEqual(acceptable_observation(item)[1], "no_photos")

    def test_bins(self):
        self.assertEqual(near_rail_bin(250), "0-250m")
        self.assertEqual(near_rail_bin(251), "250-1000m")
        self.assertEqual(near_rail_bin(1000), "250-1000m")
        self.assertEqual(near_rail_bin(1001), ">1000m")


if __name__ == "__main__":
    unittest.main()
