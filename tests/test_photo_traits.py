"""Tests use generated synthetic masks, never observational plant data."""
import sys
import unittest
from pathlib import Path

from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.photo_traits import summarize_masks


class PhotoTraitsTest(unittest.TestCase):
    def setUp(self):
        self.plant = Image.new("L", (100, 100), 0)
        painter = ImageDraw.Draw(self.plant)
        painter.rectangle((20, 20, 59, 79), fill=255)
        self.flower = Image.new("L", (100, 100), 0)
        ImageDraw.Draw(self.flower).rectangle((30, 20, 49, 29), fill=255)

    def test_solid_rectangle_solidity_and_flower_location(self):
        r = summarize_masks(self.plant, self.flower, whole_plant=True, standardized_view=True)
        self.assertEqual(r["quality_reasons"], "pass")
        self.assertTrue(r["valid_for_cross_image_shape_comparison"])
        self.assertEqual(r["plant_area_px"], 2400)
        self.assertAlmostEqual(r["plant_frame_fraction"], 0.24)
        self.assertAlmostEqual(r["plant_width_height_ratio"], 2 / 3, places=5)
        self.assertAlmostEqual(r["plant_silhouette_solidity"], 1, places=5)
        self.assertAlmostEqual(r["flower_area_over_plant_area"], 200 / 2400, places=5)
        self.assertGreater(r["flower_relative_vertical_position"], 0.9)

    def test_missing_flower_mask_is_not_absence(self):
        r = summarize_masks(self.plant, whole_plant=True, standardized_view=True)
        self.assertIsNone(r["flower_area_over_plant_area"])
        self.assertFalse(r["flower_mask_annotated"])

    def test_annotated_absence_is_zero(self):
        absent = Image.new("L", (100, 100), 0)
        r = summarize_masks(self.plant, absent, whole_plant=True, standardized_view=True)
        self.assertEqual(r["flower_area_px"], 0)
        self.assertEqual(r["flower_area_over_plant_area"], 0)
        self.assertIsNone(r["flower_relative_vertical_position"])

    def test_truncated_or_unstandardized_view_restricted(self):
        clipped = Image.new("L", (100, 100), 0)
        ImageDraw.Draw(clipped).rectangle((0, 20, 35, 50), fill=255)
        r = summarize_masks(clipped, whole_plant=True, standardized_view=False)
        self.assertFalse(r["valid_for_cross_image_shape_comparison"])
        self.assertIn("plant_touches_image_border", r["quality_reasons"])
        self.assertIn("view_not_standardized", r["quality_reasons"])

    def test_flower_must_be_within_plant(self):
        bad = Image.new("L", (100, 100), 0)
        bad.putpixel((90, 90), 255)
        with self.assertRaises(ValueError):
            summarize_masks(self.plant, bad)

    def test_dimensions_must_match(self):
        with self.assertRaises(ValueError):
            summarize_masks(self.plant, Image.new("L", (90, 100), 0))

    def test_too_small_or_empty_plant(self):
        with self.assertRaises(ValueError):
            summarize_masks(Image.new("L", (100, 100), 0))
        tiny = Image.new("L", (100, 100), 0)
        tiny.putpixel((50, 50), 255)
        r = summarize_masks(tiny, whole_plant=True, standardized_view=True)
        self.assertIn("target_roi_too_small", r["quality_reasons"])


if __name__ == "__main__":
    unittest.main()
