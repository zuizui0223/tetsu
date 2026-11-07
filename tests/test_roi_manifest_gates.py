"""End-to-end mask manifest tests guarding against patch-as-individual inference."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from PIL import Image, ImageDraw

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "extract_photo_roi.py"


class RoiManifestSafety(unittest.TestCase):
    def make_run(self, scope, reviewed):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            plant = Image.new("L", (100, 100), 0)
            ImageDraw.Draw(plant).rectangle((20, 20, 65, 75), fill=255)
            flower = Image.new("L", (100, 100), 0)
            ImageDraw.Draw(flower).rectangle((30, 22, 45, 32), fill=255)
            plant.save(root / "plant.png")
            flower.save(root / "flower.png")
            manifest = root / "manifest.csv"
            fields = [
                "image_id", "site_id", "taxon", "photo_date", "photographer_id",
                "site_class", "whole_plant", "standardized_view", "phenophase",
                "plant_mask", "flower_mask", "target_scope", "mask_reviewed"
            ]
            with manifest.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerow({
                    "image_id": "synth01", "site_id": "safe01", "taxon": "Fictitia exemplaris",
                    "photo_date": "2026-06-01", "photographer_id": "researcher",
                    "site_class": "unknown", "whole_plant": "true",
                    "standardized_view": "true", "phenophase": "flowering",
                    "plant_mask": "plant.png", "flower_mask": "flower.png",
                    "target_scope": scope, "mask_reviewed": str(reviewed).lower()
                })
            out = root / "roi.csv"
            run = subprocess.run([
                sys.executable, str(SCRIPT), "--manifest", str(manifest),
                "--root", str(root), "--output", str(out)
            ], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            with out.open(newline="") as f:
                return next(csv.DictReader(f))

    def test_reviewed_single_individual_passes(self):
        row = self.make_run("individual", True)
        self.assertEqual(row["quality_status"], "pass")
        self.assertEqual(row["valid_for_cross_image_shape_comparison"], "True")

    def test_reviewed_patch_restricted(self):
        row = self.make_run("patch", True)
        self.assertEqual(row["quality_status"], "restricted")
        self.assertEqual(row["valid_for_cross_image_shape_comparison"], "False")
        self.assertIn("not_individual_target_scope", row["quality_reasons"])
        self.assertGreater(float(row["flower_area_over_plant_area"]), 0.0)

    def test_unreviewed_individual_restricted(self):
        row = self.make_run("individual", False)
        self.assertEqual(row["quality_status"], "restricted")
        self.assertIn("mask_not_independently_reviewed", row["quality_reasons"])


if __name__ == "__main__":
    unittest.main()
