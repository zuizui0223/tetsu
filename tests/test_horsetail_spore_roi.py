"""Synthetic validated horsetail strobilus-mask tests, not ecological data."""
from pathlib import Path
import sys
import tempfile
import unittest
from PIL import Image, ImageDraw

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from extract_horsetail_spore_roi import summarize_row


class HorsetailSporeRoi(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        plant=Image.new("L",(100,100),0)
        ImageDraw.Draw(plant).rectangle((20,20,60,80),fill=255)
        plant.save(self.root/"plant.png")
        cone=Image.new("L",(100,100),0)
        ImageDraw.Draw(cone).rectangle((25,20,35,30),fill=255)
        cone.save(self.root/"cone.png")
        self.row={
            "image_id":"synthetic", "site_id":"site1",
            "taxon":"Equisetum arvense", "photo_date":"2026-04-11",
            "observer_id":"observer", "plant_mask":"plant.png",
            "spore_cone_mask":"cone.png", "mask_reviewed":"true",
            "target_scope":"individual", "whole_plant":"true",
            "standardized_view":"true", "fixed_plot_view":"false",
            "same_taxon_patch":"false", "plot_id":"",
            "spore_phenophase":"fertile_cone_visible",
        }

    def test_reproductive_labels_are_spore_not_flower(self):
        out=summarize_row(self.row,self.root)
        self.assertEqual(out["quality_status"],"pass_individual_shape")
        self.assertEqual(out["reproductive_structure"],"spore_strobilus")
        self.assertGreater(out["spore_cone_area_px"],0)
        self.assertTrue(out["not_a_flowering_plant"])
        self.assertNotIn("flower_area_px",out)
        self.assertNotIn("flower_area_over_plant_area",out)

    def test_fixed_patch_is_spore_fraction_not_individual_size(self):
        row=dict(self.row,target_scope="patch",whole_plant="false",
                 fixed_plot_view="true",same_taxon_patch="true",plot_id="P1")
        out=summarize_row(row,self.root)
        self.assertEqual(out["quality_status"],"pass_fixed_patch_spore_fraction")
        self.assertFalse(out["valid_for_cross_image_shape_comparison"])
        self.assertTrue(out["valid_for_fixed_plot_spore_cone_fraction"])

    def test_no_annotation_cannot_become_cone_absence(self):
        row=dict(self.row,spore_cone_mask="")
        out=summarize_row(row,self.root)
        self.assertEqual(out["quality_status"],"rejected")
        self.assertIn("requires independently annotated cone mask",out["quality_error"])

    def test_no_flower_plant_restriction(self):
        row=dict(self.row,taxon="Solidago altissima")
        out=summarize_row(row,self.root)
        self.assertEqual(out["quality_status"],"rejected")
        self.assertIn("only accepts Equisetum",out["quality_error"])

    def test_cone_outside_plant_rejected(self):
        bad=Image.new("L",(100,100),0)
        bad.putpixel((90,90),255)
        bad.save(self.root/"bad.png")
        row=dict(self.row,spore_cone_mask="bad.png")
        out=summarize_row(row,self.root)
        self.assertEqual(out["quality_status"],"rejected")
        self.assertIn("outside plant ROI",out["quality_error"])

    def test_no_cone_at_phase_marked_absent(self):
        empty=Image.new("L",(100,100),0)
        empty.save(self.root/"empty.png")
        row=dict(self.row,spore_cone_mask="empty.png",
                 spore_phenophase="fertile_cone_absent_visible")
        out=summarize_row(row,self.root)
        self.assertEqual(out["quality_status"],"pass_individual_shape")
        self.assertEqual(out["spore_cone_area_px"],0)
        self.assertEqual(out["spore_cone_area_over_plant_area"],0)


if __name__=="__main__":
    unittest.main()
