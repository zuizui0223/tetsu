"""Unit tests for Czech original PDF rotated-text glyph source-only audit."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
    import fitz
    from audit_Czech_original_station_to_traffic_map_glyph_identifiability_v308 import (
        glyph_pairs_for_station,number_glyphs,normalize_source_map_name,
        ORIGINAL_CZECH_BOTANICAL_STATIONS,
    )
except ModuleNotFoundError as exc:
    if exc.name not in ("fitz","openpyxl"):
        raise
    fitz=None

@unittest.skipUnless(fitz is not None,"PyMuPDF/openpyxl installed only in the dedicated first-party source CI")
class GenuineSourceGlyphOnlyTests(unittest.TestCase):
    def test_original_pub_station_list_has_39_nonduplicate_normalized_labels(self):
        self.assertEqual(len(ORIGINAL_CZECH_BOTANICAL_STATIONS),39)
        self.assertEqual(len(set(map(normalize_source_map_name,ORIGINAL_CZECH_BOTANICAL_STATIONS))),39)
    def test_90_degree_pdf_rotation_cannot_corrupt_rendered_coordinate(self):
        d=fitz.open()
        p=d.new_page(width=595,height=842)
        p.insert_text((100,200),"Lysa nad Labem",fontsize=11)
        p.insert_text((110,209),"1",fontsize=11)
        p.set_rotation(90)
        labels=glyph_pairs_for_station(p,"Lysa nad Labem",number_glyphs(p),radius_pt=25)
        self.assertGreaterEqual(len(labels),1)
        q=labels[0]
        b=q["label_bbox_rotated_points"]
        self.assertGreaterEqual(b[0],0)
        self.assertGreaterEqual(b[1],0)
        self.assertLess(b[2],p.rect.width+1)
        self.assertLess(b[3],p.rect.height+1)
        self.assertEqual(q["nearby_numeric_glyphs_UNCLASSIFIED"][0]["n"],1)
        d.close()
    def test_number_is_only_glyph_not_proven_train_traffic(self):
        d=fitz.open()
        p=d.new_page(width=595,height=842)
        p.insert_text((20,160),"Colin",fontsize=11)
        p.insert_text((22,180),"135",fontsize=11)
        q=glyph_pairs_for_station(p,"Colin",number_glyphs(p),radius_pt=25)
        self.assertGreaterEqual(len(q),1)
        self.assertNotIn("train_count",str(q))
        self.assertNotIn("traffic",str(q))
        d.close()
if __name__=="__main__":unittest.main()
