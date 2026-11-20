"""v307 official original map source gate: synthetic PDFs ONLY for unit tests."""
import sys
import tempfile
import unittest
from pathlib import Path
import fitz
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Czech_2022_2023_official_actual_train_traffic_maps_v307 import inspect, verify_pair

class TestOriginalMapSourceNoFabricatedExposure(unittest.TestCase):
    def _create(self, p, message):
        d=fitz.open()
        page=d.new_page(width=450,height=600)
        page.insert_text((40,80),message*400,fontname="helv",fontsize=11)
        # An authentic synthetic PDF source must exceed 10k bytes; add a small
        # randomly colored procedural image without exporting actual copyrighted map.
        import zlib
        import random
        rng=random.Random(337)
        arr=bytes(rng.randrange(256) for _ in range(150*150*3))
        pix=fitz.Pixmap(fitz.csRGB,150,150,arr,False)
        page.insert_image(fitz.Rect(0,150,420,570),pixmap=pix)
        d.save(str(p))
        d.close()
    def test_two_years_source_and_units(self):
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/"map22.pdf";b=Path(d)/"map23.pdf"
            self._create(a,"2022 observed trains Q1")
            self._create(b,"2023 observed trains Q1")
            x=inspect(a,"2022")
            y=inspect(b,"2023")
            z=verify_pair({"2022":x,"2023":y})
            self.assertEqual(z["status"],"PUBLIC_ORIGINAL_CZECH_MAPS_2022_AND_2023_SOURCE_FILES_AVAILABLE")
            self.assertEqual(z["n_source_botanical_station_exposures_linked"],0)
            self.assertFalse(x["freight_passenger_split_proven"])
            self.assertIn("1 January",x["published_months"])
    def test_not_pdf_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"not.pdf";p.write_bytes(b"<!DOCTYPE html>"+b"X"*22000)
            with self.assertRaises(ValueError):
                inspect(p,"2022")
    def test_duplicate_year_pdf_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"same.pdf"
            self._create(p,"2022")
            m=inspect(p,"2022")
            n=inspect(p,"2023")
            with self.assertRaises(ValueError):
                verify_pair({"2022":m,"2023":n})
    def test_missing_2023_map_fails_closed(self):
        with self.assertRaises(ValueError):
            verify_pair({"2022":{}})
if __name__=="__main__": unittest.main()
