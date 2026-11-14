"""Tests for geolocated GSI annual orthophoto availability only."""
import io
import sys
import unittest
from pathlib import Path

from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_lrt_gsi_2017_aerial_tile_coverage import (
    world_tile,image_is_nonblank
)


class GsiAerialPhotoSourceTests(unittest.TestCase):
    def test_valid_japan_2017_tile_coordinate(self):
        z,x,y=world_tile(139.92,36.56,17)
        self.assertEqual(z,17)
        self.assertTrue(0<=x<2**z)
        self.assertTrue(0<=y<2**z)
        with self.assertRaises(ValueError):
            world_tile(1.0,0,17)

    def test_bright_blank_tile_not_identified_as_aerial_photo(self):
        im=Image.new("RGBA",(256,256),(255,255,255,255))
        raw=io.BytesIO()
        im.save(raw,format="PNG")
        self.assertEqual(
            image_is_nonblank(raw.getvalue()),"solid_color_not_orthophoto"
        )

    def test_transparent_tile_not_identified_as_historical_coverage(self):
        im=Image.new("RGBA",(256,256),(0,0,0,0))
        raw=io.BytesIO()
        im.save(raw,format="PNG")
        self.assertEqual(
            image_is_nonblank(raw.getvalue()),"mostly_transparent_or_missing"
        )

    def test_modern_seamless_JPEG_positive_control_is_not_historical_photo(self):
        im=Image.new("RGB",(256,256))
        for y in range(256):
            for x in range(256):
                im.putpixel((x,y),((2*x)%256,(3*y)%256,(x+y)%256))
        raw=io.BytesIO()
        im.save(raw,format="JPEG",quality=88)
        self.assertEqual(
            image_is_nonblank(raw.getvalue(),expected_format="JPEG"),
            "photo_pixels_present"
        )
        with self.assertRaises(ValueError):
            image_is_nonblank(raw.getvalue()) # annual photos must be PNG

    def test_variable_tile_pixels_pass_source_photo_gate(self):
        im=Image.new("RGB",(256,256))
        for y in range(256):
            for x in range(256):
                im.putpixel((x,y),(x%256,y%256,(x+y)%256))
        raw=io.BytesIO()
        im.save(raw,format="PNG")
        self.assertEqual(image_is_nonblank(raw.getvalue()),"photo_pixels_present")


if __name__=="__main__":
    unittest.main()
