"""Synthetic EPSG53/54/55 interoperability tests for real Japanese COG reads.

A same-ground-coordinate observation must be reprojected when satellite
tiles span adjacent UTM zones; neither plants nor numeric measurements
in this fixture are real ecological data.
"""
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock,patch

from pyproj import Transformer
from rasterio.crs import CRS

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from measure_utsunomiya_s2_ndvi_pixels import pixel_samples


class NeighboringSentinelUTMZoneTests(unittest.TestCase):
    def test_zone53_tile_must_reproject_from_n02_zone54(self):
        true_lon,true_lat=138.4,36.5
        input_e,input_n=Transformer.from_crs(
            4326,32654,always_xy=True
        ).transform(true_lon,true_lat)
        expected_e,expected_n=Transformer.from_crs(
            4326,32653,always_xy=True
        ).transform(true_lon,true_lat)
        captured=[]
        mock=MagicMock()
        mock.count=1
        mock.width=10980
        mock.crs=CRS.from_epsg(32653)
        def capture(coords,indexes=1,masked=True):
            captured.extend(coords)
            return [[1800] for _ in coords]
        mock.sample.side_effect=capture
        mock.__enter__.return_value=mock
        with patch("measure_utsunomiya_s2_ndvi_pixels.rasterio.open",
                   return_value=mock):
            result=pixel_samples(
                {"href":"https://sentinel-cogs.s3.amazonaws.com/sample.tif"},
                {"x":(input_e,input_n)}
            )
        self.assertEqual(result["x"],1800.)
        self.assertAlmostEqual(captured[0][0],expected_e,places=3)
        self.assertAlmostEqual(captured[0][1],expected_n,places=3)
        self.assertGreater(abs(captured[0][0]-input_e),100000)

    def test_zone54_tile_keeps_projected_coordinate(self):
        coord=(440000.,4038000.)
        captured=[]
        mock=MagicMock()
        mock.count=1
        mock.width=10980
        mock.crs=CRS.from_epsg(32654)
        def capture(coords,indexes=1,masked=True):
            captured.extend(coords)
            return [[2000] for _ in coords]
        mock.sample.side_effect=capture
        mock.__enter__.return_value=mock
        with patch("measure_utsunomiya_s2_ndvi_pixels.rasterio.open",
                   return_value=mock):
            pixel_samples(
                {"href":"https://sentinel-cogs.s3.amazonaws.com/sample.tif"},
                {"x":coord}
            )
        self.assertEqual(captured,[coord])

    def test_non_utm_Japan_raster_is_rejected(self):
        mock=MagicMock()
        mock.count=1
        mock.width=10980
        mock.crs=CRS.from_epsg(3857)
        mock.__enter__.return_value=mock
        with patch("measure_utsunomiya_s2_ndvi_pixels.rasterio.open",
                   return_value=mock):
            with self.assertRaises(ValueError):
                pixel_samples(
                    {"href":"https://sentinel-cogs.s3.amazonaws.com/bad.tif"},
                    {"x":(450000.,4000000.)}
                )


if __name__=="__main__":
    unittest.main()
