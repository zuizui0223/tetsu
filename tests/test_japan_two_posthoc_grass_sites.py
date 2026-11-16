"""Synthetic posthoc candidate check and GSI source-availability fixtures."""
import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch,MagicMock
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_japan_two_posthoc_grass_sites_S2_GSI import (
    REGION_CANDIDATES,eligible_posthoc_worldcover_grass,
    geo_tile,gsi_source_tile_check
)


def site(near_road=60.):
    return {
        "near_road_m":near_road,
        "far_road_m":12.,
        "far_rail_m":1650.,
        "near_station_m":500.,
        "near_xy":(425000.,4050000.),
        "road_only_xy":(425000.,4051600.),
    }


class TwoPosthocMappedGrassSitesQA(unittest.TestCase):
    def test_original_verified_route_labels_are_frozen(self):
        self.assertEqual(set(REGION_CANDIDATES),{"niigata","nagano"})
        self.assertEqual(REGION_CANDIDATES["niigata"][0],"信越線")
        self.assertEqual(REGION_CANDIDATES["nagano"][0],"小海線")

    def test_2021_grass_and_2024_road50_required_before_2025_S2(self):
        with patch(
            "audit_japan_two_posthoc_grass_sites_S2_GSI.source_sample",
            return_value={"center_class":30}
        ):
            self.assertTrue(eligible_posthoc_worldcover_grass(site(),object()))
            self.assertFalse(eligible_posthoc_worldcover_grass(site(near_road=49),object()))
        with patch(
            "audit_japan_two_posthoc_grass_sites_S2_GSI.source_sample",
            return_value={"center_class":10}
        ):
            self.assertFalse(eligible_posthoc_worldcover_grass(site(),object()))

    def test_raster_GSI_tiles_use_valid_Japanese_ZXY(self):
        for z in (17,18,19):
            zz,x,y=geo_tile(139.12,37.7,z)
            self.assertEqual(zz,z)
            self.assertTrue(0<=x<2**z)
            self.assertTrue(0<=y<2**z)
        with self.assertRaises(ValueError):
            geo_tile(2.,1.,17)

    def test_GSI_availability_does_not_assert_at_grade_or_grass(self):
        im=Image.new("RGB",(256,256))
        for y in range(256):
            for x in range(256):
                im.putpixel((x,y),(x%256,y%256,(x+y)%256))
        data=io.BytesIO()
        im.save(data,format="JPEG",quality=85)
        response=MagicMock()
        response.status_code=200
        response.content=data.getvalue()
        response.raise_for_status.return_value=None
        client=MagicMock()
        client.get.return_value=response
        status=gsi_source_tile_check(client,139.14,37.72)
        self.assertEqual(len(status),3)
        self.assertTrue(all(v=="nonblank_current_aerial_source"
                            for v in status.values()))
        self.assertEqual(client.get.call_count,3)


if __name__=="__main__":
    unittest.main()
