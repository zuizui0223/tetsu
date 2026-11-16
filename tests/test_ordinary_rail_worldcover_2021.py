"""Synthetic ESA WorldCover source/CRS/class tests, not field vegetation."""
import sys
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin,xy
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"scripts"))

from tetsu.worldcover_landcover import (
    CLASS_LABELS,tile_name,tile_url,normalize_class,
    summarize_neighborhood,summarize_route_source
)
from audit_japan_ordinary_rail_worldcover_2021 import source_sample


class WorldCoverRailGroundVersusCanopyTests(unittest.TestCase):
    def test_real_3_by_3_degree_2021_Japan_tile(self):
        self.assertEqual(tile_name(140.6,38.3),"N36E138")
        self.assertEqual(tile_name(139.3,37.7),"N36E138")
        self.assertEqual(tile_name(138.2,36.4),"N36E138")
        self.assertEqual(tile_name(138.0,36.0),"N36E138")
        self.assertEqual(tile_name(137.999,35.999),"N33E135")
        self.assertEqual(tile_name(-0.1,-0.1),"S03W003")
        self.assertIn("2021_v200_N36E138_Map.tif",tile_url("N36E138"))
        with self.assertRaises(ValueError):
            tile_url("N36E138",base="http://wrong.example")

    def test_tree_is_not_wild_herbaceous_and_0_is_unknown(self):
        self.assertEqual(CLASS_LABELS[10],"tree_cover")
        self.assertEqual(CLASS_LABELS[30],"grassland")
        self.assertEqual(CLASS_LABELS[40],"cropland")
        self.assertIsNone(normalize_class(0))
        with self.assertRaises(ValueError):
            normalize_class(255)

    def test_window_tie_is_unknown_not_forced_trees_or_grassland(self):
        p=summarize_neighborhood([10,10,10,30,30,30,40,40,40])
        self.assertTrue(p["plurality_tied"])
        self.assertIsNone(p["plurality_class"])
        self.assertEqual(p["n_classified"],9)
        self.assertAlmostEqual(p["tree_cover_fraction_among_classified"],3/9,places=4)
        self.assertEqual(summarize_neighborhood([0]*9)["n_classified"],0)

    def test_route_aggregate_keeps_crop_tree_built_unknown_separate(self):
        r=summarize_route_source([
            {"near_class":10,"control_class":30,
             "near_window":summarize_neighborhood([10]*9)},
            {"near_class":30,"control_class":30,
             "near_window":summarize_neighborhood([30]*9)},
            {"near_class":40,"control_class":50,
             "near_window":summarize_neighborhood([40]*9)},
            {"near_class":None,"control_class":10,
             "near_window":summarize_neighborhood([0]*9)}
        ])
        self.assertEqual(r["n_rail_road_geometry_selected_pairs"],4)
        self.assertEqual(r["n_near_2021_worldcover_tree_cover"],1)
        self.assertEqual(r["n_near_2021_worldcover_grassland_or_herbaceous_wetland"],1)
        self.assertEqual(r["n_near_and_control_potential_ground_herbaceous"],1)
        self.assertEqual(r["near_2021_worldcover_class_counts"]["unknown"],1)
        self.assertEqual(r["near_2021_worldcover_class_counts"]["cropland"],1)

    def test_actual_rasterio_pixel_read_uses_utm54_to_geo4326(self):
        # Tiny in-memory 7x7 ESA-looking map, center is class 30;
        # all other pixels tree 10. The inputs are SYNTHETIC.
        image=np.full((7,7),10,dtype="uint8")
        image[3,3]=30
        trans=from_origin(140.62,38.35,1/12000,1/12000)
        lon,lat=xy(trans,3,3,offset="center")
        e,n=Transformer.from_crs(4326,32654,always_xy=True).transform(lon,lat)
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=7,height=7,
                          count=1,dtype="uint8",nodata=0,
                          crs="EPSG:4326",transform=trans) as source:
                source.write(image,1)
                out=source_sample(source,e,n)
        self.assertEqual(out["center_class"],30)
        self.assertEqual(out["window"]["n_classified"],9)
        self.assertEqual(out["window"]["plurality_class"],10)
        self.assertAlmostEqual(out["window"]["tree_cover_fraction_among_classified"],8/9,places=4)
        self.assertFalse(out["window"]["plurality_tied"])

    def test_wrong_source_CRS_fails_closed(self):
        arr=np.zeros((7,7),dtype="uint8")
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=7,height=7,
                          count=1,dtype="uint8",nodata=0,
                          crs="EPSG:3857",transform=from_origin(
                              140.62,38.35,1/12000,1/12000
                          )) as source:
                source.write(arr,1)
                with self.assertRaises(ValueError):
                    source_sample(source,450000,4000000)


if __name__=="__main__":
    unittest.main()
