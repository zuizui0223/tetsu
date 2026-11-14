"""Source-selection tests only: synthetic geometry/rasters are not biological results."""
import math
from pathlib import Path
import sys
import unittest

import numpy as np
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from pyproj import Transformer

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_japan_rail_macro_holdout_worldcover_v17 import (
    HOLDOUT, validate_holdout_frame, read_esa_point, summarize_region,
    holdout_summary
)


class MacroHoldoutTests(unittest.TestCase):
    def test_holdout_mesh_projection_and_worldcover_tiles(self):
        validate_holdout_frame()
        self.assertEqual(len(HOLDOUT),8)
        self.assertNotIn("sendai_miyagi",HOLDOUT)
        self.assertNotIn("niigata",HOLDOUT)
        self.assertNotIn("nagano",HOLDOUT)

    def test_crossing_mesh_or_utm_fails_closed(self):
        bad={"foo":{"mesh":"5339","utm":32654,"bbox":(139.2,35.8,139.6,36.1)}}
        with self.assertRaisesRegex(ValueError,"mesh crossing"):
            validate_holdout_frame(bad)
        bad2={"foo":{"mesh":"5339","utm":32653,"bbox":(139.2,35.4,139.6,35.8)}}
        with self.assertRaisesRegex(ValueError,"projected UTM"):
            validate_holdout_frame(bad2)

    def test_region_specific_pixel_transform_and_class_domain(self):
        pixels=np.full((7,7),30,dtype=np.uint8)
        pixels[3,3]=10
        pixels[3,4]=0
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",height=7,width=7,count=1,dtype="uint8",
                          transform=from_origin(139.,36.,0.0002,0.0002),
                          crs="EPSG:4326",nodata=0) as ds:
                ds.write(pixels,1)
                lon,lat=ds.xy(3,3)
                x,y=Transformer.from_crs(4326,32654,always_xy=True).transform(lon,lat)
                output=read_esa_point(ds,x,y,32654,"N33E138")
                self.assertEqual(output["class"],10)
                self.assertEqual(output["context"]["n_classified"],8)
                with self.assertRaisesRegex(ValueError,"outside prespecified"):
                    read_esa_point(ds,x,y,32654,"N36E138")

    @staticmethod
    def _record(name,nc,rc,road):
        return {
            "route_name":name,"operator_name":"example",
            "near_class":nc,"road_class":rc,
            "near_road_m":road,"far_road_m":20.0,
            "far_rail_m":1800.0,"near_station_m":550.,
        }

    def test_selection_threshold_is_not_species_effect(self):
        rows=[
            self._record("A",10,10,100),
            self._record("B",30,30,60),
            self._record("C",40,40,30),
            self._record("D",30,30,0),
        ]
        report=summarize_region(rows)
        self.assertEqual(report["n_geo_selected_rail_normals"],4)
        self.assertEqual(report["road_control_eligible_by_near_road_min_m"]["75"]["n"],1)
        self.assertEqual(report["road_control_eligible_by_near_road_min_m"]["50"]["n"],2)
        self.assertEqual(report["road_control_eligible_by_near_road_min_m"]["25"]["n"],3)
        self.assertEqual(report["road_control_eligible_by_near_road_min_m"]["75"]["near_tree_minus_all_fixed_fraction"],.75)
        self.assertEqual(report["road_control_eligible_by_near_road_min_m"]["50"]["both_near_and_road_grass_n"],1)
        self.assertTrue(report["not_a_vegetation_or_adaptation_effect"])

    def test_unknown_class_is_not_grass_and_zero_eligibility_is_null(self):
        rows=[self._record("A",None,None,0),
              self._record("B",10,30,0)]
        report=summarize_region(rows)
        self.assertEqual(report["n_classified"],1)
        self.assertEqual(report["near_2021_WorldCover_class_counts"]["unknown"],1)
        self.assertIsNone(report["road_control_eligible_by_near_road_min_m"]["75"]["near_tree_fraction"])
        self.assertIsNone(report["road_control_eligible_by_near_road_min_m"]["75"]["near_tree_minus_all_fixed_fraction"])

    def test_descriptive_aggregation_uses_regions_not_repeated_points(self):
        a={"n_geo_selected_rail_normals":25,
           "road_control_eligible_by_near_road_min_m":{
               k:{"n":3,"near_tree_minus_all_fixed_fraction":v}
               for k,v in zip(("25","50","75","100"),(.1,.2,.3,.4))}}
        b={"n_geo_selected_rail_normals":30,
           "road_control_eligible_by_near_road_min_m":{
               k:{"n":4,"near_tree_minus_all_fixed_fraction":v}
               for k,v in zip(("25","50","75","100"),(-.1,-.2,-.1,-.4))}}
        c={"n_geo_selected_rail_normals":15,
           "road_control_eligible_by_near_road_min_m":{
               k:{"n":7,"near_tree_minus_all_fixed_fraction":.9}
               for k in ("25","50","75","100")}}
        result=holdout_summary({"a":a,"b":b,"c":c})
        self.assertEqual(result["75"]["n_comparable_geographically_separate_regions"],2)
        self.assertEqual(result["75"]["n_regions_with_higher_tree_fraction"],1)
        self.assertEqual(result["75"]["n_regions_with_lower_tree_fraction"],1)
        self.assertAlmostEqual(result["75"]["median_delta_tree_fraction"],.1)
        self.assertNotIn("p_value",result["75"])


if __name__=="__main__":
    unittest.main()
