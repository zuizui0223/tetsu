"""N13 road proximity vs 2017 SCL pre-balance unit tests; SYNTHETIC only."""
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import shapefile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"scripts"))
from audit_lrt_n13_road_land_balance import only_aggregates
from tetsu.rail_road_landbalance import (
    projected_road_segments_from_zip,pair_balance_flags,count_balance,
    blocked_summary,verify_url
)


class RoadAndBaselineGates(unittest.TestCase):
    def test_official_url_allowlist(self):
        verify_url("https://nlftp.mlit.go.jp/ksj/gml/data/N13/N13-24/N13-24_5439_SHP.zip")
        with self.assertRaises(ValueError):
            verify_url("https://attacker.example/N13-24_5439_SHP.zip")

    def test_geo_shape_fixture_clipped_and_projected(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            base=folder/"roads"
            with shapefile.Writer(str(base),shapeType=shapefile.POLYLINE) as w:
                w.field("name","C")
                w.line([[(139.9,36.51),(139.91,36.51)]])
                w.record("inside")
                w.line([[(132.,31.),(132.1,31.1)]])
                w.record("outside")
            archive=folder/"test.zip"
            with zipfile.ZipFile(archive,"w") as z:
                z.write(str(base)+".shp",arcname="roads.shp")
                z.write(str(base)+".shx",arcname="roads.shx")
            geoms,meta=projected_road_segments_from_zip(archive)
            self.assertEqual(meta["all_source_lines"],2)
            self.assertEqual(meta["lines_in_source_bbox"],1)
            self.assertEqual(meta["projected_parts_in_bbox"],1)
            self.assertGreater(geoms[0].length,100)

    def test_post_road_matching_not_implied_pre_road_existence(self):
        f=pair_balance_flags(25.,35.,4,4,.51,.57)
        self.assertTrue(f["road_abs_diff_le50_both_le100"])
        self.assertTrue(f["2017_same_scl_class"])
        self.assertTrue(f["2017_both_vegetation_scl4"])
        self.assertTrue(f["road50_sameclass_preNDVI_similar"])
        self.assertEqual(
            count_balance({("n","f"):f})["n_pairs_by_gate"]["road50_sameclass_preNDVI_similar"],1
        )

    def test_2017_veg_vs_bare_cannot_be_landmatched(self):
        f=pair_balance_flags(35.,40.,4,5,.2,.2)
        self.assertTrue(f["road_abs_diff_le50_both_le100"])
        self.assertFalse(f["2017_same_scl_class"])
        self.assertFalse(f["road50_and_2017_same_class"])

    def test_road_distances_from_different_road_access_not_comparable(self):
        f=pair_balance_flags(15.,250.,4,4,.2,.2)
        self.assertFalse(f["near_road_le100_and_far_road_le100"])
        self.assertFalse(f["road_abs_diff_le50_both_le100"])

    def test_real_source_baseline_is_flat_scalar_ndvi_not_nested_dict(self):
        pairs=[("T00_L_near","T00_L_control")]
        pts={"T00_L_near":(0,0),"T00_L_control":(200,300)}
        road={"T00_L_near":14.0,"T00_L_control":30.0}
        pre={"T00_L_near":.3,"T00_L_control":.38}
        source_day={
            "2017-07-22":{
                "T00_L_near":{
                    "provider_flag_branch":.3,"scl_class":4
                },
                "T00_L_control":{
                    "provider_flag_branch":.38,"scl_class":4
                }
            }
        }
        summary,flags=only_aggregates(pts,pairs,road,pre,source_day)
        self.assertEqual(
            summary["2017_abs_ndvi_between_exposed_and_far_pairs_n"],1
        )
        self.assertTrue(
            flags[pairs[0]]["road50_sameclass_preNDVI_similar"]
        )

    def test_diff_in_diff_skips_missing_and_splits_fixed_blocks(self):
        pairs=[("T00_L_near","T00_L_control"),("T20_R_near","T20_R_control")]
        flags={p:{"road50_and_2017_same_class":(i==0)} for i,p in enumerate(pairs)}
        d={
            2017:{
                "T00_L_near":.30,"T00_L_control":.40,
                "T20_R_near":.30,"T20_R_control":.40,
            },
            2024:{
                "T00_L_near":.50,"T00_L_control":.45,
                "T20_R_near":.20,"T20_R_control":.43,
            }
        }
        result=blocked_summary(pairs,d,flags,2024,gate="road50_and_2017_same_class")
        self.assertEqual(result["paired_complete"]["n"],1)
        self.assertAlmostEqual(result["paired_complete"]["median"],.15)
        self.assertEqual(result["n_longitudinal_blocks_with_data"],1)
        result_all=blocked_summary(pairs,d,flags,2024,gate="all")
        self.assertEqual(result_all["n_longitudinal_blocks_with_data"],2)
        self.assertEqual(result_all["paired_complete"]["n"],2)


if __name__=="__main__":
    unittest.main()
