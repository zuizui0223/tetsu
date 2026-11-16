"""Synthetic N02 conventional rail and N13 independent road-control QA.

Unit tests do NOT contain real ecological observations or railway effects.
"""
import sys
import unittest
from pathlib import Path

from shapely.geometry import LineString
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"scripts"))

from tetsu.ordinary_rail_corridors import (
    is_ordinary_rail,route_normal_pairs,classify_candidate,
    corridor_aggregate,bin_distance_m
)
from audit_japan_ordinary_rail_road_separation import REGIONS,private_scope_includes


class OrdinaryRailCorridorScreens(unittest.TestCase):
    def test_ordinary_jr_and_other_vs_high_speed_and_streettram(self):
        record={
            "N02_001":"11","N02_002":"2",
            "N02_003":"信越本線","N02_004":"東日本旅客鉄道"
        }
        self.assertTrue(is_ordinary_rail(record))
        self.assertTrue(is_ordinary_rail({
            **record,"N02_001":"12","N02_003":"別の在来線"
        }))
        self.assertFalse(is_ordinary_rail({
            **record,"N02_001":"21","N02_003":"宇都宮芳賀ライトレール線"
        }))
        self.assertFalse(is_ordinary_rail({
            **record,"N02_002":"1","N02_003":"東北新幹線"
        }))
        self.assertFalse(is_ordinary_rail({
            **record,"N02_003":"仙台市地下鉄南北線"
        }))
        # Original v1.5 source audit admitted public-operated "東西線":
        # that is Sendai subway, not an exposed railway grass corridor.
        self.assertFalse(is_ordinary_rail({
            **record,"N02_001":"12","N02_002":"3",
            "N02_003":"東西線","N02_004":"仙台市"
        }))

    def test_rail_anchor_geography_from_fixed_lengths_not_plant_phenology(self):
        rail=[LineString([(2000.,10000.),(18000.,10000.)])]
        out,n=route_normal_pairs(
            rail,bbox_m=(0,0,20000,20000),n_cap=30
        )
        self.assertGreaterEqual(n,8)
        self.assertGreaterEqual(len(out),12)
        for near,far in out:
            self.assertAlmostEqual(abs(near.y-10000.),35.,places=4)
            self.assertAlmostEqual(abs(far.y-10000.),1600.,places=4)
            self.assertAlmostEqual(near.x,far.x,places=4)

    def test_distinct_surface_road_control_required(self):
        props=classify_candidate(
            None,None,near_road=125.,far_road=20.,
            nearest_any_rail_far=1700.,near_station=430.
        )
        self.assertTrue(props["strict_near_road75_far_road50_far_rail1000_station300"])
        self.assertTrue(props["near_road_ge100"])
        for variable,value in (
            ("near_road",60.),("far_road",80.),
            ("nearest_any_rail_far",950.),("near_station",290.)
        ):
            args={
                "near_road":125.,"far_road":20.,
                "nearest_any_rail_far":1700.,"near_station":430.
            }
            args[variable]=value
            x=classify_candidate(None,None,**args)
            self.assertFalse(
                x["strict_near_road75_far_road50_far_rail1000_station300"]
            )

    def test_2024_road_near_is_not_biological_rail_exposure(self):
        bad=classify_candidate(
            None,None,near_road=11.,far_road=20.,
            nearest_any_rail_far=1700.,near_station=430.
        )
        self.assertFalse(bad["near_road_ge75"])
        self.assertFalse(bad["strict_near_road75_far_road50_far_rail1000_station300"])
        grouped=corridor_aggregate([
            {"eligibility":bad,"near_road_m":11.,"far_road_m":20.}
        ])
        self.assertEqual(grouped["n_geographically_sampled_rail_sides"],1)
        self.assertTrue(grouped["each_side_is_correlated_with_rail_segment"])
        self.assertEqual(grouped["road_distance_near_bins"]["lt25m"],1)

    def test_no_missing_source_assumed_zero_road(self):
        with self.assertRaises(ValueError):
            classify_candidate(None,None,near_road=float("nan"),
                               far_road=20.,nearest_any_rail_far=1300.,
                               near_station=420.)

    def test_all_fixed_geo_scope_is_posthoc_only_not_main10_screen(self):
        strict=classify_candidate(
            None,None,near_road=125.,far_road=20.,
            nearest_any_rail_far=1700.,near_station=430.
        )
        urban=classify_candidate(
            None,None,near_road=10.,far_road=20.,
            nearest_any_rail_far=1700.,near_station=430.
        )
        self.assertTrue(private_scope_includes(strict,"strict"))
        self.assertFalse(private_scope_includes(urban,"strict"))
        self.assertTrue(private_scope_includes(strict,"all_fixed_geo"))
        self.assertTrue(private_scope_includes(urban,"all_fixed_geo"))
        with self.assertRaises(ValueError):
            private_scope_includes(strict,"cherry_pick_by_NDVI")

    def test_regions_separate_geographies_and_within_first_level_mesh(self):
        self.assertEqual(set(REGIONS),{
            "sendai_miyagi","niigata","nagano"
        })
        for entry in REGIONS.values():
            xmin,ymin,xmax,ymax=entry["bbox"]
            self.assertGreater(xmax,xmin)
            self.assertGreater(ymax,ymin)
            self.assertEqual(entry["utm"],32654)
            self.assertTrue(entry["mesh"].isdigit())


if __name__=="__main__":
    unittest.main()
