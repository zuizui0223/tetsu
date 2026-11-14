"""Synthetic POSTHOC source-frame audit tests; never biological railway effects."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"scripts"))

from audit_japan_ordinary_rail_worldcover_expanded import (
    source_comparator_eligibility,audit_region
)


class ExpandedFixedGeographyWorldCoverTests(unittest.TestCase):
    def test_near_road_thresholds_still_need_independent_road_control(self):
        sample={
            "near_road_m":63.,
            "far_road_m":25.,
            "far_rail_m":1500.,
            "near_station_m":500.
        }
        self.assertTrue(source_comparator_eligibility(sample,25))
        self.assertTrue(source_comparator_eligibility(sample,50))
        self.assertFalse(source_comparator_eligibility(sample,75))
        with self.assertRaises(ValueError):
            source_comparator_eligibility(sample,0)
        for key,val in (
            ("far_road_m",65.),("far_rail_m",950.),
            ("near_station_m",250.)
        ):
            invalid=dict(sample)
            invalid[key]=val
            self.assertFalse(source_comparator_eligibility(invalid,25))

    def test_fixed_original_81_source_sites_and_only2_strict_preserved(self):
        sites=[]
        for i in range(81):
            strict=i<2
            sites.append({
                "route_name":"合成仮想線",
                "operator_name":"合成運営者",
                "near_xy":(float(i),1.),
                "road_only_xy":(float(i),2.),
                "near_road_m":85. if strict else 20.,
                "far_road_m":15.,
                "far_rail_m":1400.,
                "near_station_m":400.,
                "source_gates":{
                    "strict_near_road75_far_road50_far_rail1000_station300":strict
                }
            })
        def fake2021(_raster,x,y):
            cl=30 if x==0 else 10
            if y==2:cl=10
            return {
                "center_class":cl,
                "window":{
                    "plurality_class":cl,
                    "tree_cover_fraction_among_classified":1. if cl==10 else 0.
                }
            }
        with patch(
            "audit_japan_ordinary_rail_worldcover_expanded.source_sample",
            side_effect=fake2021
        ):
            report=audit_region("sendai_miyagi",sites,object())
        self.assertEqual(report["n_original_GIS_fixed_source_pairs"],81)
        self.assertEqual(report["n_strict_original_road75_source_pairs"],2)
        self.assertEqual(report["all_fixed_source_grassland_or_herb_wetland_n"],1)
        self.assertEqual(
            report["2024_road_distance_threshold_sensitivities"]["75"]["n_2024_road_control_eligible_pairs"],2
        )
        self.assertEqual(
            report["2024_road_distance_threshold_sensitivities"]["25"]["n_2024_road_control_eligible_pairs"],2
        )
        self.assertTrue(report["posthoc_discovery_NOT_confirmatory"])

    def test_outcome_independent_source_cohort_drift_rejected(self):
        with self.assertRaises(ValueError):
            audit_region("sendai_miyagi",[],object())


if __name__=="__main__":
    unittest.main()
