"""Pure tests of rail-road proximity controls and descriptive odds calculators."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.transport_overlap import (road_rail_status,
                                     weighted_case_control_table,
                                     monthly_mh_ratio)


def obs(label,rail,road,w=1,month=6,block="a"):
    return {"label":label,"rail_m":rail,"road_m":road,
            "ref_weight":w,"month":month,"grid2km":block}


class TransportOverlapTests(unittest.TestCase):
    def test_road_matched_bands(self):
        self.assertEqual(road_rail_status(250,100),"road_matched_rail_near")
        self.assertEqual(road_rail_status(251,100),"road_matched_rail_middle_excluded")
        self.assertEqual(road_rail_status(1000,100),"road_matched_rail_far")
        self.assertEqual(road_rail_status(5,101),"road_distant")

    def test_invalid_geometry_rejected(self):
        with self.assertRaises(ValueError):
            road_rail_status(float("nan"),20)
        with self.assertRaises(ValueError):
            road_rail_status(0,-3)

    def test_basic_descriptive_ratio(self):
        rows=[
            obs("target",50,25),obs("target",80,50),
            obs("target",1600,30),obs("reference",75,50),
            obs("reference",1250,80),obs("reference",1550,35)
        ]
        m=weighted_case_control_table(rows)
        self.assertEqual(m["target_rail_near"],2)
        self.assertEqual(m["target_rail_far"],1)
        self.assertEqual(m["reference_weighted_rail_near"],1.0)
        self.assertEqual(m["reference_weighted_rail_far"],2.0)
        self.assertGreater(m["approx_target_to_reference_ratio_of_ratios"],1)
        self.assertEqual(m["n_grid2km_with_valid_comparison_rows"],1)

    def test_sampled_reference_weights(self):
        rows=[obs("target",30,10),obs("target",1400,10),
              obs("reference",50,10,10),obs("reference",1500,10,5)]
        m=weighted_case_control_table(rows)
        self.assertEqual(m["reference_weighted_rail_near"],10)
        self.assertEqual(m["reference_weighted_rail_far"],5)

    def test_monthwise_mh_ratio(self):
        rows=[obs("target",10,10,month=5),obs("target",1400,10,month=5),
              obs("reference",10,10,month=5),
              obs("reference",1400,10,month=5)]
        r=monthly_mh_ratio(rows)
        self.assertEqual(r["seasonal_mantel_haenszel_ratio"],1.0)
        self.assertEqual(r["months_with_both_target_and_reference"],1)


if __name__=="__main__":
    unittest.main()
