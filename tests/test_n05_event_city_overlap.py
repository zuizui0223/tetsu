"""Unit tests of source geometry overlap and sample-cohort classification."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_n05_event_city_overlap import bbox_overlap,summarize_n05_city_event_support


class FakeShape:
    shapeType=3
    def __init__(self,bounds):self.bbox=bounds


class FakeFeature:
    def __init__(self,row,bounds):
        self.record=row
        self.shape=FakeShape(bounds)


class FakeReader:
    fields=[("DeletionFlag","C",1,0),
            ("N05_004","C",10,0),
            ("N05_005b","C",10,0),
            ("N05_005e","C",10,0),
            ("N05_006","C",15,0)]
    def __init__(self,features):self.features=features
    def iterShapeRecords(self):return iter(self.features)


class EventCityAuditTests(unittest.TestCase):
    AREAS={"sample":(35.0,139.0,35.5,139.5)}

    def test_lonlat_bbox_intersection(self):
        self.assertTrue(bbox_overlap(
            (139.1,35.2,139.15,35.3), self.AREAS["sample"],padding=0
        ))
        self.assertFalse(bbox_overlap(
            (140.2,35.2,140.4,35.3), self.AREAS["sample"],padding=0
        ))

    def test_open_and_potential_closed_groups(self):
        features=[
            FakeFeature(["2020","2020","9999","new1"],
                        (139.1,35.2,139.2,35.3)),
            FakeFeature(["2000","2000","2019","active_changed"],
                        (139.1,35.2,139.2,35.3)),
            FakeFeature(["2000","2020","9999","active_changed"],
                        (139.1,35.2,139.2,35.3)),
            FakeFeature(["2011","2011","2018","closed_candidate"],
                        (139.1,35.2,139.2,35.3)),
            FakeFeature(["2012","2012","9999","not_in_area"],
                        (140.1,35.2,140.2,35.3))
        ]
        out=summarize_n05_city_event_support(FakeReader(features),self.AREAS)
        self.assertEqual(out["n_relation_groups"],4)
        c=out["city_event_candidate_groups"]["sample"]
        self.assertEqual(c["opening_known_2010_2025"],2)
        self.assertEqual(c["possible_terminal_not_yet_verified_2010_2024"],1)
        self.assertEqual(c["possible_terminal_year_distribution"],{"2018":1})


if __name__=="__main__":
    unittest.main()
