"""Synthetic historical OSM date/road network QA, not ground truth roads."""
import sys
import unittest
from pathlib import Path
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_lrt_2017_osm_road_history import (
    DATES,OSM_QUERY_FILTER,date_aoi,lines_from_osm,nearest_metric
)


class HistoricalOSMRoadQA(unittest.TestCase):
    def test_pre_2018_date_and_fixed_2024_reference(self):
        self.assertEqual(DATES["2017"],"2017-12-31")
        self.assertEqual(DATES["2024"],"2024-09-30")
        self.assertIn("highway in (",OSM_QUERY_FILTER)
        self.assertIn("type:way",OSM_QUERY_FILTER)

    def test_bbox_order_lon_lat(self):
        ll=date_aoi(140.,36.5)
        self.assertEqual(ll,(139.972,36.48,140.028,36.52))

    def test_empty_OSM_is_unknown_physical_road_not_zero_distance(self):
        out=nearest_metric({"n":(400000.,4000000.)},[])
        self.assertIsNone(out["n"])

    def test_actual_2017_line_road_distances_not_invented(self):
        rows=[
            {
                "type":"Feature",
                "geometry":{"type":"LineString","coordinates":[
                    [139.999,36.500],[140.001,36.500]
                ]},
                "properties":{}
            },
            {
                "type":"Feature",
                "geometry":{"type":"Point","coordinates":[140.,36.501]},
                "properties":{}
            }
        ]
        utm=Transformer.from_crs(4326,32654,always_xy=True)
        lines,counts=lines_from_osm(rows,utm)
        self.assertEqual(len(lines),1)
        self.assertEqual(counts["non_road_line_geometry"],1)
        near=utm.transform(140.,36.500)
        far=utm.transform(140.,36.510)
        d=nearest_metric({"near":near,"far":far},lines)
        self.assertLess(d["near"],5)
        self.assertGreater(d["far"],900)


if __name__=="__main__":
    unittest.main()
