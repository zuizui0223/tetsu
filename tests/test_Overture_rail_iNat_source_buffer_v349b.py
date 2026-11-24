#!/usr/bin/env python3
"""Synthetic original full Overture source buffer/city membership regression."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from pilot_Overture_rail_iNat_six_continent_v349b import (
    source_rail_query_bbox,SOURCE_RAIL_BUFFER_M,CITIES
)
from pilot_self_extracted_six_continent_rail_flora_v348 import bounds

class V349bSourceWindowTests(unittest.TestCase):
    def test_full_rail_source_geometry_extends_beyond_observation_window(self):
        for city in CITIES:
            name,continent,country,lat,lon=city
            northsouth=111320.0
            a,b,c,d=bounds(lat,lon)
            west,south,east,north=source_rail_query_bbox(lat,lon)
            self.assertGreaterEqual((a-south)*northsouth,SOURCE_RAIL_BUFFER_M-.01)
            self.assertGreaterEqual((north-c)*northsouth,SOURCE_RAIL_BUFFER_M-.01)
            ew=111320.0*math.cos(math.radians(lat))
            self.assertGreaterEqual((b-west)*ew,SOURCE_RAIL_BUFFER_M-.01)
            self.assertGreaterEqual((east-d)*ew,SOURCE_RAIL_BUFFER_M-.01)
    def test_population_unchanged_six_cities(self):
        self.assertEqual(len(CITIES),6)
        self.assertEqual(len({x[1] for x in CITIES}),6)

if __name__=="__main__":unittest.main()
