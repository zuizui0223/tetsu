#!/usr/bin/env python3
"""v363: exact source connector identity, deterministic unbiased 1/256 UUID selection."""
import sys,hashlib,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_global_GERS_connector_position_consistency_v363 import (
    pick_uuid,MODULUS,EXPECTED_SOURCE_ALL_RAIL,THRESHOLDS,true_WGS84_geodesic_at
)
from shapely.geometry import LineString
from pyproj import Geod
class V363GERSPhysicalLocationQATests(unittest.TestCase):
    def test_source_immutable_sample_probability_over_uuid_IDs(self):
        self.assertEqual(MODULUS,256)
        self.assertEqual(EXPECTED_SOURCE_ALL_RAIL,2054855)
        selection={i for i in range(65536) if pick_uuid("overture-connector-"+str(i))}
        self.assertGreater(len(selection),180)
        self.assertLess(len(selection),350)
        self.assertEqual(selection,{i for i in range(65536) if pick_uuid("overture-connector-"+str(i))})

    def test_same_true_WGS84_connector_at_same_location_is_zero(self):
        lineA=LineString([(14.0,46.0),(15.0,46.0)])
        lineB=LineString([(14.5,45.5),(14.5,46.5)])
        a=lineA.interpolate(.5,normalized=True)
        b=lineB.interpolate(.5,normalized=True)
        self.assertAlmostEqual(Geod(ellps="WGS84").inv(a.x,a.y,b.x,b.y)[2],0,places=6)

    def test_true_source_at_is_WGS84_geodetic_fraction_not_planar_lonlat(self):
        geod=Geod(ellps="WGS84")
        curved=LineString([(0.0,68.0),(7.0,68.0),(7.0,75.0)])
        source_pos=true_WGS84_geodesic_at(curved,0.5,geod)
        planar_pos=curved.interpolate(0.5,normalized=True)
        mismatch_m=abs(geod.inv(source_pos[0],source_pos[1],planar_pos.x,planar_pos.y)[2])
        self.assertGreater(mismatch_m,10000.0)
        a=true_WGS84_geodesic_at(curved,0.0,geod)
        b=true_WGS84_geodesic_at(curved,1.0,geod)
        self.assertAlmostEqual(a[0],0.0)
        self.assertAlmostEqual(a[1],68.0)
        self.assertAlmostEqual(b[0],7.0)
        self.assertAlmostEqual(b[1],75.0)

    def test_real_source_mismatches_should_be_measured_not_ignored(self):
        lineA=LineString([(14.0,46.0),(15.0,46.0)])
        lineB=LineString([(16.0,46.0),(17.0,46.0)])
        a=lineA.interpolate(.5,normalized=True)
        b=lineB.interpolate(.5,normalized=True)
        d=Geod(ellps="WGS84").inv(a.x,a.y,b.x,b.y)[2]
        self.assertGreater(d,1000)
        self.assertIn(100,THRESHOLDS)
        self.assertIn(1000,THRESHOLDS)
if __name__=="__main__":unittest.main()
