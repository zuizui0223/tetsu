#!/usr/bin/env python3
"""v370 original worldwide railway GIS source cell clipping and sum tests."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from aggregate_entire_world_GERS_rail_100km_200km_grid_v370 import (
  tile_index,clip_source_line_to_tiles,rollup,CELL_SIZE_M,X_START,Y_START,
  MAX_ALLOWED_MAPPED_GERS_RAIL_GRID_KM_MISMATCH
)

class V370WholeWorldRealRailGridTests(unittest.TestCase):
    def setUp(self):
        from pyproj import Transformer,Geod
        self.inv=Transformer.from_crs('EPSG:6933','EPSG:4326',always_xy=True)
        self.fwd=Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True)
        self.geod=Geod(ellps='WGS84')
    def test_full_equatorial_and_cross_dateline_source_grid_coverage(self):
        self.assertEqual(CELL_SIZE_M,100000)
        self.assertEqual(MAX_ALLOWED_MAPPED_GERS_RAIL_GRID_KM_MISMATCH,2.0)
        self.assertEqual(tile_index(X_START+123456,Y_START+123456),(1,1))
        self.assertIsNone(tile_index(float('nan'),0))
    def test_source_segment_fully_inside_one_original_100km_cell_is_not_rounded(self):
        from shapely.geometry import LineString
        x=X_START+180*100000;y=Y_START+100*100000
        a=self.inv.transform(x+12000,y+16000)
        b=self.inv.transform(x+34000,y+17000)
        src=LineString([a,b])
        L=abs(self.geod.geometry_length(src))
        result,missing,status=clip_source_line_to_tiles(
            src,L,self.fwd,self.inv,self.geod)
        self.assertEqual(status,'source_one_tile')
        self.assertEqual(len(result),1)
        self.assertAlmostEqual(sum(result.values()),L,places=4)
    def test_true_100km_cell_boundary_splits_geodesic_source_track_length(self):
        from shapely.geometry import LineString
        x=X_START+180*100000;y=Y_START+100*100000
        a=self.inv.transform(x+99000,y+35000)
        b=self.inv.transform(x+101000,y+35000)
        src=LineString([a,b])
        L=abs(self.geod.geometry_length(src))
        result,missing,status=clip_source_line_to_tiles(
            src,L,self.fwd,self.inv,self.geod)
        self.assertEqual(status,'source_multitile_clip')
        self.assertEqual(len(result),2)
        self.assertAlmostEqual(sum(result.values()),L,delta=.5)
    def test_dateline_must_not_draw_fake_35000km_railway_across_globe(self):
        from shapely.geometry import LineString
        track=LineString([(179.8,50.0),(-179.8,50.0)])
        L=abs(self.geod.geometry_length(track))
        result,missing,status=clip_source_line_to_tiles(
            track,L,self.fwd,self.inv,self.geod)
        self.assertEqual(status,'source_dateline_crossing')
        self.assertEqual(result,{})
        self.assertAlmostEqual(missing,L)
    def test_200km_rollup_matches_four_100km_original_source_landcells(self):
        cells={}
        for i in (2,3):
            for j in (4,5):
                cells[i,j]={'source_rail_length_metres':1000.,'original_source_GERS_segmentation_instances':1}
        out=rollup(cells,200000)
        self.assertEqual(len(out),1)
        self.assertEqual(out[0]['source_equal_area_grid_ix'],1)
        self.assertEqual(out[0]['source_equal_area_grid_iy'],2)
        self.assertAlmostEqual(out[0]['source_real_original_Overture_normal_rail_km_clipped_into_cell'],4.)

if __name__=='__main__':unittest.main()
