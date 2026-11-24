#!/usr/bin/env python3
"""v352: frozen six-cell selection and clipped original railway geometry tests."""
import sys,unittest,math,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_global_landcell_Overture_rail_v352 import (
    wave_a,projected_cell_square,clipped_geodesic_rail_length,
    SOURCE_FRAME_SHA,RAIL_RELEASE,CELL_M
)

FRAME=Path(__file__).resolve().parents[1]/"data/derived/v351_global_land_climate_100km_probabilistic_source_verified.json"

class V352FrozenRailSourceTests(unittest.TestCase):
    def test_original_frozen_world_frame_exactly_six_independently_of_rails(self):
        j=json.loads(FRAME.read_text())
        cells=wave_a(j)
        self.assertEqual(j["sampled_cell_unique_hash"],SOURCE_FRAME_SHA)
        self.assertEqual(len(cells),6)
        self.assertEqual(len({x["continent"] for x in cells}),6)
        self.assertEqual(len({x["source_original_fixed_grid_index"] for x in cells}),6)
        self.assertEqual(RAIL_RELEASE,"2026-11-18.0")
        for cell in cells:
            self.assertGreater(cell["source_wave_A_cell_inclusion_probability"],0)
            self.assertLessEqual(cell["source_wave_A_cell_inclusion_probability"],1)
            self.assertAlmostEqual(cell["source_wave_A_cell_inclusion_probability"],
                                   1/cell["source_parent_population_n"])

    def test_original_probability_source_mutation_rejected(self):
        j=json.loads(FRAME.read_text())
        j["sampled_cell_unique_hash"]="bad"
        with self.assertRaises(ValueError):wave_a(j)

    def test_exact_100km_cell_vs_enclosing_query_bbox(self):
        c=wave_a(json.loads(FRAME.read_text()))[0]
        poly,bbox=projected_cell_square(c)
        self.assertAlmostEqual(poly.area,CELL_M*CELL_M)
        self.assertLess(bbox[0],bbox[2])
        self.assertLess(bbox[1],bbox[3])

    def test_railway_outside_cell_although_inside_bbox_is_not_counted(self):
        from pyproj import Transformer,Geod
        from shapely.geometry import LineString
        from shapely.ops import transform
        c=wave_a(json.loads(FRAME.read_text()))[0]
        poly,_=projected_cell_square(c)
        x0,y0,x1,y1=poly.bounds
        inv=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
        fwd=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
        geod=Geod(ellps="WGS84")
        outside=LineString([(x0-2000,y0+1200),(x0-500,y0+1200)])
        within=LineString([(x0-1000,y0+2000),(x0+1000,y0+2000)])
        absent,npart=clipped_geodesic_rail_length(transform(inv.transform,outside),poly,fwd,inv,geod)
        present,npart2=clipped_geodesic_rail_length(transform(inv.transform,within),poly,fwd,inv,geod)
        self.assertAlmostEqual(absent,0.0,places=4)
        self.assertEqual(npart,0)
        self.assertGreater(present,300)
        self.assertEqual(npart2,1)

if __name__=="__main__":unittest.main()
