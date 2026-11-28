import unittest
from collections import Counter
from pyproj import Transformer,Geod
from scan_grip4_actual_road_lengths import part_km
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
class Vector(unittest.TestCase):
    def test_same_cell_and_boundary(self):
        inv=Transformer.from_crs(6933,4326,always_xy=True)
        fw=Transformer.from_crs(4326,6933,always_xy=True)
        ident=Transformer.from_crs(4326,4326,always_xy=True)
        geod=Geod(ellps="WGS84")
        x=X_START+87*200000;y=Y_START+37*200000+50000
        part=[inv.transform(x-1000,y),inv.transform(x+1000,y)]
        q=part_km(part,{(86,37),(87,37)},ident,fw,inv,geod,Counter())
        self.assertEqual(set(q),{(86,37),(87,37)})
        self.assertAlmostEqual(sum(q.values()),abs(geod.line_length(*zip(*part)))/1000,places=5)
        r=part_km(part,{(87,37)},ident,fw,inv,geod,Counter())
        self.assertEqual(set(r),{(87,37)})
        self.assertAlmostEqual(r[(87,37)],q[(87,37)])
if __name__=="__main__":unittest.main()
