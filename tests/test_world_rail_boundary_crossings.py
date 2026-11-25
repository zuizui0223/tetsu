import unittest
from aggregate_world_rail_boundary_crossings import crossing_edges
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
class TestCrossings(unittest.TestCase):
    def path(self,coords,w=100000):
        return crossing_edges([(X_START+x*w,Y_START+y*w) for x,y in coords],w)
    def test_crossing(self):
        self.assertEqual(self.path([(.5,.5),(1.5,.5)]),{((0,0),(1,0))})
        self.assertEqual(self.path([(.5,.5),(1,.5),(1.5,.5)]),{((0,0),(1,0))})
    def test_not_corner_or_along_boundary(self):
        self.assertEqual(self.path([(.5,.5),(1.5,1.5)]),set())
        self.assertEqual(self.path([(.5,1),(1.5,1)]),set())
    def test_same_cell(self):
        self.assertEqual(self.path([(.2,.2),(.8,.8)]),set())
    def test_multiple_cells(self):
        self.assertEqual(self.path([(.5,.5),(3.5,.5)]),{((0,0),(1,0)),((1,0),(2,0)),((2,0),(3,0))})
if __name__=="__main__":unittest.main()
