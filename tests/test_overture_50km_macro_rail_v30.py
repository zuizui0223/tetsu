"""Offline checks for outcome-blind 50km macro physical rail connectivity."""
import sys
from pathlib import Path
import unittest
from shapely.geometry import LineString
from shapely.wkb import dumps
from pyproj import Transformer

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_overture_50km_macro_rail_v30 import (
    macro_geometry,bbox_with_halo,query_macro_rail,
    source_rail_macro,run
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS, GRID_CRS


def geometry_for_cells(block,a,b):
    _,cells=macro_geometry(block)
    trans=Transformer.from_crs(GRID_CRS,4326,always_xy=True)
    c1=cells[a].centroid;c2=cells[b].centroid
    return dumps(LineString([trans.transform(c1.x,c1.y),
                             trans.transform(c2.x,c2.y)]))


def fake_source(block,a,b,uid):
    return (uid,"standard_gauge",[
        {"connector_id":uid+"start","at":0.0},
        {"connector_id":uid+"end","at":1.0}
    ],geometry_for_cells(block,a,b))


class MacroRail50kmTests(unittest.TestCase):
    def test_source_50km_grid_and_country_stratum(self):
        p,cells=macro_geometry("jp_tohoku")
        self.assertEqual(len(cells),16)
        self.assertAlmostEqual(p.area/1e6,40000.)
        self.assertEqual(cells["jp_tohoku_x0_y0"].area/1e6,2500.)
        with self.assertRaises(ValueError):
            macro_geometry("not_frozen")

    def test_source_halo_query_covers_macro_window(self):
        w,s,e,n=bbox_with_halo("jp_tohoku")
        self.assertLess(w,BLOCKS["jp_tohoku"][1])
        self.assertGreater(e,BLOCKS["jp_tohoku"][1])
        self.assertLess(s,BLOCKS["jp_tohoku"][2])
        self.assertGreater(n,BLOCKS["jp_tohoku"][2])
        q=query_macro_rail("jp_tohoku")
        self.assertIn("2026-10-21.0",q)
        self.assertIn("subtype='rail'",q)
        self.assertIn("bbox.xmin",q)

    def test_real_connectors_inside_50km_cells_define_physical_link(self):
        block="jp_tohoku"
        a,b="jp_tohoku_x0_y0","jp_tohoku_x1_y0"
        out=source_rail_macro(block,[fake_source(block,a,b,"one")])
        self.assertEqual(out["n_50km_equal_area_source_cells"],16)
        self.assertEqual(len(out["cell_pair_rail_connectivity"]),120)
        self.assertEqual(out["n_equal_area_cells_with_rail_passage"],2)
        self.assertEqual(out["n_rail_rail_pairs_with_both_connector"],1)
        self.assertEqual(out["n_rail_rail_pairs_physically_connected"],1)
        self.assertEqual(out["n_rail_rail_pairs_not_connected_within_halo"],0)
        self.assertFalse(out["exposure_positivity_gate_min3_each"])

    def test_two_parallel_segments_do_not_forge_same_component(self):
        block="us_midwest"
        r1=fake_source(block,"us_midwest_x0_y0","us_midwest_x1_y0","r1")
        r2=fake_source(block,"us_midwest_x0_y1","us_midwest_x1_y1","r2")
        out=source_rail_macro(block,[r1,r2])
        self.assertEqual(out["n_mapped_rail_components_in_75km_halo"],2)
        self.assertEqual(out["n_rail_rail_pairs_physically_connected"],2)
        self.assertEqual(out["n_rail_rail_pairs_not_connected_within_halo"],4)

    def test_bad_duplicate_source_does_not_get_biological_effect(self):
        block="jp_tohoku"
        one=fake_source(block,"jp_tohoku_x0_y0","jp_tohoku_x1_y0","same")
        with self.assertRaisesRegex(ValueError,"Invalid/duplicated"):
            source_rail_macro(block,[one,one])
        with self.assertRaisesRegex(ValueError,"No Overture"):
            source_rail_macro(block,[])

    def test_nonpredefined_block_selection_fails(self):
        with self.assertRaises(ValueError):
            run(["madeup"])
        with self.assertRaises(ValueError):
            run([])


if __name__=="__main__":
    unittest.main()
