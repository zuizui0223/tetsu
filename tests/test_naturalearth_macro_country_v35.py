"""Synthetic land-jurisdiction source checks for true 50km equal-area blocks."""
import sys
from pathlib import Path
from shapely.ops import unary_union
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_naturalearth_macro_country_v35 import (
    _cell_geometry,summarize_block,source_countries,COUNTRY_FRAC_MIN
)


def fake_coverage(block):
    cells={}
    for cid in _cell_geometry(block):
        cells[cid]={
            "country":"JP",
            "period_coverage":{
                "2010_2017":{"source_eligible":True},
                "2018_2025":{"source_eligible":True}
            }
        }
    return {"complete":True,"blocks":{block:{"cells":cells}}}


class MacroNaturalEarthCountryTests(unittest.TestCase):
    def test_equal_area_source_cells_preserved(self):
        cells=_cell_geometry("jp_tohoku")
        self.assertEqual(len(cells),16)
        self.assertTrue(all(round(g.area/1e6)==2500 for g in cells.values()))
        self.assertAlmostEqual(sum(x.area for x in cells.values())/1e6,40000)

    def test_source_country_polygon_fraction_without_biological_outcome(self):
        block="jp_tohoku"
        frame=_cell_geometry(block)
        covered=unary_union([g for k,g in frame.items() if "_y0" in k or "_y1" in k])
        src=fake_coverage(block)
        out=summarize_block(block,covered,src)
        self.assertEqual(out["n_source_country_fraction_ge0p80"],8)
        self.assertEqual(out["n_prior_GBIF_both_period_eligible"],16)
        self.assertEqual(out["n_country_and_GBIF_eligible"],8)
        self.assertEqual(out["n_2500km2_source_cells"],16)
        self.assertEqual(COUNTRY_FRAC_MIN,.8)

    def test_missing_GBIF_source_does_not_invent_observed_plants(self):
        block="jp_tohoku"
        frame=_cell_geometry(block)
        covered=unary_union(list(frame.values()))
        src=fake_coverage(block)
        src["blocks"][block]["cells"][f"{block}_x0_y0"]["period_coverage"]["2018_2025"]["source_eligible"]=False
        out=summarize_block(block,covered,src)
        self.assertEqual(out["n_source_country_fraction_ge0p80"],16)
        self.assertEqual(out["n_country_and_GBIF_eligible"],15)

    def test_invalid_country_zip_and_bad_country_source_fail_closed(self):
        with self.assertRaises(ValueError):
            source_countries(b"corrupt")
        with self.assertRaises(ValueError):
            summarize_block("unknown",unary_union(list(_cell_geometry("jp_tohoku").values())),fake_coverage("jp_tohoku"))

    def test_wrong_country_in_coverage_fails(self):
        block="jp_tohoku";src=fake_coverage(block)
        src["blocks"][block]["cells"][f"{block}_x0_y0"]["country"]="US"
        with self.assertRaisesRegex(ValueError,"Country mismatch"):
            summarize_block(block,unary_union(list(_cell_geometry(block).values())),src)


if __name__=="__main__":
    unittest.main()
