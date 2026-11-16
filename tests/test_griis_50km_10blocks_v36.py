"""Offline source mask provenance checks for 10-block GRIIS macro expansion."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from harvest_griis_50km_10blocks_v36 import apply_country_quality
from audit_gbif_equalarea_macroblocks_v29 import make_grid


def sources(block="jp_tohoku"):
    old={};mask={}
    for i,c in enumerate(make_grid(block)):
        cid=c["cell_id"]
        old[cid]={
            "country":"JP",
            "period_coverage":{
                "2010_2017":{"source_eligible":True},
                "2018_2025":{"source_eligible":i!=2},
            },
        }
        allowed=i not in (0,1)
        mask[cid]={
            "country_fraction_min0p80":allowed,
            "country_and_prior_GBIF_source_qualified":allowed and i!=2,
        }
    return (
        {"complete":True,"blocks":{block:{"cells":old}},"source_errors":{}},
        {"complete":True,"errors":{},"regions":{block:{"country":"JP","cells":mask}}},
    )


class CountryAndAlienMacroSourceTests(unittest.TestCase):
    def test_country_mask_does_not_reclassify_original_source_as_biological_zeros(self):
        old,mask=sources()
        new,trace=apply_country_quality(old,mask,["jp_tohoku"])
        self.assertTrue(old["blocks"]["jp_tohoku"]["cells"]["jp_tohoku_x0_y0"]["period_coverage"]["2010_2017"]["source_eligible"])
        self.assertFalse(new["blocks"]["jp_tohoku"]["cells"]["jp_tohoku_x0_y0"]["period_coverage"]["2010_2017"]["source_eligible"])
        self.assertEqual(trace["jp_tohoku"]["n_source_bins_excluded_for_country_fraction"],4)
        self.assertEqual(trace["jp_tohoku"]["n_two_period_source_bins_before_country_mask"],31)
        self.assertEqual(trace["jp_tohoku"]["n_two_period_source_bins_after_country_mask"],27)
        self.assertEqual(trace["jp_tohoku"]["n_source_cells_country_80pct_and_both_period_eligible"],13)

    def test_mask_and_source_region_alignment(self):
        old,mask=sources()
        mask["regions"]["jp_tohoku"]["cells"].pop("jp_tohoku_x0_y0")
        with self.assertRaisesRegex(ValueError,"geographic cell"):
            apply_country_quality(old,mask,["jp_tohoku"])

    def test_missing_country_source_is_not_a_real_zero(self):
        old,mask=sources();mask["complete"]=False
        with self.assertRaisesRegex(ValueError,"Incomplete independent"):
            apply_country_quality(old,mask,["jp_tohoku"])
        old,mask=sources();old["complete"]=False
        with self.assertRaisesRegex(ValueError,"Incomplete original"):
            apply_country_quality(old,mask,["jp_tohoku"])

    def test_non_predeclared_sample_cannot_be_selected_by_outcomes(self):
        old,mask=sources()
        with self.assertRaisesRegex(ValueError,"Non-frozen"):
            apply_country_quality(old,mask,["outside_region"])
        with self.assertRaisesRegex(ValueError,"Non-frozen"):
            apply_country_quality(old,mask,[])


if __name__=="__main__":
    unittest.main()
