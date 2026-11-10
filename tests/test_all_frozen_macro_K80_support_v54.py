"""Ten original country macro blocks K80 source opportunity not selected by train results."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_all_frozen_macro_K80_support_v54 import (
    eligible_count_at_depth,analyze,DEPTHS,MIN_INDEPENDENT_SOURCE_CELLS
)
from audit_gbif_equalarea_macroblocks_v29 import PILOT,BLOCKS,make_grid
from analyze_6country_50km_network_plants_v44 import PERIODS


def source_fixture():
    regions={}
    for block in PILOT:
        country=BLOCKS[block][0]
        cells={}
        for i,c in enumerate(make_grid(block)):
            n_each=(5+i)*5
            counts={str(j):n_each for j in range(5)}
            cells[c["cell_id"]]={
                "source_pre_outcome_both_period_recording_eligible":True,
                "periods":{p:{
                    "full_facet_source_pagination_completed":True,
                    "alien_observed_minimum_5_species_and_20_records":True,
                    "observed_alien_specieskey_counts":counts,
                    "n_indexed_GRIIS_introduced_plant_records":5*n_each
                } for p in PERIODS}
            }
        regions[block]={"country":country,"cells":cells}
    return {"complete":True,"source_errors":{},"regions":regions}


class AllFrozenRegionSourceDepthTests(unittest.TestCase):
    def test_8_country_10_original_geographical_blocks(self):
        self.assertEqual(len(PILOT),10)
        self.assertEqual(len({BLOCKS[b][0] for b in PILOT}),8)
        self.assertEqual(DEPTHS,(20,40,80,160,320))

    def test_deeper_source_requires_actual_more_indexed_GRIIS_records(self):
        row=source_fixture()["regions"]["us_southeast"]["cells"]["us_southeast_x0_y0"]
        self.assertTrue(eligible_count_at_depth(row,20))
        self.assertTrue(eligible_count_at_depth(row,80))
        self.assertFalse(eligible_count_at_depth(row,320))
        row["periods"]["2010_2017"]["n_indexed_GRIIS_introduced_plant_records"]=-1
        with self.assertRaisesRegex(ValueError,"count inconsistent"):
            eligible_count_at_depth(row,80)

    def test_country_block_source_nonoverlapping_and_not_pseudoreplicated(self):
        result=analyze(source_fixture())
        self.assertEqual(result["n_source_country_blocks"],10)
        self.assertEqual(result["n_distinct_source_countries"],8)
        self.assertEqual(result["total_source_geographic_cells"],160)
        self.assertEqual(len(result["regions"]),10)
        self.assertTrue(result["complete"])
        for block in PILOT:
            d=result["regions"][block]["by_depth"]
            self.assertTrue(d["20"]["n_original_source_eligible_both_periods"]>=d["80"]["n_original_source_eligible_both_periods"])
            self.assertTrue(d["80"]["n_original_source_eligible_both_periods"]>=d["320"]["n_original_source_eligible_both_periods"])
            self.assertTrue(d["80"]["max_possible_source_pair_comparisons_before_rail_or_road_route_QA"]>=0)
        self.assertEqual(result["source_summary"]["80"]["n_geographical_country_blocks_passing_pre_rail_source_minimum"],10)
        self.assertEqual(MIN_INDEPENDENT_SOURCE_CELLS,8)

    def test_missing_original_country_source_is_never_assigned_zero(self):
        prior=source_fixture()
        prior["regions"].pop("nz_north")
        with self.assertRaisesRegex(ValueError,"Pre-frozen 10"):
            analyze(prior)
        prior=source_fixture()
        prior["complete"]=False
        with self.assertRaisesRegex(ValueError,"incomplete"):
            analyze(prior)


if __name__=="__main__":
    unittest.main()
