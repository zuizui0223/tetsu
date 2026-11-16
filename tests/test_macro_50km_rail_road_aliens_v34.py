"""Synthetic numerical and provenance gates for two-block 50km real-source workflow."""
import sys
from pathlib import Path
from itertools import combinations
import unittest
import math
import copy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_macro_50km_rail_road_aliens_v34 import (
    partial_rank_r,_eligible_both_period_cells,analyze_one_block,run,
    MACRO_BLOCKS,STRATEGIES
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS


def inputs():
    alien={"complete":True,"regions":{}}
    rail={"complete":True,"regions":{}}
    road={"complete":True,"regions":{}}
    for region in MACRO_BLOCKS:
        country=BLOCKS[region][0]
        city_cells={}
        for i in range(16):
            x=i%4;y=i//4
            cid=f"{region}_x{x}_y{y}"
            count={str(100+j):15 for j in range(3)}
            count[str(200+i)]=20
            count[str(300+i%4)]=25
            per={}
            for period in ("2010_2017","2018_2025"):
                per[period]={
                    "full_facet_source_pagination_completed":True,
                    "alien_observed_minimum_5_species_and_20_records":True,
                    "observed_alien_specieskey_counts":count,
                    "n_indexed_GRIIS_introduced_plant_records":90,
                    "indexed_all_plantae_source_records":150+i*7,
                    "prior_2026_GBIF_indexed_source_record_count":150+i*7
                }
            city_cells[cid]={
                "source_pre_outcome_both_period_recording_eligible":True,
                "periods":per,
            }
        alien["regions"][region]={"country":country,"cells":city_cells}
        railpairs=[];roadpairs=[]
        for i,(a,b) in enumerate(combinations(sorted(city_cells),2)):
            geod=50*(1+i%5)
            roadratio=1.10+.025*(i%11)
            railratio=1.33+.075*(i%13)
            railpairs.append({
                "a":a,"b":b,
                "rail_to_geodesic_ratio":railratio,
                "cell_centroid_snap_augmented_rail_path_km":geod*railratio,
                "source_cell_centroid_geodesic_km":geod,
            })
            roadpairs.append({
                "a":a,"b":b,
                "road_to_geodesic_ratio":roadratio,
                "centroid_snap_augmented_road_path_km":geod*roadratio
            })
        rail["regions"][region]={"routes":{strategy:{"pair_routes":copy.deepcopy(railpairs)}
            for strategy in STRATEGIES}}
        road["regions"][region]={"routes":{strategy:{"pair_paths":copy.deepcopy(roadpairs)}
            for strategy in STRATEGIES}}
    return alien,rail,road


class EqualAreaSourceComparisonTests(unittest.TestCase):
    def test_tied_ranks_and_no_unsupported_small_sample_rho(self):
        x=list(range(40))
        y=list(reversed(range(40)))
        z=[[math.sin(i*.31),math.cos(i*.29)] for i in range(40)]
        v=partial_rank_r(x,y,z)
        self.assertLess(v,-.99)
        self.assertIsNone(partial_rank_r(x[:10],y[:10],z[:10]))

    def test_both_period_eligibility_and_drift(self):
        a,_,_=inputs();x=_eligible_both_period_cells("us_midwest",a)
        self.assertEqual(len(x),16)
        k="us_midwest_x0_y0"
        a["regions"]["us_midwest"]["cells"][k]["periods"]["2018_2025"]["indexed_all_plantae_source_records"]=500
        y=_eligible_both_period_cells("us_midwest",a)
        self.assertEqual(len(y),15)
        self.assertNotIn(k,y)

    def test_real_source_join_cardinality_and_both_anchors(self):
        a,b,c=inputs()
        out=analyze_one_block("fr_burgundy",a,b,c)
        self.assertEqual(out["n_source_qualified_two_period_plant_cells"],16)
        for period in out["periods"]:
            for strat in STRATEGIES:
                x=out["periods"][period][strat]
                self.assertEqual(x["n_plant_pair_routes_with_both_network_lengths"],120)
                self.assertEqual(x["n_plant_qualified_cells_with_complete_macro_routes"],16)
                self.assertEqual(len(x["network_rail_over_road_path_ratio_quantiles"]),3)
                self.assertNotIn("p_value",x)
        self.assertTrue(out["periods"]["2018_2025"]["nearest_in_cell"]["responses"][
            "source_expected_rarefied20_intersection_union_ratio"]["no_p_value_due_to_shared_cells_and_two_macro_blocks"])

    def test_false_imputed_network_routes_are_omitted_not_zero(self):
        a,b,c=inputs()
        r=b["regions"]["us_midwest"]["routes"]["nearest_in_cell"]["pair_routes"][0]
        r["rail_to_geodesic_ratio"]=None
        out=analyze_one_block("us_midwest",a,b,c)
        self.assertEqual(out["periods"]["2010_2017"]["nearest_in_cell"]["n_plant_pair_routes_with_both_network_lengths"],119)
        self.assertEqual(out["periods"]["2010_2017"]["nearest_largest_component_in_cell"]["n_plant_pair_routes_with_both_network_lengths"],120)

    def test_incomplete_sources_and_mismatched_pair_frame_fail_closed(self):
        a,b,c=inputs();a["complete"]=False
        with self.assertRaisesRegex(ValueError,"Incomplete required"):
            run(a,b,c)
        a,b,c=inputs();c["regions"]["us_midwest"]["routes"]["nearest_in_cell"]["pair_paths"].pop()
        with self.assertRaises(ValueError):
            analyze_one_block("us_midwest",a,b,c)

    def test_two_blocks_are_not_misrepresented_as_replication(self):
        a,b,c=inputs()
        out=run(a,b,c)
        self.assertEqual(out["n_independent_geographic_blocks"],2)
        self.assertFalse(out["independent_geographic_replication_sufficient"])
        self.assertTrue(out["NO_causal_propagule_or_evolution_inference"])
        self.assertEqual(set(out["regions"]),set(MACRO_BLOCKS))


if __name__=="__main__":
    unittest.main()
