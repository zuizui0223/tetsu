"""Synthetic source integrity and exact GBIF indexed rarefaction depth checks."""
from __future__ import annotations
import sys
from pathlib import Path
from itertools import combinations
from math import comb

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import unittest
from analyze_rarefaction_depth_6country_v47 import (
    stable_record_rarefaction_probability,stable_probabilities,
    expected_overlap_ratio,calculate_block,aggregate_country,analyze,
    DEPTHS,FIXED_COHORT_MIN_DEPTH
)
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX,PERIODS,STRATEGIES
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS


def sources():
    flora={"complete":True,"regions":{},"source_errors":{}}
    rail={"complete":True,"regions":{}}
    road={"complete":True,"regions":{}}
    for block in BLOCKS_SIX:
        cells={}
        for i in range(16):
            cid=f"{block}_x{i%4}_y{i//4}"
            scale=1 if i<4 else (2 if i<10 else 4)
            counts={
                "speciesA":(46+3*i)*scale,
                "speciesB":(30+2*i)*scale,
                "speciesC":(30+(i%5))*scale,
                "speciesD":(20+(i%3))*scale,
                f"species_{i}":(11+(i%4))*scale
            }
            total=sum(counts.values())
            per={}
            for period in PERIODS:
                per[period]={
                    "full_facet_source_pagination_completed":True,
                    "alien_observed_minimum_5_species_and_20_records":True,
                    "n_indexed_GRIIS_introduced_plant_records":total,
                    "observed_alien_specieskey_counts":counts,
                    "indexed_all_plantae_source_records":1500+i*21,
                    "prior_2026_GBIF_indexed_source_record_count":1500+i*21,
                }
            cells[cid]={
                "source_pre_outcome_both_period_recording_eligible":True,
                "periods":per,
            }
        flora["regions"][block]={"country":BLOCKS[block][0],"cells":cells}
        pairs=list(combinations(sorted(cells),2))
        rp=[];dp=[]
        for j,(a,b) in enumerate(pairs):
            d=50*(1+j%5)
            railratio=1.2+.03*(j%11)
            roadratio=1.12+.02*(j%13)
            rp.append({"a":a,"b":b,"rail_to_geodesic_ratio":railratio,
                       "source_cell_centroid_geodesic_km":d,
                       "cell_centroid_snap_augmented_rail_path_km":d*railratio})
            dp.append({"a":a,"b":b,"road_to_geodesic_ratio":roadratio,
                       "centroid_snap_augmented_road_path_km":d*roadratio})
        rail["regions"][block]={"country":BLOCKS[block][0],
                    "routes":{s:{"pair_routes":[dict(r) for r in rp]} for s in STRATEGIES}}
        road["regions"][block]={"source":{"country":BLOCKS[block][0]},
                    "routes":{s:{"pair_paths":[dict(r) for r in dp]} for s in STRATEGIES}}
    return flora,rail,road


class RarefactionDepthTests(unittest.TestCase):
    def test_stable_hypergeometric_probability_matches_combinatorics(self):
        p=stable_record_rarefaction_probability(12,3,5)
        self.assertAlmostEqual(p,1-comb(9,5)/comb(12,5),places=13)
        self.assertAlmostEqual(stable_record_rarefaction_probability(12,10,5),1)
        self.assertGreater(stable_record_rarefaction_probability(1_000_000,1,100),0.)
        self.assertAlmostEqual(stable_record_rarefaction_probability(1_000_000,1,100),.0001,delta=1e-8)

    def test_depth_monotonic_and_no_sampling_more_than_source(self):
        counts={"A":80,"B":20,"C":4}
        p20=stable_probabilities(counts,20)
        p80=stable_probabilities(counts,80)
        self.assertTrue(all(p80[k]>=p20[k] for k in counts))
        self.assertGreater(expected_overlap_ratio(p20,p20),0.)
        with self.assertRaises(ValueError):
            stable_probabilities(counts,200)
        with self.assertRaises(ValueError):
            stable_record_rarefaction_probability(50,0,20)
        with self.assertRaises(ValueError):
            stable_record_rarefaction_probability(50,20,51)

    def test_source_all_six_blocks_predeclared_not_speculative(self):
        flora,rail,road=sources()
        x=analyze(flora,rail,road)
        self.assertEqual(len(x["country_results"]),6)
        self.assertEqual(x["depths"],list(DEPTHS))
        self.assertTrue(x["not_a_causal_railway_seed_effect"])
        self.assertTrue(x["posthoc_sensitivity_not_independent_preregistered_test"])

    def test_fixed80_matches_same_cells_at_every_shallow_depth(self):
        flora,rail,road=sources()
        block=BLOCKS_SIX[0]
        result=calculate_block(block,flora,rail,road)
        for anchor in STRATEGIES:
            fixed=result["results"][anchor]["fixed80"]
            cells=[fixed[str(k)]["by_period"]["2018_2025"]["n_eligible_cells"] for k in (20,40,80)]
            self.assertEqual(cells,[16,16,16])
            self.assertTrue(fixed["160"]["unsupported_above_fixed80_source_depth"])
            self.assertTrue(fixed["320"]["unsupported_above_fixed80_source_depth"])
        self.assertEqual(result["n_original_qualifying_cells_K20"],16)
        self.assertEqual(result["n_fixed_K80_qualifying_cells"],16)

    def test_dynamic_depth_can_change_support_and_never_invent_nulls(self):
        flora,rail,road=sources()
        row=calculate_block(BLOCKS_SIX[0],flora,rail,road)
        dyn=row["results"]["nearest_in_cell"]["depth_specific"]
        n20=dyn["20"]["by_period"]["2018_2025"]["n_eligible_cells"]
        n320=dyn["320"]["by_period"]["2018_2025"]["n_eligible_cells"]
        self.assertEqual(n20,16)
        self.assertLess(n320,n20)
        self.assertEqual(dyn["320"]["by_period"]["2018_2025"]["n_unique_cells_with_at_least_one_network_pair"],n320)

    def test_country_equal_weight_replicate_not_pairwise_pseudoreplicates(self):
        flora,rail,road=sources()
        result=analyze(flora,rail,road)
        a=result["country_equal_weight"]["2018_2025"]["nearest_in_cell"]["fixed80"]["20"]
        self.assertEqual(a["n_independent_supported_country_blocks"],6)
        self.assertEqual(len(a["country_rhos"]),6)
        self.assertTrue(a["min_5_country_source_coverage_gate"])
        self.assertNotIn("p_value",a)
        self.assertTrue(a["no_p_value_or_train_propagule_inference"])

    def test_missing_network_or_source_coverage_fails_closed(self):
        flora,rail,road=sources()
        road["regions"].pop("ca_ontario")
        with self.assertRaises(ValueError):
            analyze(flora,rail,road)
        flora,rail,road=sources()
        flora["regions"]["us_midwest"]["cells"]["us_midwest_x0_y0"]["periods"]["2018_2025"]["n_indexed_GRIIS_introduced_plant_records"]=-10
        with self.assertRaises(ValueError):
            analyze(flora,rail,road)


if __name__=="__main__":
    unittest.main()
