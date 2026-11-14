"""Synthetic country GRIIS × source dataset overlap × physical RAIL/ROAD path tests."""
from __future__ import annotations
import sys
from pathlib import Path
from itertools import combinations
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_dataset_overlap_adjusted_rail_plants_v51 import (
    dataset_overlap,merge_dataset_chunks,valid_dataset_profiles,
    compute_country,analyze,RESPONSES,STRATEGIES,PERIODS,BLOCKS_SIX
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS


def fake_sources():
    flora={"complete":True,"regions":{},"source_errors":{}}
    rail={"complete":True,"regions":{}}
    road={"complete":True,"regions":{}}
    chunks=[]
    for block in BLOCKS_SIX:
        country=BLOCKS[block][0]
        pflora={};dataset_cells={}
        for i in range(16):
            cid=f"{block}_x{i%4}_y{i//4}"
            counts={str(101+k):30+i%7 for k in range(3)}
            counts[str(800+i)]=25
            counts[str(900+i%5)]=30
            sums=sum(counts.values())
            frow={
                "source_pre_outcome_both_period_recording_eligible":True,
                "country":country,
                "cell_id":cid,
                "periods":{p:{
                    "observed_alien_specieskey_counts":counts,
                    "full_facet_source_pagination_completed":True,
                    "alien_observed_minimum_5_species_and_20_records":True,
                    "n_indexed_GRIIS_introduced_plant_records":sums,
                    "indexed_all_plantae_source_records":1000+i*20,
                    "prior_2026_GBIF_indexed_source_record_count":1000+i*20
                } for p in PERIODS}}
            pflora[cid]=frow
            pa=700-i*10
            pb=250+i*5
            pc=1000-pa-pb
            d={
                "dataset_key_record_counts":{"dataset_A":pa,"dataset_B":pb,f"dataset_{i%4}":pc},
                "source_indexed_all_Plantae_records":1000,
                "source_profile_valid_for_recording_bias_covariates":True,
                "top_datasetKey_source_share":pa/1000,
                "top_datasetKey_source_id":"dataset_A",
                "source_dataset_facet_fully_paged":True,
            }
            dataset_cells[cid]={
                "country":country,
                "original_GRIIS_two_period_source_eligible":True,
                "periods":{p:dict(d) for p in PERIODS}
            }
        flora["regions"][block]={"country":country,"cells":pflora}
        railpairs=[];roadpairs=[]
        for j,(a,b) in enumerate(combinations(sorted(pflora),2)):
            geo=50+10*(j%5)
            rail_ratio=1.21+.06*(j%13)
            road_ratio=1.07+.03*(j%11)
            railpairs.append({"a":a,"b":b,
                "rail_to_geodesic_ratio":rail_ratio,
                "source_cell_centroid_geodesic_km":geo,
                "cell_centroid_snap_augmented_rail_path_km":geo*rail_ratio})
            roadpairs.append({"a":a,"b":b,
                "road_to_geodesic_ratio":road_ratio,
                "centroid_snap_augmented_road_path_km":geo*road_ratio})
        rail["regions"][block]={"country":country,
            "routes":{s:{"pair_routes":[dict(p) for p in railpairs]} for s in STRATEGIES}}
        road["regions"][block]={"source":{"country":country},
            "routes":{s:{"pair_paths":[dict(p) for p in roadpairs]} for s in STRATEGIES}}
        chunks.append({"complete":True,"errors":{},"block":block,
                       "country":country,"source":"synthetic dataset source",
                       "cells":dataset_cells})
    return flora,rail,road,chunks


class DatasetOverlapAdjustTests(unittest.TestCase):
    def test_weighted_shared_dataset_fraction_not_taxon_presence(self):
        a={"dataset_key_record_counts":{"x":70,"y":30},"source_indexed_all_Plantae_records":100}
        b={"dataset_key_record_counts":{"x":50,"z":50},"source_indexed_all_Plantae_records":100}
        self.assertAlmostEqual(dataset_overlap(a,b),.5)
        c={"dataset_key_record_counts":{"y":30},"source_indexed_all_Plantae_records":100}
        self.assertAlmostEqual(dataset_overlap(a,c),.3)

    def test_six_actual_source_territories_are_fixed_not_selected_by_plant_sign(self):
        flora,rail,road,chunks=fake_sources()
        joined=merge_dataset_chunks(chunks)
        self.assertEqual(set(joined["regions"]),set(BLOCKS_SIX))
        self.assertEqual(joined["n_unique_country_blocks_with_true_source_dataset_composition"],6)
        with self.assertRaisesRegex(ValueError,"duplicate"):
            merge_dataset_chunks(chunks+[chunks[0]])

    def test_invalid_dataset_source_does_not_get_code_as_zero_plants(self):
        flora,rail,road,chunks=fake_sources()
        block="us_midwest"
        cid=f"{block}_x0_y0"
        chunks[0]["cells"][cid]["periods"]["2018_2025"]["source_profile_valid_for_recording_bias_covariates"]=False
        joined=merge_dataset_chunks(chunks)
        cells,ex=valid_dataset_profiles(block,flora,joined)
        self.assertEqual(len(cells),15)
        self.assertIn(cid,ex)
        summary=compute_country(block,flora,rail,road,joined)
        self.assertEqual(summary["source_cells_excluded_due_to_dataset_marginal_quality"],1)
        self.assertEqual(summary["periods"]["2018_2025"]["nearest_in_cell"]["n_dependent_50km_cell_pair_rows_with_rail_road_and_dataset"],105)

    def test_reuses_same_pairs_for_before_after_source_covariate(self):
        flora,rail,road,chunks=fake_sources()
        joined=merge_dataset_chunks(chunks)
        row=compute_country("fr_burgundy",flora,rail,road,joined)
        self.assertEqual(row["n_true_dataset_quality_qualified_cells"],16)
        for p in PERIODS:
            for st in STRATEGIES:
                result=row["periods"][p][st]
                self.assertEqual(result["n_dependent_50km_cell_pair_rows_with_rail_road_and_dataset"],120)
                for resp in RESPONSES:
                    vals=result["responses"][resp]
                    self.assertIn("partial_rho_add_allPlantae_dataset_composition_overlap",vals)
                    self.assertIn("partial_rho_original_controls_on_IDENTICAL_dataset_supported_cell_pairs",vals)
                    self.assertTrue(vals["NO_causal_source_bias_adjustment_or_valid_pairwise_p_value"])

    def test_country_equal_block_aggregate_not_fake_pair_sample(self):
        flora,rail,road,chunks=fake_sources()
        joined=merge_dataset_chunks(chunks)
        report=analyze(flora,rail,road,joined)
        self.assertEqual(report["n_predeclared_country_macro_blocks"],6)
        self.assertTrue(report["source_dataset_marginal_ALL_Plantae_NOT_only_aliens"])
        main=report["country_equal_aggregates"]["2018_2025"]["nearest_in_cell"][RESPONSES[0]]
        self.assertIn("partial_rho_add_allPlantae_dataset_composition_overlap",main)
        self.assertEqual(main["partial_rho_add_allPlantae_dataset_composition_overlap"]["n_independent_qualifying_country_blocks"],
                         len(main["partial_rho_add_allPlantae_dataset_composition_overlap"]["regional_coefficients"]))
        self.assertTrue(main["partial_rho_add_allPlantae_dataset_composition_overlap"]["not_pooled_nonindependent_pair_pvalue"])

    def test_missing_real_country_dataset_source_fails_without_inventing_overlap(self):
        flora,rail,road,chunks=fake_sources()
        with self.assertRaises(ValueError):
            merge_dataset_chunks(chunks[:-1])
        data=merge_dataset_chunks(chunks)
        data["regions"]["ca_ontario"]["cells"].pop("ca_ontario_x0_y0")
        with self.assertRaisesRegex(ValueError,"cell geometry mismatch"):
            analyze(flora,rail,road,data)


if __name__=="__main__":
    unittest.main()
