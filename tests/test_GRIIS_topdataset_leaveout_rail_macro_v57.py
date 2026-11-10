"""Test same-pair GBIF source dominant-dataset removal in six-country mapped RAIL/ROAD study."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[0]))
from analyze_GRIIS_topdataset_leaveout_rail_macro_v57 import (
    merge_species_sources,source_eligible_cell_pairs,compute_block,run,
    DEPTHS,COMPARISONS,STRATEGIES,PERIODS
)
from test_dataset_overlap_adjusted_rail_plants_v51 import fake_sources
from analyze_dataset_overlap_adjusted_rail_plants_v51 import merge_dataset_chunks
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX


def fake():
    flora,rail,road,raw_ds=fake_sources()
    by_country=[]
    for block in BLOCKS_SIX:
        cells={}
        for cid,old in flora["regions"][block]["cells"].items():
            source={}
            for period in PERIODS:
                full=old["periods"][period]["observed_alien_specieskey_counts"]
                other={k:max(1,int(v*.8)) for k,v in full.items()}
                source[period]={
                    "full_source_current_country_GRIIS_aliens":dict(full),
                    "remaining_datasets_country_GRIIS_aliens":other,
                    "dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys":True,
                }
            cells[cid]={
                "country":BLOCKS[block][0],
                "original_two_period_country_GRIIS_eligible":True,
                "periods":source
            }
        by_country.append({
            "complete":True,"source_errors":{},
            "source_geo_block":block,"source_country":BLOCKS[block][0],
            "source_cells":cells,
        })
    return flora,rail,road,merge_species_sources(by_country),merge_dataset_chunks(raw_ds)


class AlienDominantDatasetRemovalTests(unittest.TestCase):
    def test_same_source_country_eligible_cells_and_pair_universe(self):
        flora,rail,road,split,ds=fake()
        result=compute_block("us_midwest",flora,rail,road,split,ds)
        for k in DEPTHS:
            for period in PERIODS:
                for anchor in STRATEGIES:
                    x=result["depths"][str(k)][period][anchor]
                    self.assertEqual(x["n_source_cells_with_full_and_without_top_dataset_K_records_in_BOTH_periods"],16)
                    self.assertEqual(x["n_matched_RAIL_ROAD_independently_sourced_pairs_NOT_independent"],120)
                    self.assertEqual(set(x["partial_rank_same_pair_flora_full_and_after_dominant_dataset_removal"]),set(COMPARISONS))
                    self.assertTrue(x["NO_pvalue_or_causal_seed_transport"])

    def test_minK_source_cohort_changes_both_outcomes_together(self):
        flora,rail,road,split,ds=fake()
        cell="us_midwest_x0_y0"
        for period in PERIODS:
            split["regions"]["us_midwest"]["source_cells"][cell]["periods"][period]["remaining_datasets_country_GRIIS_aliens"]={
                "101":50,"102":40
            }
        result=compute_block("us_midwest",flora,rail,road,split,ds)
        x=result["depths"]["20"]["2018_2025"]["nearest_in_cell"]
        self.assertEqual(x["n_source_cells_with_full_and_without_top_dataset_K_records_in_BOTH_periods"],15)
        self.assertEqual(x["n_matched_RAIL_ROAD_independently_sourced_pairs_NOT_independent"],105)
        self.assertEqual(x["n_source_cells_with_K_records_all_datasets_before_top_source_exclusion"],16)

    def test_full_six_country_outputs_keep_absence_of_pairwise_significance(self):
        flora,rail,road,split,ds=fake()
        report=run(flora,rail,road,split,ds)
        self.assertEqual(len(report["regions"]),6)
        self.assertTrue(report["source_data_subtraction_is_not_observer_effort_or_true_plant_absence"])
        for k in DEPTHS:
            row=report["equal_country_summaries"][str(k)]["2018_2025"]["nearest_in_cell"]
            self.assertEqual(set(row),set(COMPARISONS))
            self.assertTrue(row["without_dominant_dataset"]["NO_pairwise_significance_or_causal_seed_claim"])
            self.assertEqual(row["without_dominant_dataset"]["n_source_supported_independent_country_blocks"],
                             len(row["without_dominant_dataset"]["regions"]))

    def test_missing_country_source_or_duplicate_country_block_fails(self):
        flora,rail,road,split,ds=fake()
        split["regions"].pop("fr_burgundy")
        with self.assertRaisesRegex(ValueError,"six-country geographic"):
            run(flora,rail,road,split,ds)
        flora,rail,road,split,ds=fake()
        fragments=[{"complete":True,"source_errors":{},"source_geo_block":b,
            "source_country":BLOCKS[b][0],"source_cells":split["regions"][b]["source_cells"]}
            for b in BLOCKS_SIX]
        with self.assertRaisesRegex(ValueError,"Unexpected/duplicate"):
            merge_species_sources(fragments+[fragments[0]])

    def test_incomplete_country_period_taxon_source_fails_closed(self):
        flora,rail,road,split,ds=fake()
        split["regions"]["ca_ontario"]["source_cells"]["ca_ontario_x0_y0"]["periods"]["2018_2025"].pop(
            "dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys")
        with self.assertRaisesRegex(ValueError,"partition missing"):
            source_eligible_cell_pairs("ca_ontario",flora,split,ds,20)


if __name__=="__main__":
    unittest.main()
