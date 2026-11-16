"""v7.6 fixed US train-count × GRIIS-GIFT macro source regression tests.

All fixtures are synthetic. No seed-vector, freight-tonne, mowing, wind,
species establishment, or evolutionary effects are inferred.
"""
import sys, math, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_US_train_counts_GRIIS_traits_v76 import (
    BLOCKS, PERIOD, K, MIN_COUNTRY_CELLS, split_counts, species_metrics,
    build_block, partial_rank_cells, leave_one_cell_out_diagnostic,
)

TRAITS={"11":"herbaceous","12":"herbaceous","13":"herbaceous",
        "21":"woody","22":"woody","23":"woody"}

def fake_sources():
    flora={"country":"US","cells":{}}
    fra={"source_open_closed_balance_verified":True,
         "by_original_50km_cells":{}}
    for y in range(4):
        for x in range(4):
            cid=f"us_midwest_x{x}_y{y}"
            counts={"11":25,"12":20,"13":15,"21":24,"22":18,"23":14}
            flora["cells"][cid]={
                "source_pre_outcome_both_period_recording_eligible":True,
                "periods":{PERIOD:{
                  "full_facet_source_pagination_completed":True,
                  "alien_observed_minimum_5_species_and_20_records":True,
                  "indexed_all_plantae_source_records":1000 + 31*x + 17*y + 4*x*y,
                  "observed_alien_specieskey_counts":counts}}}
            fra["by_original_50km_cells"][cid]={
                "train_opportunity_metrics_USE_confirmed_OPEN_crossings_ONLY":True,
                "n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields":10,
                "n_crossings_confirmed_open_in_current_source":10+x+y,
                "median_reported_total_daily_through_trains_at_eligible_crossings":2+x+3*y+(x*y)%3}
    return flora,fra

class CountAndTraitSafety(unittest.TestCase):
    def test_fixed_two_us_blocks_and_2018_2025_source_period(self):
        self.assertEqual(BLOCKS,("us_midwest","us_southeast"))
        self.assertEqual(PERIOD,"2018_2025")
        self.assertEqual(K,20)
        self.assertEqual(MIN_COUNTRY_CELLS,10)
    def test_missing_trait_is_unknown_not_herbaceous(self):
        allcounts={"11":30,"21":35,"999":40}
        groups=split_counts(allcounts,TRAITS)
        self.assertEqual(groups["all_GRIIS_introduced"],allcounts)
        self.assertEqual(groups["herbaceous"],{"11":30})
        self.assertEqual(groups["woody"],{"21":35})
        self.assertNotIn("999",groups["herbaceous"])
        self.assertNotIn("999",groups["woody"])
    def test_rarefied_richness_not_abundance(self):
        counts={"11":12,"12":12,"13":12}
        m=species_metrics(counts)
        self.assertEqual(m["raw_recorded_species_richness"],3)
        self.assertTrue(0<m["rarefied20_expected_species_richness"]<=3)
    def test_source_cells_are_fixed_16_not_pseudoreplicates(self):
        flora,fra=fake_sources()
        r=build_block("us_midwest",flora,fra,TRAITS)
        self.assertEqual(r["n_original_cells"],16)
        self.assertEqual(r["n_eligible_cells_for_source_observations"],16)
        self.assertTrue(r["within_region_cells_NOT_independent_country_replicates"])
        self.assertEqual(set(r["descriptive_partial_rank_results"]),{
            "all_GRIIS_introduced","herbaceous","woody"})
    def test_closed_crossing_sensitive_source_is_rejected(self):
        flora,fra=fake_sources()
        cid="us_midwest_x1_y1"
        fra["by_original_50km_cells"][cid]["train_opportunity_metrics_USE_confirmed_OPEN_crossings_ONLY"]=False
        with self.assertRaisesRegex(ValueError,"CLOSED"):
            build_block("us_midwest",flora,fra,TRAITS)
    def test_2018to2026_reporting_year_is_not_historical_train_effect(self):
        flora,fra=fake_sources()
        for r in fra["by_original_50km_cells"].values():
            r["n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields"]=0
        out=build_block("us_midwest",flora,fra,TRAITS)
        self.assertEqual(out["n_eligible_cells_for_source_observations"],0)
        self.assertIsNone(out["descriptive_partial_rank_results"]["woody"]
                          ["rarefied20_expected_species_richness"]
                          ["gbif_and_crossing_opportunity"]["rho"])
    def test_no_spurious_effect_if_traffic_has_no_variation(self):
        flora,fra=fake_sources()
        for r in fra["by_original_50km_cells"].values():
            r["median_reported_total_daily_through_trains_at_eligible_crossings"]=5
        out=build_block("us_midwest",flora,fra,TRAITS)
        x=out["descriptive_partial_rank_results"]["all_GRIIS_introduced"]
        self.assertIsNone(x["raw_recorded_species_richness"]
                          ["gbif_and_crossing_opportunity"]["rho"])

    def test_cell_jackknife_does_not_create_country_replicates(self):
        flora,fra=fake_sources()
        out=build_block("us_midwest",flora,fra,TRAITS)
        test=out["leave_one_cell_out_sensitivity_gbif_and_crossing_adjusted"]
        self.assertEqual(test["woody"]["rarefied20_expected_species_richness"]
                         ["complete_n_withheld_cells"],16)
        self.assertTrue(test["woody"]["rarefied20_expected_species_richness"]
                        ["no_country_replication_or_statistical_significance_from_withheld_cells"])
        self.assertEqual(test["woody"]["raw_recorded_species_richness"]
                         ["baseline_partial_rho"],
                         out["descriptive_partial_rank_results"]["woody"]
                         ["raw_recorded_species_richness"]
                         ["gbif_and_crossing_opportunity"]["rho"])
    def test_no_claim_of_operating_train_gradient_if_positive_support_too_small(self):
        flora,fra=fake_sources()
        # Simulate a region with 8 running-report cells and 8 0-train cells.
        for index,k in enumerate(sorted(fra["by_original_50km_cells"])):
            if index>=8:
                fra["by_original_50km_cells"][k]["median_reported_total_daily_through_trains_at_eligible_crossings"]=0
        out=build_block("us_midwest",flora,fra,TRAITS)
        self.assertEqual(out["n_nonzero_reported_train_50km_cells"],8)
        stats=out["positive_reported_train_only_sensitivity"]["woody"]["rarefied20_expected_species_richness"]
        self.assertIsNone(stats["gbif_and_crossing_opportunity"]["rho"])
        self.assertFalse(stats["gbif_and_crossing_opportunity"]["support"])
if __name__=="__main__":unittest.main()
