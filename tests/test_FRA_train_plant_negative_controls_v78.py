"""v7.8 synthetic guardrails for current FRA/GBIF/GIFT negative controls."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_FRA_train_plant_negative_controls_v78 import (
    BLOCKS,PERIOD,CONTROL_OUTCOMES,EXPOSURES,source_rows,source_rho,one_block)
from analyze_US_train_counts_GRIIS_traits_v76 import build_block,COVARIATES

TRAITS={"11":"herbaceous","12":"herbaceous","13":"herbaceous",
        "21":"woody","22":"woody","23":"woody"}

def fixture():
    flora={"country":"US","cells":{}}
    fra={"source_open_closed_balance_verified":True,
         "by_original_50km_cells":{}}
    for x in range(4):
        for y in range(4):
            cid=f"us_midwest_x{x}_y{y}"
            counts={"11":25,"12":20,"13":15,"21":24,"22":18,"23":14}
            flora["cells"][cid]={
                "source_pre_outcome_both_period_recording_eligible":True,
                "periods":{PERIOD:{
                    "full_facet_source_pagination_completed":True,
                    "alien_observed_minimum_5_species_and_20_records":True,
                    "indexed_all_plantae_source_records":1000+31*x+17*y+4*x*y,
                    "n_all_observed_species_key_categories":40+x+y,
                    "observed_alien_specieskey_counts":counts}}}
            fra["by_original_50km_cells"][cid]={
                "train_opportunity_metrics_USE_confirmed_OPEN_crossings_ONLY":True,
                "n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields":10,
                "n_crossings_confirmed_open_in_current_source":10+x+y,
                "n_crossings_confirmed_closed_excluded_from_operating_exposure":4+x,
                "n_distinct_crossing_inventory_records":14+2*x+y,
                "median_reported_total_daily_through_trains_at_eligible_crossings":2+x+3*y+(x*y)%3,
                "reported_train_count_year_histogram":{"2023":10,"2026":2}}
    return flora,fra

class SourceFalsificationControls(unittest.TestCase):
    def test_two_preselected_us_macro_blocks_not_six_fake_countries(self):
        self.assertEqual(BLOCKS,("us_midwest","us_southeast"))
        self.assertIn("non_GRIIS_speciesKey_categories_UNKNOWN_nativity",CONTROL_OUTCOMES)
        self.assertNotIn("native_species_count",CONTROL_OUTCOMES)
    def test_same_v76_eligible_cell_cohort_and_unknown_native_status(self):
        flora,fra=fixture()
        rows,exc=source_rows("us_midwest",flora,fra,TRAITS)
        self.assertEqual(len(rows),16)
        self.assertEqual(exc,{})
        self.assertEqual(rows[0]["responses"]["non_GRIIS_speciesKey_categories_UNKNOWN_nativity"],34)
        self.assertEqual(len(EXPOSURES),2)
    def test_closed_fraction_not_called_train_count(self):
        flora,fra=fixture()
        rows,exc=source_rows("us_midwest",flora,fra,TRAITS)
        self.assertTrue(all(0<=r["_closed_crossing_fraction"]<=1 for r in rows))
        self.assertNotEqual(rows[0]["_closed_crossing_fraction"],
                            rows[0]["reported_train_exposure_log1p_median_day_night_through_trains"])
    def test_source_plant_count_comparator_does_not_invent_species(self):
        flora,fra=fixture()
        cid="us_midwest_x0_y0"
        flora["cells"][cid]["periods"][PERIOD]["n_all_observed_species_key_categories"]=2
        with self.assertRaisesRegex(ValueError,"inconsistent"):
            source_rows("us_midwest",flora,fra,TRAITS)
    def test_actual_prior_coefficient_requires_exact_reproduction(self):
        flora,fra=fixture()
        prev={"regions":{"us_midwest":build_block("us_midwest",flora,fra,TRAITS)}}
        v=one_block("us_midwest",flora,fra,TRAITS,prev)
        same=prev["regions"]["us_midwest"]["descriptive_partial_rank_results"][
            "woody"]["rarefied20_expected_species_richness"]["gbif_and_crossing_opportunity"]["rho"]
        self.assertEqual(v["reproduced_original_v76_woody_rho"],same)
        self.assertEqual(v["n_same_eligible_50km_source_cells_as_v76"],16)
        self.assertEqual(v["positive_through_train_only_same_source_cell_sensitivity"][
            "woody_GRIIS_rarefied20_expected_richness"]["n_50km_cells"],16)
    def test_changed_baseline_fails_closed(self):
        flora,fra=fixture()
        prev={"regions":{"us_midwest":build_block("us_midwest",flora,fra,TRAITS)}}
        prev["regions"]["us_midwest"]["descriptive_partial_rank_results"]["woody"][
            "rarefied20_expected_species_richness"]["gbif_and_crossing_opportunity"]["rho"]=0.999
        with self.assertRaisesRegex(ValueError,"did not replicate exactly"):
            one_block("us_midwest",flora,fra,TRAITS,prev)
    def test_unmeasured_closed_source_does_not_enter_traffic(self):
        flora,fra=fixture()
        fra["by_original_50km_cells"]["us_midwest_x0_y0"][
            "train_opportunity_metrics_USE_confirmed_OPEN_crossings_ONLY"]=False
        with self.assertRaisesRegex(ValueError,"closed"):
            source_rows("us_midwest",flora,fra,TRAITS)

if __name__=="__main__":unittest.main()
