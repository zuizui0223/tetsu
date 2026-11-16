"""Six-country fixed K80 same-cell dataset source + rarefaction numerical QA."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[0]))
from analyze_joint_rarefaction_dataset_overlaps_v52 import (
    block_source,aggregate,run,DEPTHS,FIXED_RECORD_COHORT,PERIODS,STRATEGIES
)
from analyze_dataset_overlap_adjusted_rail_plants_v51 import merge_dataset_chunks
from test_dataset_overlap_adjusted_rail_plants_v51 import fake_sources


class JointRecordingSensitivityTests(unittest.TestCase):
    def test_exact_fixed_K80_sample_country_source(self):
        f,r,d,raw=fake_sources()
        source=merge_dataset_chunks(raw)
        block=block_source("fr_burgundy",f,r,d,source)
        self.assertEqual(block["n_cells_fixed_at_K80_both_periods"],16)
        self.assertEqual(FIXED_RECORD_COHORT,80)
        for period in PERIODS:
            for strategy in STRATEGIES:
                out=block["periods"][period][strategy]
                self.assertEqual(set(out),{"20","40","80"})
                self.assertEqual(len({out[str(k)]["n_nonindependent_source_cell_pairs"] for k in DEPTHS}),1)
                self.assertEqual(out["20"]["n_nonindependent_source_cell_pairs"],120)
                self.assertIn("dataset_overlap_adjusted_rho",out["20"]["coefficients"])
                self.assertTrue(out["80"]["no_railway_seed_transport_or_pairwise_p_value"])

    def test_insufficient_K80_source_cell_excluded_without_setting_zero_effect(self):
        f,r,d,raw=fake_sources()
        block="us_midwest";cid=f"{block}_x0_y0"
        for period in PERIODS:
            cells=f["regions"][block]["cells"][cid]["periods"][period]
            cells["observed_alien_specieskey_counts"]={str(i):10 for i in range(5)}
            cells["n_indexed_GRIIS_introduced_plant_records"]=50
        src=merge_dataset_chunks(raw)
        outcome=block_source(block,f,r,d,src)
        self.assertEqual(outcome["n_cells_fixed_at_K80_both_periods"],15)
        self.assertEqual(outcome["periods"]["2018_2025"]["nearest_in_cell"]["20"]["n_nonindependent_source_cell_pairs"],105)

    def test_country_equivalent_weight_no_pairwise_pvalue(self):
        f,r,d,raw=fake_sources()
        result=run(f,r,d,merge_dataset_chunks(raw))
        self.assertEqual(len(result["results"]),6)
        self.assertEqual(result["n_frozen_country_blocks"],6)
        self.assertTrue(result["posthoc_noncausal_method_sensitivity"])
        for period in PERIODS:
            for depth in DEPTHS:
                a=result["equal_country_blocks"][period]["nearest_in_cell"][str(depth)]
                self.assertIn("dataset_overlap_adjusted_rho",a)
                self.assertTrue(a["dataset_overlap_adjusted_rho"]["not_cell_pair_p_value_or_causal_train_seed_transfer"])
        self.assertIn("dataset-source",result["limits"][0].lower())

    def test_bad_dataset_country_source_not_silently_filled(self):
        f,r,d,raw=fake_sources()
        ds=merge_dataset_chunks(raw)
        ds["regions"].pop("de_saxony")
        with self.assertRaisesRegex(ValueError,"Six-country source geography changed"):
            run(f,r,d,ds)
        ds=merge_dataset_chunks(raw)
        ds["complete"]=False
        with self.assertRaisesRegex(ValueError,"One original source not complete"):
            run(f,r,d,ds)


if __name__=="__main__":
    unittest.main()
