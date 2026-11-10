"""Frozen geographically untouched US Southeast K80 source-holdout support gates."""
from __future__ import annotations
import json
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[0]))
from test_dataset_overlap_adjusted_rail_plants_v51 import fake_sources
from evaluate_us_southeast_geographic_holdout_v53 import (
    evaluate,prepare_original_holdout_sources,HOLDOUT_BLOCK,PRIMARY_FIELD,
    PRIMARY_K,MIN_K80_CELLS,MIN_MAPPED_SOURCE_PAIRS
)


def mock_holdout():
    flora,rails,roads,chunks=fake_sources()
    base="us_midwest"
    def transform(x):
        return json.loads(json.dumps(x).replace(base,HOLDOUT_BLOCK))
    f=transform(flora)
    rail={"complete":True,"regions":{HOLDOUT_BLOCK:frozen}
          for frozen in []}
    rail={"complete":True,"regions":{
        HOLDOUT_BLOCK:transform(rails["regions"][base])
    }}
    road={"complete":True,"regions":{
        HOLDOUT_BLOCK:transform(roads["regions"][base])
    }}
    dataset=transform(chunks[0])
    dataset["block"]=HOLDOUT_BLOCK
    return f,rail,road,dataset


class HeldoutCountryRegionTests(unittest.TestCase):
    def test_predeclared_location_response_and_source_gates(self):
        self.assertEqual(HOLDOUT_BLOCK,"us_southeast")
        self.assertEqual(PRIMARY_K,"80")
        self.assertEqual(PRIMARY_FIELD,"dataset_overlap_adjusted_rho")
        self.assertEqual(MIN_K80_CELLS,8)
        self.assertEqual(MIN_MAPPED_SOURCE_PAIRS,20)

    def test_same_country_different_geography_does_not_invent_independent_country(self):
        f,r,d,ds=mock_holdout()
        out=evaluate(f,r,d,ds)
        self.assertEqual(out["source_geography"],HOLDOUT_BLOCK)
        self.assertEqual(out["country"],"US")
        self.assertTrue(out["not_a_globally_independent_country_holdout"])
        self.assertTrue(out["not_historical_active_train_freight_or_seed_transport_causality"])
        self.assertEqual(out["primary_precommitted_metric"]["rarefaction_indexed_RECORDS"],80)
        self.assertEqual(out["primary_precommitted_metric"]["n_valid_50km_source_cells"],16)
        self.assertTrue(out["primary_precommitted_metric"]["no_pairwise_pvalue_or_independent_country_replication_claim"])

    def test_source_quality_ineligible_block_not_falsely_zero(self):
        f,r,d,ds=mock_holdout()
        cells=f["regions"][HOLDOUT_BLOCK]["cells"]
        for idx,cid in enumerate(sorted(cells)):
            if idx>=7:
                cells[cid]["source_pre_outcome_both_period_recording_eligible"]=False
        out=evaluate(f,r,d,ds)
        self.assertEqual(out["test_status"],"NOT_ESTIMABLE")
        self.assertEqual(out["primary_precommitted_metric"]["n_valid_50km_source_cells"],7)
        self.assertIsNone(out["primary_precommitted_metric"]["identified_descriptive_partial_rho"])
        self.assertEqual(out["primary_precommitted_metric"]["direction"],"non_estimable")

    def test_missing_GBIF_or_actual_road_source_fails_closed(self):
        f,r,d,ds=mock_holdout()
        d["regions"].pop(HOLDOUT_BLOCK)
        with self.assertRaises(ValueError):
            evaluate(f,r,d,ds)
        f,r,d,ds=mock_holdout()
        ds["complete"]=False
        with self.assertRaises(ValueError):
            prepare_original_holdout_sources(f,r,d,ds)

    def test_new_source_country_mismatch_is_rejected(self):
        f,r,d,ds=mock_holdout()
        ds["country"]="FR"
        with self.assertRaises(ValueError):
            prepare_original_holdout_sources(f,r,d,ds)


if __name__=="__main__":
    unittest.main()
