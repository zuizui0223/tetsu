#!/usr/bin/env python3
"""v5.3 held-out US Southeast physical RAIL vs ROAD and introduced flora.

Protocol docs/holdout_us_southeast_v53_protocol.md was committed
BEFORE downloading its real mapped 50km physical networks or dataset
composition. The posthoc-chosen K80 recording-depth main diagnostic is
locked; new source data are untouched out-of-analysis-region geography,
not an independently randomized causal intervention.

Uses source-identical K20/40/80 fixed-K80-cell cohort; all-plants
GBIF source dataset overlap adjustment, and the 2018–2025 period with
nearest physical RAIL/ROAD connector in each original 2500 km2 cell.
Physical link geometry from 2026 not actual historical freight services.

Returns NOT ESTIMABLE if <8 source-qualified K80 cells, <20 mapped
RAIL/ROAD cell pairs, or source-rank degeneracy. No pairwise p value.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from analyze_joint_rarefaction_dataset_overlaps_v52 import (
    block_source,FIXED_RECORD_COHORT
)
from analyze_6country_50km_network_plants_v44 import PERIODS,STRATEGIES
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS

HOLDOUT_BLOCK="us_southeast"
PRIMARY_PERIOD="2018_2025"
PRIMARY_SOURCE_ANCHOR="nearest_in_cell"
PRIMARY_K="80"
PRIMARY_FIELD="dataset_overlap_adjusted_rho"
MIN_K80_CELLS=8
MIN_MAPPED_SOURCE_PAIRS=20


def prepare_original_holdout_sources(flora,rail,road,dataset):
    if HOLDOUT_BLOCK not in BLOCKS or BLOCKS[HOLDOUT_BLOCK][0]!="US":
        raise ValueError("Geographic source holdout country frame has changed")
    if not flora.get("complete") or flora.get("source_errors"):
        raise ValueError("Original eight-country GRIIS flora source invalid")
    if HOLDOUT_BLOCK not in flora.get("regions",{}) or flora["regions"][HOLDOUT_BLOCK].get("country")!="US":
        raise ValueError("Original US Southeast 50km GRIIS source unavailable")
    for name,v in (("rail",rail),("road",road)):
        if v.get("complete") is not True or set(v.get("regions",{}))!={HOLDOUT_BLOCK}:
            raise ValueError(f"Source-geographically untouched holdout {name} incomplete")
    if rail["regions"][HOLDOUT_BLOCK].get("country")!="US":
        raise ValueError("Holdout physical rail source country mismatch")
    if road["regions"][HOLDOUT_BLOCK].get("source",{}).get("country")!="US":
        raise ValueError("Holdout physical road source country mismatch")
    if dataset.get("complete") is not True or dataset.get("errors") or dataset.get("block")!=HOLDOUT_BLOCK or dataset.get("country")!="US":
        raise ValueError("New actual GBIF dataset source holdout missing or incomplete")
    if len(dataset.get("cells",{}))!=16:
        raise ValueError("Holdout 50km data-source cell coverage incomplete")
    datasets={"complete":True,"regions":{
        HOLDOUT_BLOCK:{"country":"US","cells":dataset["cells"],
                       "recording_source":dataset["source"]}
    }}
    return datasets


def evaluate(flora,rail,road,dataset):
    joined=prepare_original_holdout_sources(flora,rail,road,dataset)
    src=block_source(HOLDOUT_BLOCK,flora,rail,road,joined)
    value=src["periods"][PRIMARY_PERIOD][PRIMARY_SOURCE_ANCHOR][PRIMARY_K]
    n=src["n_cells_fixed_at_K80_both_periods"]
    pairs=value["n_nonindependent_source_cell_pairs"]
    coef=value["coefficients"][PRIMARY_FIELD]
    eligible=(n>=MIN_K80_CELLS and pairs>=MIN_MAPPED_SOURCE_PAIRS and coef is not None)
    result={
        "utc":datetime.now(timezone.utc).isoformat(),
        "protocol_frozen_before_US_Southeast_source_rail_road_effect_check":
            "docs/holdout_us_southeast_v53_protocol.md",
        "source_geography":HOLDOUT_BLOCK,
        "country":"US",
        "test_status":"ESTIMABLE_EXPLORATORY_REGION_DIAGNOSTIC" if eligible else "NOT_ESTIMABLE",
        "primary_precommitted_metric":{
            "period":PRIMARY_PERIOD,
            "rarefaction_indexed_RECORDS":FIXED_RECORD_COHORT,
            "source_anchor":PRIMARY_SOURCE_ANCHOR,
            "observed_introduced_specieskey_response":"ratio of expected rarefied20? NO; K80 expected intersection/expected union",
            "partial_rank_route_predictor":"log physical RAIL distance / physical ROAD distance",
            "source_covariates":"geography, ROAD shortest paths, GBIF record scale/balance, ALL-Plantae DATASET_KEY overlap",
            "negative_sign_compatible_with_railway_source_hypothesis":True,
            "n_valid_50km_source_cells":n,
            "n_nonindependent_plant_rail_road_source_pairs":pairs,
            "source_positive_support_min8_cells_min20_pairs":n>=MIN_K80_CELLS and pairs>=MIN_MAPPED_SOURCE_PAIRS,
            "identified_descriptive_partial_rho":coef if eligible else None,
            "direction":"negative_hypothesis_compatible" if eligible and coef<0 else (
                "positive_opposite_direction" if eligible and coef>0 else (
                    "zero" if eligible else "non_estimable")),
            "no_pairwise_pvalue_or_independent_country_replication_claim":True
        },
        "all_precommitted_sensitivities_and_source_gates":src,
        "source_is_previous_frozen_US_Southeast_GRIIS_and_new_2026_OvERTURE":True,
        "not_historical_active_train_freight_or_seed_transport_causality":True,
        "not_a_globally_independent_country_holdout":True,
        "limits":[
            "An outcome-unseen US Southeast GEOGRAPHIC region was selected after observing K80 posthoc six-country results; this is not independently preregistered research or a new country.",
            "New independent geographic block within the US shares GRIIS taxon status and other national recording systematics with the original US Midwest block.",
            "All plant data remain positive indexed GBIF taxon records and all-plant dataset source fractions, not DOI event rows or equal-effort vegetation surveys.",
            "Sample size is one 200km block; within-block species similarities share 50km cells and have no valid independent pair-level significance.",
            "Physical Overture 2026 rail and road graphs do not establish actual freight train operation in 2010-2025 or seed/spore transport.",
            "A missing data-source or network support gate is NOT a biological zero or contradictory plant effect.",
        ]
    }
    print(json.dumps({
        "holdout":HOLDOUT_BLOCK,
        "source_cells_K80":n,
        "source_pairs_K80":pairs,
        "status":result["test_status"],
        "primary_K80_dataset_overlap_adjusted_rho":result["primary_precommitted_metric"]["identified_descriptive_partial_rho"],
        "direction":result["primary_precommitted_metric"]["direction"],
        "K20_dataset_adjusted":src["periods"][PRIMARY_PERIOD][PRIMARY_SOURCE_ANCHOR]["20"]["coefficients"][PRIMARY_FIELD],
        "K40_dataset_adjusted":src["periods"][PRIMARY_PERIOD][PRIMARY_SOURCE_ANCHOR]["40"]["coefficients"][PRIMARY_FIELD]
    },ensure_ascii=False),flush=True)
    return result


def main():
    ap=argparse.ArgumentParser()
    for name in ("flora","rail","road","dataset"):
        ap.add_argument("--"+name,required=True)
    ap.add_argument("--output",default="results/audit/holdout_US_southeast_GRIIS_K80_route_v53.json")
    args=ap.parse_args()
    files=[json.loads(Path(getattr(args,k)).read_text(encoding="utf-8"))
           for k in ("flora","rail","road","dataset")]
    result=evaluate(*files)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"completed":True,"source_support":result["test_status"],
                      "no_caused_seed_transport":True}),flush=True)


if __name__=="__main__":
    main()
