#!/usr/bin/env python3
"""v26.2: prespecified ridge regularization sensitivity for source-verified
v26.1 Kashmir taxon identity plus original 2014 within-plot abundance rank.

Compare fixed C=0.2, 1.0, 5.0 with SAME original 27 site leave-one-site-out
folds. N source studies remains 1. Not model/hyperparameter selection:
never choose a C after observing test performance and claim a new p-value.

All original publisher SHA verification and temporal site strict source
cohorts inherited from v26.1. Raw site x taxon and original abundance
values never serialized.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import mean

from validate_Kashmir_abundance_after_species_site_holdout_v261 import (
    BASE_RIDGE_C, MODELS, load_authentic_source, fold_evaluation
)

C_VALUES=(0.2,1.0,5.0)
EXPECTED_SOURCE_SHA="529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474"
MAX_SOLVER_SCORE_NUMERIC_DIFFERENCE=1e-8  # liblinear optimizer stopping variation, never source drift
EXPECTED_MAIN_TRACK_M2_M3={
    "station":0.015090379941841148,
    "between":0.040527147537008146,
}


def compare(fitted:dict[float,dict]) -> dict:
    if set(fitted)!=set(C_VALUES):
        raise ValueError("All three previously fixed shrinkage settings required")
    out={}
    for group in ("station","between"):
        zone="Track"
        rows={C:fitted[C]["source_predictions_by_group_and_habitat"][group][zone]
              for C in C_VALUES}
        out[group]={
            "source_independent_sites":rows[1.0]["n_original_source_heldout_site_folds"],
            "source_2014_risk_events":rows[1.0]["source_2014_at_risk_taxon_site_events"],
            "source_2017_rerecorded":rows[1.0]["source_2017_re_recorded_original_taxon_site_events"],
            "rank_added_after_species_logloss_gain_by_fixed_C":{
                str(C):rows[C]["M2_minus_M3_rank_added_after_taxon_logloss_gain"]
                for C in C_VALUES},
            "rank_and_track_interaction_added_after_species_logloss_gain_by_fixed_C":{
                str(C):rows[C]["M2_minus_M4_rank_and_track_interaction_added_after_taxon_logloss_gain"]
                for C in C_VALUES},
            "rank_added_after_species_Brier_gain_by_fixed_C":{
                str(C):rows[C]["M2_minus_M3_rank_added_after_taxon_Brier_gain"]
                for C in C_VALUES},
            "sign_consistent_in_all_three_fixed_C":{
                metric:all(rows[C][metric]>0 for C in C_VALUES)
                for metric in ("M2_minus_M3_rank_added_after_taxon_logloss_gain",
                               "M2_minus_M3_rank_added_after_taxon_Brier_gain")
            },
            "source_site_folds_where_M3_beats_M2_by_C":{
                str(C):rows[C]["n_source_site_folds_M3_improves_M2"]
                for C in C_VALUES},
        }
        if abs(rows[1.0]["M2_minus_M3_rank_added_after_taxon_logloss_gain"]-
               EXPECTED_MAIN_TRACK_M2_M3[group])>MAX_SOLVER_SCORE_NUMERIC_DIFFERENCE:
            raise ValueError(f"Frozen v26.1 C=1 group={group} Track score did not reproduce: got={rows[1.0]['M2_minus_M3_rank_added_after_taxon_logloss_gain']:.15g} expected={EXPECTED_MAIN_TRACK_M2_M3[group]:.15g}")
    return out


def run(path:Path)->dict:
    rows,qc=load_authentic_source(path)
    if qc["original_Dryad_publisher_SHA256"]!=EXPECTED_SOURCE_SHA:
        raise ValueError("Original Dryad publisher source SHA drift")
    scores={C:fold_evaluation(rows,C) for C in C_VALUES}
    main=compare(scores)
    other={}
    for group in ("station","between"):
        other[group]={}
        for z in ("Embankment","Interior Plot 1","Interior Plot 2"):
            other[group][z]={
                "rank_after_species_gain_by_fixed_C":{
                    str(C):scores[C]["source_predictions_by_group_and_habitat"][group][z][
                        "M2_minus_M3_rank_added_after_taxon_logloss_gain"]
                    for C in C_VALUES}
            }
    return {
        "version":"v26.2",
        "source_provenance":qc,
        "all_fixed_C_values":list(C_VALUES),
        "M2_species_only_vs_M3_species_plus_rank_track_sensitivity":main,
        "adjacent_source_habitats_not_independent_railway_studies":other,
        "source_ridge_C1_v261_score_reproduced_within_numeric_optimizer_tolerance_1e_minus_8":True,
        "source_ridge_C1_numeric_tolerance_nats":MAX_SOLVER_SCORE_NUMERIC_DIFFERENCE,
        "no_selection_of_C_by_heldout_outcomes":True,
        "no_causal_survival_mowing_train_dispersal_or_adaptation_conclusions":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out=run(a.source)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
        "version":"v26.2","primary_track_ridge_C_sensitivity":
            out["M2_species_only_vs_M3_species_plus_rank_track_sensitivity"]
    },ensure_ascii=False))


if __name__=="__main__":
    main()
