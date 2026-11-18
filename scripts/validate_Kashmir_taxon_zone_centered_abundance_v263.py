#!/usr/bin/env python3
"""v26.3: stricter authentic Kashmir species identity × habitat controls.

Freeze same 2014 presence-risk set, original 2017 re-recording target,
same 27 whole-spatial-site source heldout folds as v26.1.

M2Z: source group×zone intercepts + TAXON×ZONE categorical indicators
M3Z: M2Z + common within-2014-plot numeric abundance midrank
M4Z: M2Z + within-2014-plot rank CENTERED by each source taxon×zone's
     training-sites-only mean, with never-seen (taxon,zone) residual set
     to zero rather than imputing actual source test outcomes.

Taxon×zone modeling is more flexible and sparse than the original taxon
model, so no claim that one model is biologically correct. This sensitivity
explicitly probes whether v26.1 added abundance prediction was merely
picking up unmodeled taxon-specific preference for particular habitats.
Never create per-site, per-taxon, or GPS-containing outputs.

No observation of real individual survival, seed transport, grass cutting,
train traffic, wind exposure, or genetic adaptation. Two source years of ONE
non-random longitudinal field survey are not independent research studies.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean

import numpy as np

from validate_Kashmir_abundance_after_species_site_holdout_v261 import (
    GROUPS,ZONES,N_TAXA,MODELS,BASE_RIDGE_C,
    load_authentic_source,make_design,fit_predict,row_losses,
    fold_evaluation
)

LABELS=("M2Z_species_by_zone",
        "M3Z_species_by_zone_plus_rank",
        "M4Z_species_by_zone_plus_training_species_zone_centered_rank")


def build_species_zone_design(rows:list[dict])->tuple[np.ndarray,np.ndarray,np.ndarray,dict]:
    y,base,by_site=make_design(rows)
    n=len(rows)
    gz=base[MODELS[0]]
    ztax=np.zeros((n,N_TAXA*len(ZONES)),dtype=float)
    ranks=np.zeros(n,dtype=float)
    species_by_zone=[]
    for i,r in enumerate(rows):
        z=ZONES.index(r["zone"])
        t=r["taxon"]
        ztax[i,z*N_TAXA+t]=1.
        ranks[i]=r["x2014_rank"]
        species_by_zone.append((t,z))
    return y,np.column_stack((gz,ztax)),ranks,{"by_site":by_site,"keys":species_by_zone}


def train_only_centered_rank(
    train_indices:np.ndarray,test_indices:np.ndarray,
    keys:list[tuple[int,int]],ranks:np.ndarray
)->tuple[np.ndarray,np.ndarray,dict]:
    means=defaultdict(list)
    for i in train_indices:
        means[keys[i]].append(ranks[i])
    train_av={k:mean(v) for k,v in means.items()}
    center_train=np.array([ranks[i]-train_av[keys[i]] for i in train_indices],dtype=float)
    center_test=np.array(
        [ranks[i]-train_av[keys[i]] if keys[i] in train_av else 0.
         for i in test_indices],dtype=float)
    support={
        "n_train_source_taxon_zone_categories":len(train_av),
        "n_test_events_with_train_same_taxon_zone":sum(keys[i] in train_av for i in test_indices),
        "n_test_events_without_train_same_taxon_zone":sum(keys[i] not in train_av for i in test_indices),
        "test_label_or_future_record_never_used_for_centering":True
    }
    if abs(np.mean(center_train))>1e-10:
        raise AssertionError("Source train-only species-zone rank centering failed")
    return center_train,center_test,support


def audit(rows:list[dict],C:float=BASE_RIDGE_C)->dict:
    y,base,ranks,source=build_species_zone_design(rows)
    n=len(rows);bysite=source["by_site"];keys=source["keys"]
    models={m:np.zeros(n,dtype=float) for m in LABELS}
    support={site:None for site in bysite}
    for site in sorted(bysite):
        test=np.array(bysite[site],dtype=int)
        mask=np.ones(n,dtype=bool);mask[test]=False
        train=np.flatnonzero(mask)
        if any(rows[i]["site"]==site for i in train):
            raise AssertionError("Entire original site must be excluded from training")
        tcenter,vcenter,group_support=train_only_centered_rank(train,test,keys,ranks)
        support[site]=group_support
        Xs={
            LABELS[0]:(base[train],base[test]),
            LABELS[1]:(
                np.column_stack((base[train],ranks[train])),
                np.column_stack((base[test],ranks[test]))
            ),
            LABELS[2]:(
                np.column_stack((base[train],tcenter)),
                np.column_stack((base[test],vcenter))
            )
        }
        for key,(Xtrain,Xtest) in Xs.items():
            p,_=fit_predict(Xtrain,y[train],Xtest,C)
            models[key][test]=p
    if any(np.any(models[m]<=0) or np.any(models[m]>=1) for m in LABELS):
        raise RuntimeError("Not every original spatial-site row has heldout probability")
    scores={m:row_losses(y,models[m]) for m in LABELS}
    per_group={}
    for group in GROUPS:
        per_group[group]={}
        for zone in ZONES:
            by=[]
            for site,ids in bysite.items():
                if rows[ids[0]]["group"]!=group:continue
                ix=[i for i in ids if rows[i]["zone"]==zone]
                if not ix:raise ValueError("Original strict site missing habitat")
                by.append({
                    "n_events":len(ix),
                    "n_rerecorded":int(y[ix].sum()),
                    "n_test_source_taxon_zone_seen_in_other_site":sum(
                        keys[i] in set(keys[j] for j in range(n) if rows[j]["site"]!=site) for i in ix),
                    "ll":{m:float(np.mean(scores[m][0][ix])) for m in LABELS},
                    "brier":{m:float(np.mean(scores[m][1][ix])) for m in LABELS}
                })
            if len(by)!=(15 if group=="station" else 12):
                raise ValueError("Cohort expected 15 or 12 independent source locations")
            avg={m:{"source_equal_site_logloss":mean(x["ll"][m] for x in by),
                    "source_equal_site_Brier":mean(x["brier"][m] for x in by)}
                 for m in LABELS}
            per_group[group][zone]={
                "source_site_heldout_folds":len(by),
                "2014_risk_species_site_events":sum(x["n_events"] for x in by),
                "2017_rerecorded_species_site_events":sum(x["n_rerecorded"] for x in by),
                "test_species_zone_known_on_other_training_sites_events":sum(x["n_test_source_taxon_zone_seen_in_other_site"] for x in by),
                "models":avg,
                "M2Z_minus_M3Z_logloss_gain":avg[LABELS[0]]["source_equal_site_logloss"]-avg[LABELS[1]]["source_equal_site_logloss"],
                "M2Z_minus_M4Z_centered_logloss_gain":avg[LABELS[0]]["source_equal_site_logloss"]-avg[LABELS[2]]["source_equal_site_logloss"],
                "M2Z_minus_M3Z_Brier_gain":avg[LABELS[0]]["source_equal_site_Brier"]-avg[LABELS[1]]["source_equal_site_Brier"],
                "M2Z_minus_M4Z_centered_Brier_gain":avg[LABELS[0]]["source_equal_site_Brier"]-avg[LABELS[2]]["source_equal_site_Brier"],
                "n_site_folds_M3Z_better_than_M2Z":sum(x["ll"][LABELS[1]]<x["ll"][LABELS[0]] for x in by),
                "n_site_folds_M4Z_better_than_M2Z":sum(x["ll"][LABELS[2]]<x["ll"][LABELS[0]] for x in by),
                "folds_not_independent_biological_experiments":True,
            }
    return {
        "same_27_original_spatial_sites_source_holdout":True,
        "regularization_C_fixed":C,
        "model_names":LABELS,
        "source_species_zone_feature_count":N_TAXA*len(ZONES),
        "source_taxon_identity_not_genetic_individual_identity":True,
        "source_group_zone_habitat_results":per_group,
        "total_training_species_zone_support_by_original_test_site":{
            "min":min(v["n_train_source_taxon_zone_categories"] for v in support.values()),
            "median":float(np.median([v["n_train_source_taxon_zone_categories"] for v in support.values()])),
            "max":max(v["n_train_source_taxon_zone_categories"] for v in support.values()),
        },
        "n_original_holdout_source_events_unseen_taxon_zone_by_training":sum(
            v["n_test_events_without_train_same_taxon_zone"] for v in support.values()),
        "training_only_rank_centering_never_leaks_test_site":True,
        "no_strict_causal_or_fitness_inference":True
    }


def run(source:Path,C:float=BASE_RIDGE_C)->dict:
    rows,source_qc=load_authentic_source(source)
    result=audit(rows,C)
    return {
        "version":"v26.3",
        "source":source_qc,
        "result":result,
        "status":"AUTHENTIC_ORIGINAL_SOURCE_IF_CI_SUCCESS",
        "one_original_repeated_Kashmir_field_study_only":True,
        "not_outcome_tuned_regularization":True,
        "no_original_taxon_by_site_rows_disclosed":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    out=run(a.source)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    for group in GROUPS:
        for zone in ZONES:
            r=out["result"]["source_group_zone_habitat_results"][group][zone]
            print(json.dumps({
                "version":"v26.3","original_site_group":group,"source_habitat":zone,
                "n_sites":r["source_site_heldout_folds"],
                "source_risk_events":r["2014_risk_species_site_events"],
                "M2Z_M3Z_logloss_gain":r["M2Z_minus_M3Z_logloss_gain"],
                "M2Z_M4Z_centered_logloss_gain":r["M2Z_minus_M4Z_centered_logloss_gain"],
                "M2Z_M3Z_Brier_gain":r["M2Z_minus_M3Z_Brier_gain"],
                "M2Z_M4Z_centered_Brier_gain":r["M2Z_minus_M4Z_centered_Brier_gain"],
                "n_site_folds_M3Z_improves":r["n_site_folds_M3Z_better_than_M2Z"]
            },ensure_ascii=False))
if __name__=="__main__":
    main()
