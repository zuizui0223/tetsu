#!/usr/bin/env python3
"""v26.7 original MIREN ROAD: leave an entire biological REGION code out.

Independent generalization stress test of v26.6 initial abundance rank ->
later *re-recording*. MIREN mountain ROADS, NOT railway/ballast observations.

Reuses exact author-original 2012+ first-survey 251 spatial sites, 20 roads,
7 Region codes, 14,682 early taxon×plot events from the source-authenticated
v26.6 reader; verifies the original Zenodo CSV MD5. Each heldout REGION
removes all its roads, sites, years, plots, species observations and outcomes
from training. Pooled L2 logistic models are fixed C=1, no tuning.

Because REGION is wholly unseen, unlike v26.6 road-holdout, the model's
features deliberately OMIT a region-indicator/intercept. All predictors are
baseline plot, initial survey year, year interval, species identity, taxon×
plot identity and within-plot original Abundance 1–3 midrank. These are
baseline source *metadata* not estimated survival or true arrival.

M2 = plot/year/gap/taxon identity.
M3 = M2 + initial within-plot ordinal abundance rank.
M2Z = plot/year/gap/taxon×plot identity.
M3Z = M2Z + initial within-plot abundance rank.
M4Z = M2Z + taxon×plot rank deviation, mean estimated only on training
regions; never use the heldout region in rank centering.

Evaluate original-site-equal scores, then road-equal, then REGION-code-equal,
not taxon-events-as-independent-replicates. Also separate taxon seen vs
unseen among other *training regions*, and taxon×habitat seen vs unseen.
A test-only taxon has zero training coefficient in L2 logistic; this is
out-of-region transfer, not a novel gene or tested dispersal.

Source MIREN repeated protocol versions differ, vegetation detection
imperfect, and even a positive result implies neither train transport nor
rail-specific selection. Raw site×taxon matrix/GPS never exported.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean

import numpy as np
from scipy.sparse import csr_matrix, hstack
from sklearn.feature_extraction import DictVectorizer

from audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 import source_read, evaluate_pairs
from analyze_MIREN_road_abundance_after_taxon_out_of_road_v266 import (
    risk_set, centered_rank_for_fold, fit_prediction
)

LABELS=("M2_taxon", "M3_taxon_rank", "M2Z_taxon_zone",
        "M3Z_taxon_zone_rank", "M4Z_taxon_zone_training_centered_rank")
PLOTS=("1","3")
C=1.0
EXPECTED={"sites":251,"roads":20,"regions":7,"events":14682,"taxa":1925}
PRIMARY=("M2_to_M3","M2Z_to_M3Z","M2Z_to_M4Z")


def features_without_region(rows:list[dict])->tuple[dict,np.ndarray,dict,dict,list[tuple[str,str]]]:
    """Fixed metadata feature dictionary; no REGION leakage into predictions."""
    base=DictVectorizer(sparse=True)
    species=DictVectorizer(sparse=True)
    taxonzone=DictVectorizer(sparse=True)
    B=base.fit_transform([{
        "plot="+r["plot"]:1.0,
        "first_year="+str(r["initial_year"]):1.0,
        "year_gap="+str(r["year_gap"]):1.0
    } for r in rows])
    T=species.fit_transform([{"taxon="+r["taxon"]:1.0} for r in rows])
    TZ=taxonzone.fit_transform([{"taxon_plot="+r["taxon"]+"|"+r["plot"]:1.0} for r in rows])
    rank=np.asarray([r["xrank"] for r in rows],dtype=float)
    R=csr_matrix(rank.reshape(-1,1))
    design={
        LABELS[0]:hstack([B,T],format="csr"),
        LABELS[1]:hstack([B,T,R],format="csr"),
        LABELS[2]:hstack([B,TZ],format="csr"),
        LABELS[3]:hstack([B,TZ,R],format="csr")
    }
    y=np.asarray([r["re_recorded"] for r in rows],dtype=int)
    byregion=defaultdict(list)
    keys=[]
    for i,r in enumerate(rows):
        byregion[r["region"]].append(i)
        keys.append((r["taxon"],r["plot"]))
    return design,y,dict(byregion),{
        "base_plot_year_interval_feature_count":B.shape[1],
        "source_distinct_species":T.shape[1],
        "source_species_by_plot_categories":TZ.shape[1],
        "no_region_or_road_one_hot_features":True,
        "test_only_species_columns_have_zero_weight_after_fit":True
    },keys


def score_one_model(y:np.ndarray,p:np.ndarray)->tuple[np.ndarray,np.ndarray]:
    if len(y)!=len(p) or np.any(p<=0) or np.any(p>=1):
        raise ValueError("Original source heldout probability invalid")
    return -(y*np.log(p)+(1-y)*np.log1p(-p)),(y-p)**2


def averages_from_rows(rows:list[dict], test:list[int], ll:dict, brier:dict,
                       mask:list[bool], region:str, plot:str)->dict:
    """Equal original site, then road, before any Region-level aggregation."""
    grouped=defaultdict(list)
    for local_i,global_i in enumerate(test):
        r=rows[global_i]
        if r["region"]!=region:
            raise AssertionError("Heldout region contained training-source observations")
        if r["plot"]==plot and mask[local_i]:
            grouped[(r["road"],r["site"])].append(local_i)
    byroad=defaultdict(list)
    for (road,site),ix in grouped.items():
        byroad[road].append({
            "site":site,
            "n_events":len(ix),
            "n_rerecorded":int(sum(rows[test[k]]["re_recorded"] for k in ix)),
            "logloss":{m:float(np.mean(ll[m][ix])) for m in LABELS},
            "brier":{m:float(np.mean(brier[m][ix])) for m in LABELS}
        })
    if not byroad:
        return {"identified":False,"original_sites":0,"source_taxon_site_events":0}
    one_road=[]
    for road,sites in sorted(byroad.items()):
        one_road.append({
            "source_sites":len(sites),
            "source_events":sum(x["n_events"] for x in sites),
            "scores":{metric:{m:mean(s[metric][m] for s in sites) for m in LABELS}
                      for metric in ("logloss","brier")}
        })
    out={
        "identified":True,
        "source_roads":len(one_road),
        "source_sites":sum(x["source_sites"] for x in one_road),
        "source_taxon_site_events":sum(x["source_events"] for x in one_road),
        "models":{
            metric:{m:mean(road["scores"][metric][m] for road in one_road)
                    for m in LABELS}
            for metric in ("logloss","brier")
        }
    }
    out["gain"]={
        metric:{
            "M2_to_M3":out["models"][metric][LABELS[0]]-out["models"][metric][LABELS[1]],
            "M2Z_to_M3Z":out["models"][metric][LABELS[2]]-out["models"][metric][LABELS[3]],
            "M2Z_to_M4Z":out["models"][metric][LABELS[2]]-out["models"][metric][LABELS[4]]
        } for metric in ("logloss","brier")
    }
    return out


def assess(rows:list[dict],cohort:dict|None=None)->dict:
    X,y,byregion,feature_meta,keys=features_without_region(rows)
    region_names=sorted(byregion)
    if len(region_names)<3:raise ValueError("Require >=3 source regions to identify region holdout")
    all_rank=np.asarray([r["xrank"] for r in rows])
    N=len(rows)
    outcomes={}
    for region in region_names:
        test=np.asarray(byregion[region],dtype=int)
        train_mask=np.ones(N,dtype=bool)
        train_mask[test]=False
        train=np.flatnonzero(train_mask)
        if any(rows[i]["region"]==region for i in train):
            raise AssertionError("Leakage of heldout region into source training data")
        train_species={rows[i]["taxon"] for i in train}
        train_zones={keys[i] for i in train}
        cent_train,cent_test=centered_rank_for_fold(all_rank,keys,train,test)
        predictions={}
        for name in LABELS[:4]:
            predictions[name]=fit_prediction(X[name][train],X[name][test],y[train],C)
        predictions[LABELS[4]]=fit_prediction(
            hstack([X[LABELS[2]][train],csr_matrix(cent_train.reshape(-1,1))],format="csr"),
            hstack([X[LABELS[2]][test],csr_matrix(cent_test.reshape(-1,1))],format="csr"),
            y[train],C
        )
        loss={};brier={}
        for name in LABELS:
            loss[name],brier[name]=score_one_model(y[test],predictions[name])
        groups={
            "all":[True]*len(test),
            "known_taxon":[rows[i]["taxon"] in train_species for i in test],
            "unseen_taxon":[rows[i]["taxon"] not in train_species for i in test],
            "known_taxon_zone":[keys[i] in train_zones for i in test],
            "unseen_taxon_zone":[keys[i] not in train_zones for i in test]
        }
        views={}
        for groupname,mask in groups.items():
            views[groupname]={
                plot:averages_from_rows(rows,test.tolist(),loss,brier,mask,region,plot)
                for plot in PLOTS
            }
        for plot in PLOTS:
            if not views["all"][plot]["identified"]:
                raise ValueError("Source region holdout lacks a complete paired plot")
        outcomes[region]={
            "source_region_code":region,
            "original_roads":views["all"]["1"]["source_roads"],
            "original_spatial_sites":views["all"]["1"]["source_sites"],
            "source_baseline_taxon_plot_events":
                views["all"]["1"]["source_taxon_site_events"]+
                views["all"]["3"]["source_taxon_site_events"],
            "heldout_taxon_species_seen_on_other_regions_events":sum(groups["known_taxon"]),
            "heldout_taxon_species_unseen_elsewhere_events":sum(groups["unseen_taxon"]),
            "heldout_taxon_zone_seen_on_other_regions_events":sum(groups["known_taxon_zone"]),
            "heldout_taxon_zone_unseen_elsewhere_events":sum(groups["unseen_taxon_zone"]),
            "by_training_support":views,
            "no_source_region_data_in_fit":True,
        }
    def summary_gain(view:str,metric:str,k:str)->dict:
        eligible=[r for r in outcomes.values() if all(r["by_training_support"][view][p]["identified"] for p in PLOTS)]
        if not eligible:
            return {"identified":False,"n_regions":0}
        v1=[r["by_training_support"][view]["1"]["gain"][metric][k] for r in eligible]
        v3=[r["by_training_support"][view]["3"]["gain"][metric][k] for r in eligible]
        diffs=[a-b for a,b in zip(v1,v3)]
        n=len(eligible)
        return {
            "identified":True,
            "source_region_codes_equal_weight_n":n,
            "mean_region_equal_roadside_gain":mean(v1),
            "mean_region_equal_interior_gain":mean(v3),
            "mean_roadside_minus_interior_gain":mean(diffs),
            "region_roadside_gain_positive":sum(x>0 for x in v1),
            "region_interior_gain_positive":sum(x>0 for x in v3),
            "region_roadside_gain_exceeds_interior":sum(x>0 for x in diffs),
            "leave_one_region_out_difference_min_max":[
                min((sum(diffs)-x)/(n-1) for x in diffs),
                max((sum(diffs)-x)/(n-1) for x in diffs)
            ] if n>1 else None,
            "not_a_random_sampling_country_CI":True
        }
    summary={
        view:{
            metric:{label:summary_gain(view,metric,label) for label in PRIMARY}
            for metric in ("logloss","brier")
        } for view in ("all","known_taxon","unseen_taxon","known_taxon_zone","unseen_taxon_zone")
    }
    return {
        "source_region_code_heldout_count":len(region_names),
        "source_heldout_folds_all_roads_and_sites_in_one_region_excluded":True,
        "source_feature_dimensions":feature_meta,
        "source_region_specific_results":outcomes,
        "equal_region_predictive_gain":summary,
        "all_source_2012plus_and_repeated_survey_protocol_limits_retained":True
    }


def run(source_path:Path)->dict:
    source,srcqc=source_read(source_path)
    pairs,protocol=evaluate_pairs(source)
    rows,cohort=risk_set(source,pairs)
    qa={
        "sites":cohort["source_original_field_sites"],
        "roads":cohort["source_true_independent_roads"],
        "regions":cohort["source_region_codes"],
        "events":cohort["source_early_positive_taxon_plot_risk_events"],
        "taxa":cohort["source_taxa_early_positive"]
    }
    if qa!=EXPECTED:
        raise ValueError(f"Frozen v26.6 authentic 2012+ ROAD source cohort drift: {qa}")
    result=assess(rows,cohort)
    return {
        "version":"v26.7",
        "source_original_MIREN_QA":srcqc,
        "original_protocol_transition_QA":protocol,
        "frozen_source_cohort":cohort,
        "result":result,
        "source_authentic_original_data_DOI":"10.5281/zenodo.15749536",
        "actual_study":"global MIREN mountain ROADS, never new railway evidence",
        "original_source_treatment_NOT_random_or_infrastructure_specific":True,
        "no_regional_future_target_outcome_in_any_training_fold":True,
        "source_Brier_logloss_are_recording_predictions_not_fitness":True,
        "no_train_frequency_or_seed_causality":True,
        "no_raw_taxon_site_location_records_output":True
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--original-csv",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    output=run(args.original_csv)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    for view in ("all","known_taxon","unseen_taxon","known_taxon_zone"):
        print(json.dumps({
            "version":"v26.7","view":view,
            "source_sites":output["frozen_source_cohort"]["source_original_field_sites"],
            "region_loho_n":output["result"]["source_region_code_heldout_count"],
            "mean_heldout_gain":output["result"]["equal_region_predictive_gain"][view]
        },ensure_ascii=False))


if __name__=="__main__":
    main()
