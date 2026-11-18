#!/usr/bin/env python3
"""v26.1: author-original Kashmir initial numeric abundance rank versus 2017 re-recording,
conditional on source taxon identity, with complete spatial-site holdout.

SOURCE: Rashid et al 2021, Dryad official version 99735 Matrix_plot.xlsx,
SHA256 529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474.
The v25.4 authentic ZIP/XML reader must pass source bytes/counts before this
analysis. Same 27 source locations (15 station / 12 between), each with
2014/2017 in all four 100m2 zones, 123 source taxon columns. Entire spatial
site (all zones, both dates) withheld as a fold: NO y leakage between zones.

Prediction population: source 2014 PRESENCE (numeric >0) in each plot. Target
is same source species re-RECORDED (numeric >0) in corresponding 2017 plot.
True biological mortality, individual survival, recruitment, propagule source,
mowing, wind and train transport are NOT observed. 2014 abundance used only
for WITHIN original plot midrank among recorded 2014 taxa, [-1,+1]. Numeric
units are not asserted as canopy cover percentages.

Five pre-fixed penalized logistic models (identical C and site folds):
  M0 = eight original group x zone categorical intercepts only.
  M1 = M0 + normalized 2014 within-plot abundance rank.
  M2 = M0 + original taxon-identity indicators (L2 shrunk).
  M3 = M2 + 2014 rank.
  M4 = M3 + rank x track zone interaction (sensitivity, larger model).
Each complete original spatial site is held out in turn. All models have the
same already-known 2014 risk set and 2017 labels. Fixed taxon universe from
2014+2017 publisher column metadata (not labels). Diagnostic results
summarized by original station/between and habitat; no site/taxon rows exported.
Site-equal heldout log-loss and Brier (positive gain = better prediction).
Fold gains are NOT independent field or country replication.

No hyperparameter tuning on heldout outcomes, no search for favorable C.
Posterior/source coefficients not causal fitness effects. Pre-frozen C=1.0;
C=0.2/5 optional exploratory regularization sensitivity (not primary).
"""
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

import numpy as np
from sklearn.linear_model import LogisticRegression

from analyze_Kashmir_same_sites_reference_year_v254 import (
    parse_xlsx_matrix, matching_sites
)
from audit_Kashmir_2014_species_cover_persistence_2017_v259 import early_numeric_cover

GROUPS=("station", "between")
ZONES=("Track", "Embankment", "Interior Plot 1", "Interior Plot 2")
MODELS=("M0_group_zone", "M1_group_zone_rank", "M2_group_zone_taxon",
        "M3_group_zone_taxon_rank","M4_group_zone_taxon_rank_track_interaction")
GZ_SIZE=len(GROUPS)*len(ZONES)
N_TAXA=123
SEED=261118
BASE_RIDGE_C=1.0


def midranks(values:dict[int,float]) -> dict[int,float]:
    """Midrank scaled [-1,1], defined for ties and singleton plots."""
    if not values: raise ValueError("No genuinely recorded 2014 taxa")
    taxa=sorted(values)
    if len(taxa)==1: return {taxa[0]:0.0}
    ordered=sorted((values[t],t) for t in taxa)
    out={}
    start=0
    while start<len(ordered):
        end=start+1
        while end<len(ordered) and ordered[end][0]==ordered[start][0]:
            end+=1
        mid=(start+end-1)/2
        scaled=2*mid/(len(ordered)-1)-1
        for _,taxon in ordered[start:end]:
            out[taxon]=scaled
        start=end
    if set(out)!=set(values) or any(not -1-1e-9<=v<=1+1e-9 for v in out.values()):
        raise AssertionError("Midrank bookkeeping failed")
    return out


def load_authentic_source(path:Path) -> tuple[list[dict],dict]:
    presence,source=parse_xlsx_matrix(path)
    original_cover,audit=early_numeric_cover(path,presence)
    sites={g:matching_sites(presence,g,ZONES) for g in GROUPS}
    if len(sites["station"])!=15 or len(sites["between"])!=12:
        raise ValueError("Original 27 complete spatial-site cohort changed")
    records=[]
    for group in GROUPS:
        for site in sites[group]:
            for zone in ZONES:
                early=original_cover.get((site,2014,zone))
                future=presence.get((site,2017,zone))
                if early is None or future is None:raise ValueError("Original field data incomplete for strict cohort")
                ranks=midranks(early)
                for taxon in sorted(early):
                    records.append({
                        "site":site, "group":group, "zone":zone,
                        "taxon":taxon,
                        "x2014_rank":ranks[taxon],
                        "x2014_raw":early[taxon],
                        "re_recorded_2017":int(taxon in future)
                    })
    count_by_zone={}
    for group in GROUPS:
        count_by_zone[group]={}
        for zone in ZONES:
            d=[x for x in records if x["group"]==group and x["zone"]==zone]
            if len(d)<20:raise ValueError("Insufficient original taxon-site risk events")
            count_by_zone[group][zone]={
                "original_sites":len(sites[group]),
                "source_2014_species_site_risk_events":len(d),
                "source_2017_re_recorded_2014_species_site_events":sum(x["re_recorded_2017"] for x in d),
                "source_event_fraction_re_recorded":sum(x["re_recorded_2017"] for x in d)/len(d)
            }
    return records,{
        "original_Dryad_publisher_SHA256":source["original_file_sha256"],
        "original_source_structural_QA":source,
        "source_2014_numeric_risk_matrix_reconciliation":audit,
        "original_source_sites_per_group":{g:len(sites[g]) for g in GROUPS},
        "original_records_in_source_2014_risk_sets":len(records),
        "source_support_by_group_and_zone":count_by_zone,
        "2014_known_species_universe_original_column_count":N_TAXA,
        "original_site_fold_unit_all_zones_and_years":True
    }


def make_design(rows:list[dict]) -> tuple[np.ndarray,np.ndarray,dict]:
    n=len(rows)
    base=np.zeros((n,GZ_SIZE),dtype=float)
    tax=np.zeros((n,N_TAXA),dtype=float)
    rank=np.zeros((n,1),dtype=float)
    inter=np.zeros((n,1),dtype=float)
    y=np.zeros(n,dtype=int)
    bysite=defaultdict(list)
    for i,r in enumerate(rows):
        site=r["site"]
        bysite[site].append(i)
        grp=GROUPS.index(r["group"])
        zone=ZONES.index(r["zone"])
        base[i,grp*len(ZONES)+zone]=1.
        taxon=r["taxon"]
        if not isinstance(taxon,int) or not 0<=taxon<N_TAXA:
            raise ValueError("Original published species column outside fixed source universe")
        tax[i,taxon]=1.
        rank[i,0]=r["x2014_rank"]
        inter[i,0]=rank[i,0]*(zone==0)
        y[i]=r["re_recorded_2017"]
        if y[i] not in (0,1):
            raise ValueError("Original source re-recording must be binary")
    matrices={
        MODELS[0]:base,
        MODELS[1]:np.column_stack((base,rank)),
        MODELS[2]:np.column_stack((base,tax)),
        MODELS[3]:np.column_stack((base,tax,rank)),
        MODELS[4]:np.column_stack((base,tax,rank,inter)),
    }
    if len(bysite)!=27 or any(len(bysite[s])<4 for s in bysite):
        raise ValueError("Source original spatial site count wrong")
    return y,matrices,dict(bysite)


def fit_predict(xtrain:np.ndarray,ytrain:np.ndarray,xhold:np.ndarray,C:float) -> tuple[np.ndarray,np.ndarray]:
    if len(set(ytrain.tolist()))!=2:
        raise ValueError("Original train fold has no persistence outcome variation")
    clf=LogisticRegression(
        C=C,penalty="l2",solver="liblinear",fit_intercept=True,
        max_iter=1500,tol=1e-8,random_state=SEED
    )
    clf.fit(xtrain,ytrain)
    if not clf.n_iter_[0]<clf.max_iter:
        raise RuntimeError("Original logistic solver failed convergence gate")
    p=np.clip(clf.predict_proba(xhold)[:,1],1e-10,1-1e-10)
    return p,clf.coef_.ravel()


def row_losses(y:np.ndarray,p:np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    if len(y)!=len(p) or np.any((p<=0)|(p>=1)):
        raise ValueError("Expected genuine probability predictions")
    return -(y*np.log(p)+(1-y)*np.log1p(-p)),(y-p)**2


def summarize_site_mean(site_scores:list[dict],model:str)->dict:
    if not site_scores:return {"identified":False,"n_spatial_sites":0}
    vec=[x["model_scores"][model] for x in site_scores]
    return {
        "mean_heldout_site_logloss_nats":mean(x["logloss"] for x in vec),
        "mean_heldout_site_Brier":mean(x["brier"] for x in vec)
    }


def fold_evaluation(rows:list[dict],ridge_C:float=BASE_RIDGE_C) -> dict:
    y,design,by_site=make_design(rows)
    N=len(rows)
    predictions={m:np.zeros(N,dtype=float) for m in MODELS}
    train_species_heldout={}
    for site in sorted(by_site):
        test_ix=np.array(by_site[site],dtype=int)
        keep=np.ones(N,dtype=bool)
        keep[test_ix]=False
        train_ix=np.flatnonzero(keep)
        if len(train_ix)+len(test_ix)!=N or any(rows[i]["site"]==site for i in train_ix):
            raise AssertionError("Original entire-spatial-site holdout train contamination")
        seen={rows[i]["taxon"] for i in train_ix}
        train_species_heldout[site]={i:(rows[i]["taxon"] in seen) for i in test_ix}
        for model in MODELS:
            p,_=fit_predict(design[model][train_ix],y[train_ix],
                           design[model][test_ix],ridge_C)
            predictions[model][test_ix]=p
    if any(np.any(predictions[m]<=0) or np.any(predictions[m]>=1) for m in MODELS):
        raise ValueError("Some original source observations not predicted out-of-site")

    events={}
    for m in MODELS:
        logloss,brier=row_losses(y,predictions[m])
        events[m]=(logloss,brier)
    per_site_zone=[]
    for site,indices in by_site.items():
        for zone in ZONES:
            ids=[i for i in indices if rows[i]["zone"]==zone]
            if not ids:raise ValueError("Original complete site missing zone risk set")
            r=rows[ids[0]]
            if any(rows[i]["group"]!=r["group"] for i in ids):
                raise AssertionError("Original site group incoherent")
            per_site_zone.append({
                "group":r["group"],"zone":zone,
                "n_2014_at_risk":len(ids),
                "n_re_recorded_2017":int(y[ids].sum()),
                "n_taxa_seen_on_other_training_sites":sum(train_species_heldout[site][i] for i in ids),
                "model_scores":{
                    m:{"logloss":float(np.mean(events[m][0][ids])),
                       "brier":float(np.mean(events[m][1][ids]))} for m in MODELS
                }
            })
    if len(per_site_zone)!=27*4:raise AssertionError("Must evaluate 27 heldout original sites x 4 zones")

    contrasts={}
    for group in GROUPS:
        contrasts[group]={}
        for zone in ZONES:
            site_scores=[r for r in per_site_zone if r["group"]==group and r["zone"]==zone]
            avg={m:summarize_site_mean(site_scores,m) for m in MODELS}
            def gain(a,b,metric="logloss"):
                return avg[a]["mean_heldout_site_"+metric+("_nats" if metric=="logloss" else "")] - \
                    avg[b]["mean_heldout_site_"+metric+("_nats" if metric=="logloss" else "")]
            contrasts[group][zone]={
                "n_original_source_heldout_site_folds":len(site_scores),
                "source_2014_at_risk_taxon_site_events":sum(x["n_2014_at_risk"] for x in site_scores),
                "source_2017_re_recorded_original_taxon_site_events":sum(x["n_re_recorded_2017"] for x in site_scores),
                "source_2014_risk_taxa_previously_seen_on_another_training_site_events":
                    sum(x["n_taxa_seen_on_other_training_sites"] for x in site_scores),
                "models":avg,
                "M0_minus_M1_rank_over_zone_logloss_gain":gain(MODELS[0],MODELS[1]),
                "M0_minus_M2_identity_over_zone_logloss_gain":gain(MODELS[0],MODELS[2]),
                "M2_minus_M3_rank_added_after_taxon_logloss_gain":gain(MODELS[2],MODELS[3]),
                "M2_minus_M4_rank_and_track_interaction_added_after_taxon_logloss_gain":gain(MODELS[2],MODELS[4]),
                "M2_minus_M3_rank_added_after_taxon_Brier_gain":gain(MODELS[2],MODELS[3],"brier"),
                "M2_minus_M4_rank_and_track_interaction_added_after_taxon_Brier_gain":gain(MODELS[2],MODELS[4],"brier"),
                "n_source_site_folds_M3_improves_M2":sum(
                    x["model_scores"][MODELS[3]]["logloss"]<x["model_scores"][MODELS[2]]["logloss"]
                    for x in site_scores),
                "n_source_site_folds_M4_improves_M2":sum(
                    x["model_scores"][MODELS[4]]["logloss"]<x["model_scores"][MODELS[2]]["logloss"]
                    for x in site_scores),
                "not_independent_fold_counts_or_p_values":True
            }
    # Taxon unseen before: M2 and M3 cannot learn species coefficients for the
    # heldout species, but pooled rank predictive value can transfer.
    source_seen={}
    for seen in (True,False):
        idx=np.array([i for site,items in train_species_heldout.items()
                      for i,is_seen in items.items()
                      if is_seen==seen and rows[i]["zone"]=="Track"],dtype=int)
        if len(idx)==0:
            source_seen["seen" if seen else "unseen"]={"identified":False}
            continue
        score={}
        for m in (MODELS[2],MODELS[3],MODELS[4]):
            score[m]={"logloss_mean_per_taxon_site_event":float(np.mean(events[m][0][idx])),
                      "brier_mean_per_taxon_site_event":float(np.mean(events[m][1][idx]))}
        source_seen["seen" if seen else "unseen"]={
            "identified":True,
            "source_event_count_NOT_independent_sites":len(idx),
            "scores":score,
            "event_weighted_not_primary_site_weighted":True
        }
    return {
        "C_regularization_fixed_without_hyperparameter_search":ridge_C,
        "leave_entire_original_spatial_site_out_27_fold":True,
        "models":list(MODELS),
        "source_predictions_by_group_and_habitat":contrasts,
        "Track_species_seen_vs_unseen_in_other_sites":source_seen,
        "primary_metric":"Original site-equal out-of-site logloss (nat), separately by group x habitat; gain = identity-only M2 minus +rank M3",
        "secondary_metric":"Original site-equal heldout Brier; rank x Track interaction M4 is sensitivity",
        "site_folds_are_overlapping_train_sets_not_independent_field_experiments":True,
        "initial_2014_numeric_abundance_as_within_plot_rank_not_absolute_percent_cover":True,
        "no_train_seed_delivery_individual_survival_or_genetic_adaptation":True
    }


def run(path:Path,C:float=BASE_RIDGE_C)->dict:
    if C<=0:raise ValueError("Ridge C must be positive")
    rows,audit=load_authentic_source(path)
    primary=fold_evaluation(rows,C)
    return {
        "version":"v26.1",
        "source_authentic_original_field_data":audit,
        "result":primary,
        "analysis_scope":"Single genuine Kashmir longitudinal railway study, not global meta-analysis",
        "evidence_level":"conditional out-of-spatial-site predictive information, not biological fitness or causal selection",
        "predefined_site_holdout_no_outcome_tuning":True,
        "label_not_colonization_nor_survival":"presence in 2017 original source after being recorded in 2014",
        "no_original_taxon_site_abundance_or_GPS_export":True
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--C",default=BASE_RIDGE_C,type=float)
    args=parser.parse_args()
    out=run(args.source,args.C)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    for g,zones in out["result"]["source_predictions_by_group_and_habitat"].items():
        for z,r in zones.items():
            print(json.dumps({
                "version":"v26.1","source_group":g,"source_habitat":z,
                "n_source_spatial_sites":r["n_original_source_heldout_site_folds"],
                "n_risk_events":r["source_2014_at_risk_taxon_site_events"],
                "n_rerecorded":r["source_2017_re_recorded_original_taxon_site_events"],
                "M0_logloss":r["models"][MODELS[0]]["mean_heldout_site_logloss_nats"],
                "M2_taxon_only_logloss":r["models"][MODELS[2]]["mean_heldout_site_logloss_nats"],
                "M3_taxon_rank_logloss":r["models"][MODELS[3]]["mean_heldout_site_logloss_nats"],
                "M4_taxon_rank_track_logloss":r["models"][MODELS[4]]["mean_heldout_site_logloss_nats"],
                "M2_minus_M3":r["M2_minus_M3_rank_added_after_taxon_logloss_gain"],
                "M2_minus_M4":r["M2_minus_M4_rank_and_track_interaction_added_after_taxon_logloss_gain"],
                "positive_M3_folds":r["n_source_site_folds_M3_improves_M2"]
            },ensure_ascii=False))
if __name__=="__main__":
    main()
