#!/usr/bin/env python3
"""v26.6 — ORIGINAL MIREN mountain ROAD (NOT RAILWAY), source-verified
leave-whole-original-road-out abundance rank incremental prediction.

Primary cohort = v26.5 genuine 2012+ first/last repeat sites, strict source
Plot1 roadside + Plot3 interior, 2+ roads per Region code with 4+ valid
paired sites per road. Same original plant records, author source MD5.
The original 2014-style early source risk set is all detected taxa in
the first visit; target whether the SAME accepted taxon is re-recorded
at that source plot at the last eligible visit.

Fit five L2 logistic models (C=1.0 fixed; no heldout hyper-tuning):
 M2   source Region+Plot+initial YEAR+interval+taxon identity
 M3   M2 + original 1..3 abundance WITHIN-PLOT midrank [-1,+1]
 M2Z  M2 with taxon×Plot rather than taxon alone
 M3Z  M2Z + same within-plot abundance rank
 M4Z  M2Z + source rank centered by that taxon×plot mean ONLY among
       other source ROADS in training (no heldout road outcomes).
One fold = ALL source sites and both plots of one whole original road
held out; never train on same road/site/2017 outcome. Score original
sites equally within road, roads equally within Region, Regions equally.

This does not convert 25 MIREN mountain-ROAD regions into railroad
replications, does not measure actual mortality, fitness, train effects,
or rail mowing/wind selection. Earliest/latest Protocol versions differ
in all genuine repeated sites, another unavoidable observation caveat.
No original taxon, site, GPS, or cover values written to output.
"""
from __future__ import annotations
import argparse, json, math
from collections import defaultdict
from pathlib import Path
from statistics import mean
import numpy as np
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
from scipy.sparse import hstack, csr_matrix

from audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 import (
    source_read,evaluate_pairs,aggregate,PRIMARY_MIN_YEAR,
    MIN_SITES_PER_ROAD,MIN_ROADS_PER_REGION
)

LABELS=("M2_taxon","M3_taxon_plus_rank",
        "M2Z_taxon_by_zone","M3Z_taxon_by_zone_plus_rank",
        "M4Z_taxon_by_zone_training_centered_rank")
C_FIXED=1.0
SEED=266118
PLOTS=("1","3")


def midrank(values:dict[str,float])->dict[str,float]:
    if len(values)<2:raise ValueError("2012+ source risk set needs multiple taxa")
    v=sorted(values.values())
    n=len(v)
    out={}
    for taxon,value in values.items():
        lo=sum(x<value for x in v)
        eq=sum(x==value for x in v)
        out[taxon]=2*(lo+(eq-1)/2)/(n-1)-1
    if not all(-1-1e-8<=x<=1+1e-8 for x in out.values()):
        raise ValueError("Source plot ordinal abundance midranks out of range")
    return out


def risk_set(source:dict,pairs:list[dict])->tuple[list[dict],dict]:
    main=aggregate(pairs,"Abundance","post2012")
    regions=main["original_region_code_aggregate_only_not_country_replicates"]
    allowed_region_roads={}
    for region in regions:
        counts=defaultdict(int)
        for item in pairs:
            if item["region"]!=region or item["initial_year"]<PRIMARY_MIN_YEAR:
                continue
            if all(item[(z,"Abundance")]["identified"] for z in PLOTS):
                counts[item["road"]]+=1
        allowed_region_roads[region]={road for road,n in counts.items() if n>=MIN_SITES_PER_ROAD}
        if len(allowed_region_roads[region])<MIN_ROADS_PER_REGION:
            raise ValueError("Inconsistent original road and source regional support")
    plotdata=source["plants"]
    rows=[]
    bysite=defaultdict(set)
    for item in pairs:
        region,road,site,yr0,yr1=(
            item["region"],item["road"],item["site"],
            item["initial_year"],item["final_year"])
        if yr0<PRIMARY_MIN_YEAR:
            continue
        if road not in allowed_region_roads.get(region,set()):
            continue
        if not all(item[(z,"Abundance")]["identified"] for z in PLOTS):
            continue
        for plot in PLOTS:
            early=plotdata[(region,road,site,yr0,plot)]
            later=plotdata[(region,road,site,yr1,plot)]
            raw={taxon: metrics[0] for taxon,metrics in early.items()}
            if any(x is None for x in raw.values()):
                raise ValueError("MIREN 2012+ numeric early abundance must be complete")
            ranks=midrank(raw)
            for taxon in sorted(raw):
                rows.append({
                    "region":region,"road":road,"site":site,"plot":plot,
                    "initial_year":yr0,"year_gap":yr1-yr0,
                    "taxon":taxon,"xrank":ranks[taxon],
                    "re_recorded":int(taxon in later)
                })
            bysite[region,road,site].add(plot)
    if not rows or any(plots!=set(PLOTS) for plots in bysite.values()):
        raise ValueError("Original source road site paired plot cohort violated")
    if (len(bysite)!=main["source_sites_in_supported_equal_region_summary"]
        or len({(r,d) for r,d,s in bysite})!=main["source_roads_in_supported_equal_region_summary"]):
        raise ValueError(f"v26.5 cohort gate: expected site={main['source_sites_in_supported_equal_region_summary']} road={main['source_roads_in_supported_equal_region_summary']}; got site={len(bysite)} road={len({(r,d) for r,d,s in bysite})}; by region={dict(sorted(__import__('collections').Counter(r for r,d,s in bysite).items()))}")
    qa={
        "verified_original_MIREN_MD5":"b640195cc44f252ccb6ddb28eb6fff13",
        "source_original_field_sites":len(bysite),
        "source_true_independent_roads":len({(r,d) for r,d,s in bysite}),
        "source_region_codes":len({r for r,d,s in bysite}),
        "source_early_positive_taxon_plot_risk_events":len(rows),
        "source_taxa_early_positive":len({r["taxon"] for r in rows}),
        "source_original_first_year_min":min(r["initial_year"] for r in rows),
        "source_original_last_year_max":max(r["initial_year"]+r["year_gap"] for r in rows),
        "no_original_railway_sites":True
    }
    return rows,qa


def make_design(rows:list[dict])->tuple[dict,np.ndarray,dict,dict,list[tuple[str,str]]]:
    base=DictVectorizer(sparse=True)
    taxa=DictVectorizer(sparse=True)
    zonetaxa=DictVectorizer(sparse=True)
    B=base.fit_transform([{
        "region:"+r["region"]:1.0,
        "plot:"+r["plot"]:1.0,
        "initial_year:"+str(r["initial_year"]):1.0,
        "elapsed:"+str(r["year_gap"]):1.0
    } for r in rows])
    T=taxa.fit_transform([{"taxon:"+r["taxon"]:1.0} for r in rows])
    TZ=zonetaxa.fit_transform([{"taxon_zone:"+r["taxon"]+"|"+r["plot"]:1.0} for r in rows])
    rank=np.array([r["xrank"] for r in rows],dtype=float)
    R=csr_matrix(rank.reshape(-1,1))
    X={
        LABELS[0]:hstack([B,T],format="csr"),
        LABELS[1]:hstack([B,T,R],format="csr"),
        LABELS[2]:hstack([B,TZ],format="csr"),
        LABELS[3]:hstack([B,TZ,R],format="csr")
    }
    y=np.array([r["re_recorded"] for r in rows],dtype=int)
    roads=defaultdict(list)
    for i,r in enumerate(rows):
        roads[(r["region"],r["road"])].append(i)
    return X,y,dict(roads),{"base":B.shape[1],
         "distinct_source_taxa":T.shape[1],
         "distinct_source_taxon_zone":TZ.shape[1],
         "source_first_year_and_interval_categories_are_training_known_or_zero_for_unseen":True},[
             (r["taxon"],r["plot"]) for r in rows]


def centered_rank_for_fold(rank:np.ndarray,keys:list[tuple[str,str]],
                           train:np.ndarray,test:np.ndarray)->tuple[np.ndarray,np.ndarray]:
    by=defaultdict(list)
    for i in train:by[keys[i]].append(rank[i])
    mu={k:mean(xs) for k,xs in by.items()}
    tr=np.array([rank[i]-mu[keys[i]] for i in train])
    te=np.array([rank[i]-mu[keys[i]] if keys[i] in mu else 0. for i in test])
    if abs(np.mean(tr))>1e-7:
        raise ValueError("Training-only taxon-zone centering did not remove overall rank mean")
    return tr,te


def fit_prediction(Xtrain,Xtest,ytrain:np.ndarray,C:float=C_FIXED)->np.ndarray:
    if set(ytrain.tolist())!={0,1}:
        raise ValueError("Source original training risk set lacks both outcomes")
    clf=LogisticRegression(C=C,solver="liblinear",fit_intercept=True,
                           max_iter=1000,tol=1e-6,random_state=SEED)
    clf.fit(Xtrain,ytrain)
    if clf.n_iter_[0]>=clf.max_iter:
        raise RuntimeError("Fixed ridge model did not converge")
    return np.clip(clf.predict_proba(Xtest)[:,1],1e-10,1-1e-10)


def run(path:Path)->dict:
    source,source_qc=source_read(path)
    pairs,protocol=evaluate_pairs(source)
    rows,cohort=risk_set(source,pairs)
    X,y,roads,design,keys=make_design(rows)
    n=len(rows)
    results=[]
    ordinal=np.array([r["xrank"] for r in rows],dtype=float)
    from scipy.sparse import hstack
    for region,road in sorted(roads):
        test=np.array(roads[region,road],dtype=int)
        mask=np.ones(n,dtype=bool);mask[test]=False
        train=np.flatnonzero(mask)
        if any(rows[i]["region"]==region and rows[i]["road"]==road for i in train):
            raise AssertionError("Source road fold leaked into model training")
        pred={}
        for lab in LABELS[:4]:
            pred[lab]=fit_prediction(X[lab][train],X[lab][test],y[train])
        trc,tec=centered_rank_for_fold(ordinal,keys,train,test)
        pred[LABELS[4]]=fit_prediction(
            hstack([X[LABELS[2]][train],csr_matrix(trc.reshape(-1,1))],format="csr"),
            hstack([X[LABELS[2]][test],csr_matrix(tec.reshape(-1,1))],format="csr"),
            y[train])
        ll={lab:-(y[test]*np.log(p)+(1-y[test])*np.log1p(-p))
            for lab,p in pred.items()}
        br={lab:(y[test]-p)**2 for lab,p in pred.items()}
        grouped=defaultdict(list)
        for i,j in enumerate(test):
            grouped[(rows[j]["site"],rows[j]["plot"])].append(i)
        for plot in PLOTS:
            same=[ids for (site,z),ids in grouped.items() if z==plot]
            if not same:raise ValueError("True road cannot omit original paired plot")
            # First equal source plot species risk set -> original site equal,
            # then original road equal -> original Region equal.
            site_ll={lab:mean(float(np.mean(ll[lab][q])) for q in same) for lab in LABELS}
            site_br={lab:mean(float(np.mean(br[lab][q])) for q in same) for lab in LABELS}
            results.append({
                "region":region,"road":road,"plot":plot,"source_site_count":len(same),
                "source_taxon_plot_risk_events":sum(map(len,same)),
                "source_re_recorded":int(sum(y[test[q]].sum() for q in same)),
                "loss":site_ll,"brier":site_br
            })
    reg=defaultdict(lambda:defaultdict(list))
    for q in results:reg[q["region"]][q["plot"]].append(q)
    supported={}
    for region,grp in sorted(reg.items()):
        if len(grp["1"])<2 or len(grp["3"])<2:
            raise ValueError("Source original macro region must have 2+ genuine roads")
        def combine(metric:str,plot:str)->dict:
            rs=grp[plot]
            m2=mean(x[metric][LABELS[0]] for x in rs)
            m3=mean(x[metric][LABELS[1]] for x in rs)
            m2z=mean(x[metric][LABELS[2]] for x in rs)
            m3z=mean(x[metric][LABELS[3]] for x in rs)
            m4z=mean(x[metric][LABELS[4]] for x in rs)
            return {"source_m2":m2,"source_m3":m3,"source_m2z":m2z,
                    "source_m3z":m3z,"source_m4z":m4z,
                    "gain_M2_minus_M3":m2-m3,
                    "gain_M2Z_minus_M3Z":m2z-m3z,
                    "gain_M2Z_minus_M4Z_training_centered":m2z-m4z}
        supported[region]={
            "source_roads":len(grp["1"]),
            "original_spatial_sites":sum(x["source_site_count"] for x in grp["1"]),
            "original_taxon_plot_risk_events":sum(
                x["source_taxon_plot_risk_events"] for x in grp["1"]+grp["3"]),
            "roadside_logloss":combine("loss","1"),
            "interior_logloss":combine("loss","3"),
            "roadside_Brier":combine("brier","1"),
            "interior_Brier":combine("brier","3"),
        }
    fields=("gain_M2_minus_M3","gain_M2Z_minus_M3Z",
            "gain_M2Z_minus_M4Z_training_centered")
    summary={}
    for metric in ("logloss","Brier"):
        for k in fields:
            edge=[v["roadside_"+metric][k] for v in supported.values()]
            inside=[v["interior_"+metric][k] for v in supported.values()]
            delta=[a-b for a,b in zip(edge,inside)]
            summary[metric+"_"+k]={
                "mean_region_equal_roadside_gain":mean(edge),
                "mean_region_equal_interior_gain":mean(inside),
                "mean_region_equal_roadside_minus_interior_gain":mean(delta),
                "regions_rank_gain_positive_roadside":sum(x>0 for x in edge),
                "regions_rank_gain_positive_interior":sum(x>0 for x in inside),
                "regions_roadside_gain_exceeds_interior":sum(x>0 for x in delta),
                "regions_in_support":len(delta),
                "leave_one_region_out_diff_range":[
                    min((sum(delta)-v)/(len(delta)-1) for v in delta),
                    max((sum(delta)-v)/(len(delta)-1) for v in delta)] if len(delta)>1 else None
            }
    return {
        "version":"v26.6",
        "status":"AUTHOR_ORIGINAL_MIREN_ROAD_SOURCE_VERIFIED_IF_WORKFLOW_SUCCESS",
        "source_original_MIREN_QA":source_qc,
        "source_original_protocol_QA":protocol,
        "source_cohort":cohort,"source_model_feature_dimensions":design,
        "original_source_roadheldout_macro_prediction":summary,
        "region_code_source_results":supported,
        "primary_metric":"source site-mean out-of-road logloss, then equal roads, equal region codes",
        "original_2012_plus_earliest_source_ONLY":True,
        "source_midrank_categories_are_1_3_individual_count_classes_not_percent_cover":True,
        "original_source_all_road_protocol_versions_change_between_first_last":True,
        "not_railway_and_not_a_genuine_railway_effect_replication":True,
        "no_causal_survival_mowing_slipstream_or_plant_fitness_identified":True,
        "no_source_taxon_or_plot_level_data_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=run(a.original_csv)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
        "version":result["version"],
        "source_cohort":result["source_cohort"],
        "source_prediction":result["original_source_roadheldout_macro_prediction"]
    },ensure_ascii=False))
if __name__=="__main__":main()
