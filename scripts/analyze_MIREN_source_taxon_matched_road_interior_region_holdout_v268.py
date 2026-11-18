#!/usr/bin/env python3
"""v26.8 genuine MIREN ROAD: match exactly the same taxon at roadside and
interior within each original physical site, then compare its out-of-REGION
re-record prediction improvement from initial abundance rank.

Why: v26.7 shows that ~70% of heldout taxon×plot observations have species
not trained in any other region. A stronger roadside–interior comparator is
the SAME originally observed taxon in BOTH original source Plot1 and Plot3
at the SAME spatial Site/Year, rather than separate floras.

Train all models on fixed 14,682 baseline source risk records, same 7 REGION
holdouts as v26.7 (no region indicator, no heldout 2017 data). Test ONLY
same-site taxon pairs observed in BOTH plot types in first survey, and
compare the heldout logloss gain in Plot1 vs Plot3 for exactly those taxa.

Scores weighted per matched taxon within original site, then site within
road, road within Region, equal genuine Region codes. Sensitivity:
require >=2 or >=4 matched-taxonomic sites on >=2 roads per region.
Source plant pair is not a randomized road intervention or railway study;
selection of taxa present in both plot types conditions on habitat affinity,
site-specific source pool and different detection opportunities.

Source data only as authenticated by v26.5–26.7, MIREN original 2025 v3
MD5 b640195cc44f252ccb6ddb28eb6fff13; no raw plant-by-site exported.
"""
from __future__ import annotations
import argparse, json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

import numpy as np
from scipy.sparse import csr_matrix,hstack

from audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 import source_read,evaluate_pairs
from analyze_MIREN_road_abundance_after_taxon_out_of_road_v266 import risk_set, centered_rank_for_fold, fit_prediction
from analyze_MIREN_source_region_holdout_abundance_v267 import (
    features_without_region,score_one_model,LABELS,EXPECTED,C,PLOTS
)

SCENARIOS=("M2_to_M3","M2Z_to_M3Z","M2Z_to_M4Z")
THRESHOLDS=(2,4)


def matched_taxa_at_same_original_source_site(rows:list[dict])->tuple[dict,dict]:
    idx=defaultdict(dict)
    for i,r in enumerate(rows):
        key=(r["region"],r["road"],r["site"],r["taxon"])
        if r["plot"] in idx[key]:
            raise ValueError("Original site/year/plot taxon was counted twice")
        idx[key][r["plot"]]=i
    same=defaultdict(list)
    for key,byplot in idx.items():
        if set(byplot)==set(PLOTS):
            sitekey=key[:3]
            i,j=(byplot[p] for p in PLOTS)
            if (rows[i]["initial_year"]!=rows[j]["initial_year"] or
                rows[i]["year_gap"]!=rows[j]["year_gap"]):
                raise ValueError("Same original taxon compared across nonmatching survey years")
            same[sitekey].append((i,j))
    if not same:
        raise ValueError("No original early source taxon in both paired habitats")
    no_pairs={k[:3] for k in idx}-set(same)
    qa={
        "original_sites_with_at_least_one_2012plus_matched_taxon":len(same),
        "source_same_site_taxon_recorded_in_both_habitats_pairs":sum(map(len,same.values())),
        "original_2012plus_spatial_sites_with_zero_shared_taxa":len(no_pairs),
        "original_unique_early_source_site_count":len({k[:3] for k in idx}),
        "no_extra_source_species_or_sites_created":True,
    }
    return dict(same),qa


def evaluate(rows:list[dict],pairs:dict)->dict:
    X,y,byreg,features,keys=features_without_region(rows)
    source_grouped=defaultdict(list)
    for site in pairs:
        source_grouped[site[0]].append(site)
    all_risks=len(rows)
    rank=np.asarray([r["xrank"] for r in rows],dtype=float)
    results={}
    for region in sorted(byreg):
        test=np.asarray(byreg[region],dtype=int)
        mask=np.ones(all_risks,dtype=bool);mask[test]=False
        train=np.flatnonzero(mask)
        if any(rows[i]["region"]==region for i in train):
            raise AssertionError("Heldout region re-recording target leaked")
        center_train,center_test=centered_rank_for_fold(rank,keys,train,test)
        p={}
        for label in LABELS[:4]:
            p[label]=fit_prediction(X[label][train],X[label][test],y[train],C)
        p[LABELS[4]]=fit_prediction(
            hstack([X[LABELS[2]][train],csr_matrix(center_train.reshape(-1,1))],format="csr"),
            hstack([X[LABELS[2]][test],csr_matrix(center_test.reshape(-1,1))],format="csr"),
            y[train],C
        )
        ll={};brier={}
        for label in LABELS:
            ll[label],brier[label]=score_one_model(y[test],p[label])
        test_index={int(g):i for i,g in enumerate(test)}
        source_sites={}
        for site in source_grouped[region]:
            paired=pairs[site]
            if any(i not in test_index or j not in test_index for i,j in paired):
                raise AssertionError("A shared taxon pair not held out together")
            gain={}
            for metric,loss in (("logloss",ll),("brier",brier)):
                diffs={name:[] for name in SCENARIOS}
                for i,j in paired:
                    a,b=test_index[i],test_index[j]
                    byhab={
                        "M2_to_M3":(
                            loss[LABELS[0]][a]-loss[LABELS[1]][a],
                            loss[LABELS[0]][b]-loss[LABELS[1]][b]),
                        "M2Z_to_M3Z":(
                            loss[LABELS[2]][a]-loss[LABELS[3]][a],
                            loss[LABELS[2]][b]-loss[LABELS[3]][b]),
                        "M2Z_to_M4Z":(
                            loss[LABELS[2]][a]-loss[LABELS[4]][a],
                            loss[LABELS[2]][b]-loss[LABELS[4]][b])
                    }
                    for name,(ra,rb) in byhab.items():
                        diffs[name].append((float(ra),float(rb)))
                gain[metric]={
                    name:{"roadside":mean(a for a,b in vals),
                          "interior":mean(b for a,b in vals),
                          "roadside_minus_interior":mean(a-b for a,b in vals)}
                    for name,vals in diffs.items()
                }
            outcomes=Counter(
                ("both" if y[i] and y[j] else
                 "only_roadside" if y[i] else
                 "only_interior" if y[j] else "neither")
                for i,j in paired
            )
            source_sites[site]={
                "n_source_species_matched_in_both_2012plus_habitats":len(paired),
                "source_2017_rerecord_status_pair_counts":dict(outcomes),
                "mean_original_roadside_minus_interior_early_rank":
                    mean(rank[i]-rank[j] for i,j in paired),
                "gain":gain
            }
        byroad=defaultdict(list)
        for (reg,road,site),d in source_sites.items():
            byroad[road].append(d)
        results[region]={
            "original_region_code":region,
            "source_roads_with_at_least_one_matched_site":len(byroad),
            "original_spatial_sites_with_shared_taxa":len(source_sites),
            "original_2012plus_source_same_taxon_both_zones_pairs":sum(
                d["n_source_species_matched_in_both_2012plus_habitats"] for d in source_sites.values()
            ),
            "original_roads_with_matched_sites":{road:len(sites) for road,sites in sorted(byroad.items())},
            "region_source_site_results_internal":source_sites,
            "no_individual_source_species_taxon_exported":True
        }
    return results


def summarize(results:dict,threshold:int)->dict:
    if threshold not in THRESHOLDS:raise ValueError("Unexpected threshold")
    region_values={}
    for region,rep in results.items():
        byroad=defaultdict(list)
        for (reg,road,site),data in rep["region_source_site_results_internal"].items():
            byroad[road].append(data)
        roads={road:records for road,records in byroad.items() if len(records)>=threshold}
        if len(roads)<2:continue
        metric_out={}
        for metric in ("logloss","brier"):
            metric_out[metric]={}
            for name in SCENARIOS:
                metric_out[metric][name]={
                    component:mean(mean(
                        site["gain"][metric][name][component] for site in sites)
                        for road,sites in roads.items())
                    for component in ("roadside","interior","roadside_minus_interior")
                }
        region_values[region]={
            "original_roads_in_support":len(roads),
            "original_sites_in_support":sum(len(s) for s in roads.values()),
            "original_matched_species_site_pairs_in_support":sum(
                x["n_source_species_matched_in_both_2012plus_habitats"]
                for sites in roads.values() for x in sites
            ),
            "original_macro_equal_road_scores":metric_out
        }
    summary={}
    for metric in ("logloss","brier"):
        summary[metric]={}
        for name in SCENARIOS:
            vals=[d["original_macro_equal_road_scores"][metric][name] for d in region_values.values()]
            if not vals:
                summary[metric][name]={"identified":False,"source_region_codes":0}
                continue
            diffs=[x["roadside_minus_interior"] for x in vals]
            n=len(diffs)
            summary[metric][name]={
                "identified":True,
                "source_region_codes":n,
                "mean_region_equal_roadside_gain":mean(x["roadside"] for x in vals),
                "mean_region_equal_interior_gain":mean(x["interior"] for x in vals),
                "mean_region_equal_roadside_minus_interior_gain":mean(diffs),
                "regions_roadside_minus_interior_gain_positive":sum(x>0 for x in diffs),
                "regions_roadside_minus_interior_gain_negative":sum(x<0 for x in diffs),
                "delete_one_region_range":[min((sum(diffs)-x)/(n-1) for x in diffs),
                                          max((sum(diffs)-x)/(n-1) for x in diffs)] if n>1 else None,
                "not_field_randomization_CI":True
            }
    return {
        "source_min_original_matched_taxon_sites_per_road":threshold,
        "source_min_roads_per_region":2,
        "source_region_codes_identified":len(region_values),
        "source_roads_in_support":sum(v["original_roads_in_support"] for v in region_values.values()),
        "source_spatial_sites_in_support":sum(v["original_sites_in_support"] for v in region_values.values()),
        "source_shared_taxon_site_pairs_in_support":sum(
            v["original_matched_species_site_pairs_in_support"] for v in region_values.values()),
        "region_code_results":region_values,
        "equal_region_results":summary
    }


def run(source_file:Path)->dict:
    source,qc=source_read(source_file)
    original,protocol=evaluate_pairs(source)
    rows,cohort=risk_set(source,original)
    qa={
        "sites":cohort["source_original_field_sites"],"roads":cohort["source_true_independent_roads"],
        "regions":cohort["source_region_codes"],"events":cohort["source_early_positive_taxon_plot_risk_events"],
        "taxa":cohort["source_taxa_early_positive"]
    }
    if qa!=EXPECTED:raise ValueError(f"True MIREN 2012+ 7Region cohort drift: {qa}")
    paired,pairqa=matched_taxa_at_same_original_source_site(rows)
    results=evaluate(rows,paired)
    summary={str(threshold):summarize(results,threshold) for threshold in THRESHOLDS}
    if summary["2"]["source_region_codes_identified"]<3:
        raise ValueError("Insufficient source macro regions with matched source taxa in both original plots")
    # Do not emit original site IDs or taxa; only original region/road counts
    # and source macro summaries may leave this process.
    regional_support={
        region:{
            "original_roads_with_matched_sites":len(v["original_roads_with_matched_sites"]),
            "original_spatial_sites_with_shared_taxa":v["original_spatial_sites_with_shared_taxa"],
            "original_species_site_matched_pairs":v["original_2012plus_source_same_taxon_both_zones_pairs"]
        } for region,v in results.items()
    }
    return {
        "version":"v26.8",
        "original_MIREN_source_QA":qc,
        "original_source_protocol_QA":protocol,
        "frozen_2012plus_original_ROAD_cohort":cohort,
        "source_same_early_site_taxon_in_two_habitats_QA":pairqa,
        "region_support_before_macro_threshold":regional_support,
        "equal_road_equal_region_macro_sensitivity":summary,
        "sampling_condition":"Originally recorded same taxon present in both roadside and interior first-year plots in same Site",
        "whole_original_REGION_out_no_source_re_recording_target_leakage":True,
        "same_taxon_both_habitat_comparison_not_randomized_plot_allocation":True,
        "original_Site_ROAD_region_and_taxon_events_NOT_independent_replicates":True,
        "no_source_taxon_original_GPS_or_plot_rows_exported":True,
        "no_causal_rail_specific_mowing_wind_succession_or_survival_identified":True,
        "source_original_MIREN_ROAD_is_not_a_second_railway_study":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out=run(a.original_csv)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({
        "version":"v26.8",
        "support":out["source_same_early_site_taxon_in_two_habitats_QA"],
        "threshold_summary":{k:{
            "region_codes":v["source_region_codes_identified"],
            "roads":v["source_roads_in_support"],
            "sites":v["source_spatial_sites_in_support"],
            "matched_taxon_site_pairs":v["source_shared_taxon_site_pairs_in_support"],
            "macro_logloss_gain":v["equal_region_results"]["logloss"]
        } for k,v in out["equal_road_equal_region_macro_sensitivity"].items()}
    },ensure_ascii=False))


if __name__=="__main__":main()
