#!/usr/bin/env python3
"""v345: same-author and same S2 block strict-vs-broad label selection audit.

Do not overinterpret different source plot cohorts as a randomized change.
Full CzechVeg+EIVE original hashes inherited from v344 source_context.
Original paper IDs, raw vegetation rows and GPS never exported.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from statistics import mean

from analyze_CzechVeg_strict_HABITAT_EIVE_v344 import source_context
from analyze_CzechVeg_original_study_EIVE_thermal_v343 import (
    frozen_frame,report_stats,trait_data_to_positions
)

PANELS=((3,.5),(3,.7))
ZONES=("rail","road")


def restrict_papers(selected,papers):
    return {k:v for k,v in selected.items() if k[0] in papers}


def restrict_keys(selected,keys):
    return {k:v for k,v in selected.items() if k in keys}


def paper_comparison(strict,strict_scores,broad,broad_scores,traits):
    tr=trait_data_to_positions(traits)
    a,paper_a=report_stats(strict,strict_scores,tr)
    b,paper_b=report_stats(broad,broad_scores,tr)
    paper_a=paper_a or {}
    paper_b=paper_b or {}
    both=sorted(set(paper_a)&set(paper_b))
    diffs=[paper_b[p]-paper_a[p] for p in both]
    return {"strict":{"n_source_papers":a["n_original_studies"],
                       "n_source_blocks":a.get("n_original_matched_blocks",0),
                       "mean_EIVE":a.get("mean_study_equal_RAIL_minus_ROAD_EIVE"),
                       "n_positive":a.get("n_original_studies_positive_RAIL_minus_ROAD")},
            "broad":{"n_source_papers":b["n_original_studies"],
                     "n_source_blocks":b.get("n_original_matched_blocks",0),
                     "mean_EIVE":b.get("mean_study_equal_RAIL_minus_ROAD_EIVE"),
                     "n_positive":b.get("n_original_studies_positive_RAIL_minus_ROAD")},
            "source_papers_present_in_both_estimates":len(both),
            "same_paper_equal_mean_broad_minus_strict":
                mean(diffs) if diffs else None,
            "same_paper_broad_minus_strict_positive_count":sum(x>0 for x in diffs),
            "same_paper_broad_minus_strict_negative_count":sum(x<0 for x in diffs),
            "same_paper_difference_range":[min(diffs),max(diffs)] if diffs else None,
            "different_original_source_plot_subsets_not_causal_label_effect":True}


def same_original_plots(strict,broad,strict_scores,broad_scores,traits):
    common=set(strict)&set(broad)
    matched={}
    for key in common:
        arm={}
        for zone in ZONES:
            shared=set(strict[key][zone]) & set(broad[key][zone])
            if len(shared)<2:
                break
            arm[zone]=tuple(sorted(shared))
        if len(arm)==2:
            matched[key]=arm
    a,paper_a=report_stats(matched,strict_scores,trait_data_to_positions(traits))
    b,paper_b=report_stats(matched,broad_scores,trait_data_to_positions(traits))
    if a.get("n_original_matched_blocks",0)!=b.get("n_original_matched_blocks",0):
        raise AssertionError("Identical original source plot frames differ between labels")
    if paper_a!=paper_b:
        raise AssertionError("Exact same source plant plots must have exactly same temperature CWM")
    return {"source_S2_blocks_both_labels_have_identical_RAIL_ROAD_plot_IDs":len(matched),
            "n_original_papers_with_same_source_plots":a["n_original_studies"],
            "n_original_source_plots_RAIL":sum(len(v["rail"]) for v in matched.values()),
            "n_original_source_plots_ROAD":sum(len(v["road"]) for v in matched.values()),
            "same_actual_plant_plot_temperature_estimates_exactly_identical":True,
            "identical_plot_mean_EIVE":a.get("mean_study_equal_RAIL_minus_ROAD_EIVE"),
            "unavailable_if_original_plot_overlap_below_two_per_arm":True}


def analyze(zenodo,headers,species,eive):
    broad,strict,plots,traits,qa=source_context(zenodo,headers,species,eive)
    out={"version":"v345-original-CzechVeg-same-paper-source-label-audit",
         "original_MD5_source_QA":{
             "headers":qa["genuine_original_czechveg_headers"],
             "species":qa["genuine_original_czechveg_species"],
             "temperature":qa["genuine_original_EIVE"]
         },"panels":[],
         "one_country_observational_source_study_collection":True,
         "not_same_site_rail_natural_control":True}
    for n_systems,coverage in PANELS:
        sf,ss,sqa=frozen_frame(strict,plots,traits,n_systems,coverage)
        bf,bs,bqa=frozen_frame(broad,plots,traits,n_systems,coverage)
        strict_papers={k[0] for k in sf}
        strict_broad_papers=restrict_papers(bf,strict_papers)
        keys=set(sf)&set(bf)
        strict_common=restrict_keys(sf,keys)
        broad_common=restrict_keys(bf,keys)
        pre={"EIVE_min_source_systems":n_systems,
             "source_min_trait_coverage":coverage,
             "unrestricted_original_strict_papers":len(strict_papers),
             "unrestricted_original_strict_blocks":len(sf),
             "unrestricted_original_broad_papers":len({k[0] for k in bf}),
             "unrestricted_original_broad_blocks":len(bf),
             "source_strict_paper_broad_matched_blocks":len(strict_broad_papers),
             "source_strict_and_broad_S2_block_keys_in_common":len(keys),
             "matched_original_papers_source_both_label_modes":
                 paper_comparison(sf,ss,strict_broad_papers,bs,traits),
             "exact_same_original_S2_block_keys":
                 paper_comparison(strict_common,ss,broad_common,bs,traits),
             "identical_original_source_plant_plot_support":
                 same_original_plots(strict_common,broad_common,ss,bs,traits),
             "low_power_if_shared_source_publications_under5":True,
             "historical_HABITAT_AND_LOCALITY_source_labels_not_active_trackbed":True,
        }
        out["panels"].append(pre)
    out["limitations"]=[
       "v345 designed AFTER observing v343 and v344, exploratory post hoc source-cohort audit",
       "Same source publication is NOT the same vegetation plot, locality or randomized treatment",
       "Same source study region decade area block still may use different rail or road plot IDs",
       "Identical plot test only validates statistical formula, cannot adjudicate causal temperature filtering",
       "No original publication title, plot ID, source species x plot records or GPS in results",
       "European EIVE is an ordinal indicator not a microclimate measurement",
       "No train seed transport, genuine independent matched nonrail nature, or global replication"
    ]
    return out


def main():
    parser=argparse.ArgumentParser()
    for f in ("zenodo-record","headers","species","eive","output"):
        parser.add_argument("--"+f,required=True)
    a=parser.parse_args()
    out=analyze(a.zenodo_record,a.headers,a.species,a.eive)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
