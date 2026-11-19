#!/usr/bin/env python3
"""v29.0: Kashmir railway life-cycle selection reference and taxon leverage audit.

Scientifically separate: within-original-field-site annual-fraction contrasts,
effect of reference zone, true source original-site leave-one-out, and effect of
dropping a SOURCE taxon everywhere (which changes species fraction denominator).
Any unknown original species receives ONE global annual/nonannual value across
all site/year/zone observations. This does not identify train-mediated seed
dispersal, adaptation, or continent-generalized fitness.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_Kashmir_annual_microsites_v285 import original_taxon_headers, life_span_crosswalk
from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_Kew_curated_strict_annual_v289 import curated_merge
from analyze_Kashmir_global_unknown_trait_partial_id_v288 import analyze, common_source_sites, WEIGHTS

FOCUS_CONTRASTS = (
    ("Track", "Interior Plot 1"),
    ("Track", "Interior Plot 2"),
    ("Track", "Embankment"),
    ("Embankment", "Interior Plot 1"),
    ("Embankment", "Interior Plot 2"),
)
FOCUS_TAXA = (
    "Cynodon dactylon (L.) Pers.",
    "Breea arvensis (L.) Less.",
    "Taraxacum officinale Weber",
    "Lactuca serriola L.",
    "Senecio chrysanthemoides DC.",
    "Lotus corniculatus L.",
)
SOURCE_STUDY = "Rashid_2021_India_Kashmir_railway_2014_2017"

def summary_contrasts(full: dict) -> dict:
    table = {}
    for name, block in full["comparisons"].items():
        table[name] = {}
        for group in ("all","station","between"):
            table[name][group] = {}
            for time in ("2014", "2017", "two_year_mean", "2017_minus_2014"):
                r = block[group][time]
                table[name][group][time] = {
                    "original_sites": r["source_sites"],
                    "lower": r["sharp_lower"],
                    "upper": r["sharp_upper"],
                    "sign_identified": r["sign_identified"],
                    "scenario_all_remaining_nonannual": r.get("all_unlabelled_nonannual_scenario"),
                    "scenario_all_remaining_annual": r.get("all_unlabelled_annual_scenario"),
                }
    return table


def leave_one_original_site(plots: dict, classified: dict[int,str], taxon_count: int) -> dict:
    sites=common_source_sites(plots)
    source=analyze(plots,classified,taxon_count)
    raw=source["comparisons"]
    diagnostics={}
    for group in ("all","station","between"):
        group_sites=[s for s in sites if group=="all" or
                     (s.startswith("Between ") if group=="between"
                      else not s.startswith("Between "))]
        values={}
        for holdout in group_sites:
            without={key:val for key,val in plots.items() if key[0]!=holdout}
            run=analyze(without,classified,taxon_count)
            for zone,ref in FOCUS_CONTRASTS:
                name=zone+"__minus__"+ref
                r=run["comparisons"][name][group]["two_year_mean"]
                values.setdefault(name,[]).append({
                    "holdout":holdout,
                    "lower":r["sharp_lower"],
                    "upper":r["sharp_upper"],
                    "sign":r["sign_identified"]})
        diagnostics[group]={}
        for key, sampled in values.items():
            obs=raw[key][group]["two_year_mean"]
            diagnostics[group][key] = {
                "original_eligible_sites":len(group_sites),
                "original_lower":obs["sharp_lower"],
                "original_upper":obs["sharp_upper"],
                "LOSO_total_folds":len(sampled),
                "LOSO_negative_upper_folds":sum(v["upper"]<0 for v in sampled),
                "LOSO_positive_lower_folds":sum(v["lower"]>0 for v in sampled),
                "LOSO_min_upper":min(v["upper"] for v in sampled),
                "LOSO_max_upper":max(v["upper"] for v in sampled),
                "LOSO_min_lower":min(v["lower"] for v in sampled),
                "LOSO_max_lower":max(v["lower"] for v in sampled),
                "withheld_site_identifiers_not_exported":True,
                "LOSO_is_within_one_independent_original_study":True,
            }
    return diagnostics


def exclude_taxon_everywhere(plots:dict,
                             classified:dict[int,str],
                             original_names:list[str],
                             dropped_name:str)->dict:
    if dropped_name not in original_names:
        raise ValueError("Taxon is not an original published Dryad column")
    source_sites=common_source_sites(plots)
    index=original_names.index(dropped_name)
    remaining={key:(None if species is None else species-{index})
               for key,species in plots.items()}
    eligible=[s for s in source_sites
              if all(remaining.get((s,y,z)) for y in (2014,2017)
                     for z in ("Track","Embankment","Interior Plot 1","Interior Plot 2"))]
    removed_empty=[s for s in source_sites if s not in eligible]
    valid_eligible=set(eligible)
    filtered={key:s for key,s in remaining.items()
              if key[0] in valid_eligible}
    if not eligible:
        return {"dropped_original_taxon":dropped_name,"original_sites":len(source_sites),
                "remaining_complete_source_sites":0,"status":"NO_PAIRED_SOURCE_SUPPORT"}
    check=analyze(filtered,classified,len(original_names))
    effects={}
    for zone,ref in FOCUS_CONTRASTS:
        name=zone+"__minus__"+ref
        effects[name]={group: {
            "lower":check["comparisons"][name][group]["two_year_mean"]["sharp_lower"],
            "upper":check["comparisons"][name][group]["two_year_mean"]["sharp_upper"],
            "sign":check["comparisons"][name][group]["two_year_mean"]["sign_identified"],
            "original_sites":check["comparisons"][name][group]["two_year_mean"]["source_sites"],
        } for group in ("all","station","between")}
    return {"dropped_original_taxon":dropped_name,
            "original_sites":len(source_sites),
            "remaining_complete_source_sites":len(eligible),
            "sites_with_zero_original_species_after_drop":len(removed_empty),
            "dropped_taxon_was_life_span_known":index in classified,
            "comparison_denominators_are_reduced_by_one_taxon_where_observed":True,
            "status":"COMPUTED_WITH_NEW_SPECIES_DENOMINATOR",
            "two_year_mean":effects}


def label_uncertainty_ablation(plots:dict,classified:dict[int,str],
                               original_names:list[str])->dict:
    source=analyze(plots,classified,len(original_names))
    base=source["comparisons"]["Track__minus__Interior Plot 1"]["all"]["two_year_mean"]
    out={"curated_baseline_lower":base["sharp_lower"],
         "curated_baseline_upper":base["sharp_upper"]}
    for name in ("Cynodon dactylon (L.) Pers.","Breea arvensis (L.) Less."):
        idx=original_names.index(name)
        if idx not in classified:
            raise ValueError("Predeclared Kew taxon missing from curated label source")
        changed={k:v for k,v in classified.items() if k!=idx}
        r=analyze(plots,changed,len(original_names))
        rr=r["comparisons"]["Track__minus__Interior Plot 1"]["all"]["two_year_mean"]
        out[name]={"refit_if_life_span_unresolved_lower":rr["sharp_lower"],
                   "refit_if_life_span_unresolved_upper":rr["sharp_upper"],
                   "sign_identified":rr["sign_identified"],
                   "drop_a_confirmed_label_not_drop_the_species":True}
    return out


def produce(plots:dict, classified:dict[int,str], original_names:list[str])->dict:
    all_full=analyze(plots,classified,len(original_names))
    study_sites=common_source_sites(plots)
    return {
       "version":"v29.0",
       "original_field_study":SOURCE_STUDY,
       "distinct_field_studies":1,
       "source_site_count":len(study_sites),
       "source_years":[2014,2017],
       "taxa_original_headers":len(original_names),
       "known_life_span_taxa":len(classified),
       "unknown_life_span_taxa":len(original_names)-len(classified),
       "all_zone_paired_references":summary_contrasts(all_full),
       "whole_original_site_leave_one_out":leave_one_original_site(plots,classified,len(original_names)),
       "original_species_omission_sensitivity": {
           name:exclude_taxon_everywhere(plots,classified,original_names,name)
           for name in FOCUS_TAXA
       },
       "label_certainty_critical_taxa":label_uncertainty_ablation(plots,classified,original_names),
       "scientific_limits":{
          "structurally_different_reference_zones":True,
          "one_independent_source_study_not_globally_replicated":True,
          "site_selection_and_trait_distribution_not_randomized":True,
          "site_leave_one_out_not_independent_country_holdout":True,
          "extreme_trait_bounds_are_not_population_CI":True,
          "plant_fitness_genetic_adaptation_and_train_seed_traffic_not_observed":True,
       },
    }

def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rail",required=True,type=Path)
    p.add_argument("--global-traits",required=True,type=Path)
    p.add_argument("--curated-traits",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    original_names=original_taxon_headers(a.rail)
    plots,railmetadata=parse_xlsx_matrix(a.rail)
    merged,curated_info=curated_merge(original_names,a.global_traits,a.curated_traits)
    idx={name:i for i,name in enumerate(original_names)}
    classified={idx[r["original_taxon"]]:r["life_span"] for r in merged}
    report=produce(plots,classified,original_names)
    report["author_raw_source_sha256"]=railmetadata["original_file_sha256"]
    report["global_life_span_crosswalk_sha256"]=hashlib.sha256(a.global_traits.read_bytes()).hexdigest()
    report["Kew_curation_sha256"]=curated_info["independent_curation_file_sha256"]
    report["site_taxon_observations_not_redistributed"]=True
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    # print human-reviewable summary without original site identifiers
    table=report["all_zone_paired_references"]
    main={k:{
       "all_two_year":v["all"]["two_year_mean"],
       "station_two_year":v["station"]["two_year_mean"],
       "between_two_year":v["between"]["two_year_mean"],
       "annual_2014":v["all"]["2014"],
       "annual_2017":v["all"]["2017"],
       "LOSO":report["whole_original_site_leave_one_out"]["all"][k]
       } for k,v in table.items()}
    print("V290_SOURCE_DIAGNOSTICS "+json.dumps({
        "source_original_sites":report["source_site_count"],
        "known_taxa":report["known_life_span_taxa"],
        "unknown_taxa":report["unknown_life_span_taxa"],
        "comparisons":main,
        "taxon_omission":report["original_species_omission_sensitivity"],
        "curation_dependency":report["label_certainty_critical_taxa"]
    },ensure_ascii=False))


if __name__=="__main__":
    main()
