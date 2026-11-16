#!/usr/bin/env python3
"""v9.4 MACRO observational source sensitivity: REAL original same 2022
FAF6 crop RAIL/TRUCK directional state OD versus 2018–21 origin / 2023–
25 destination GRIIS-US-country-listed GBIF 'first-recorded' taxa.

Decompose GBIF indexed speciesKey record-count facets into exact
iNaturalist Research-grade Observations provider 50c9509d... and all
OTHER providers by PER-SPECIES-KEY subtraction (all - iNat). Qualify
the same 40 U.S. Census states with ≥100 other-source all-plant
records, ≥5 source taxa and ≥20 other-source GRIIS indexed records
in EACH period, prior to inspecting any freight/plant correlation.

Fit identical original geographic freight OD pair cohort across all/
iNaturalist/other provider response definitions. Compare both the
original all-source effort covariates held FIXED and source-specific
effort covariates. No dyadic pseudo-independent significance tests.

Plant 'newly recorded' is NOT local biological immigration, state-
specific alien species status, actual train carriage or adaption to
mowing/wind/gust disturbances. Positive GBIF record counts are not
surveys of absence.
"""
from __future__ import annotations
import argparse,copy,json
from datetime import datetime,timezone
from pathlib import Path
from hashlib import sha256
import requests

from analyze_US_FAF6_directional_source_taxa_v86 import (
    COHORTS,RESPONSES,construct_directions,signed_results,
    partial_rank,cohort)
from analyze_FAF6_GRIIS_US_state_OD_alignment_v83 import (
    CENSUS_TIGER_URL,parse_census_centroids)

GROUPS=("all_GBIF_original","iNaturalist_only","all_OTHER_GBIF_except_iNaturalist")
WINDOWS=("before_2022","after_2022")
MIN_NON_INAT_SOURCE_RECORDS=20
MIN_NON_INAT_SOURCE_ALLPLANT=100
MIN_NON_INAT_SOURCE_TAXA=5
MIN_SAME_PAIR_QUALITY_STATES=40
MIN_SHARED_DIRECTIONAL_PAIRS=750
MAIN_RESPONSE="expected_20record_source_linked_first_record_fraction"


def split_state_source(baseline,inat_source,code,window):
    orig=baseline["states"][code]["source_windows"][window]
    i=inat_source["states"][code]["source_result"][
       "two_GBIF_year_source_periods"][window]
    all_counts=orig["US_country_GRIIS_positive_speciesKey_record_counts"]
    i_counts=i["iNaturalist_source_country_GRIIS_taxon_keys"]
    if not set(i_counts).issubset(all_counts) or any(i_counts[k]>all_counts[k] for k in i_counts):
        raise ValueError("Original source speciesKey record counts don't conserve under provider exclusion")
    other={k:n-i_counts.get(k,0) for k,n in all_counts.items()
           if n-i_counts.get(k,0)>0}
    ia=i["iNaturalist_original_all_Plantae_record_source_count"]
    total=orig["n_all_Plantae_indexed_records"]
    other_all=total-ia
    if total<=0 or other_all<0 or ia<0:
        raise ValueError("Original all-Plantae source record count changed")
    if (sum(other.values())!=orig["n_US_country_GRIIS_indexed_records"]
        -i["iNaturalist_source_country_GRIIS_record_count"]):
        raise ValueError("Country GRIIS source speciesKey indexed records did not conserve")
    return {
       "all_GBIF_original":(dict(all_counts),total),
       "iNaturalist_only":(dict(i_counts),ia),
       "all_OTHER_GBIF_except_iNaturalist":(other,other_all),
    }


def source_quality_states(temporal,inat):
    if (temporal.get("complete") is not True or
        inat.get("complete") is not True or
        inat.get("source_original_GBIF_qualified_states")!=48 or
        len(inat.get("states",{}))!=48):
        raise ValueError("Original 48 state source provider identity not independently validated")
    allstates={}
    eligible=[]
    reason={}
    for code,p in sorted(inat["states"].items()):
        if (not p["complete"] or
            temporal["states"][code]["two_year_source"]["two_period_quality"] is not True):
            raise ValueError("Original prequalified same provider 48-state source mismatch")
        bywin={win:split_state_source(temporal,inat,code,win) for win in WINDOWS}
        allstates[code]=bywin
        failed=[]
        for win in WINDOWS:
            counts,n_all=bywin[win]["all_OTHER_GBIF_except_iNaturalist"]
            if n_all<MIN_NON_INAT_SOURCE_ALLPLANT:failed.append(win+" original source <100 all-Plantae records")
            if len(counts)<MIN_NON_INAT_SOURCE_TAXA:failed.append(win+" <5 other provider country-GRIIS recorded taxa")
            if sum(counts.values())<MIN_NON_INAT_SOURCE_RECORDS:failed.append(win+" <20 other provider GRIIS indexed records")
        if failed:reason[code]=failed
        else:eligible.append(code)
    if len(eligible)<MIN_SAME_PAIR_QUALITY_STATES:
        return allstates,eligible,reason
    return allstates,eligible,reason


def variant_temporal(temporal,allstates,eligible,variant):
    """Only source-qualified common cohort 40 states; never choose
    source/geographic units by source freight partial rho."""
    if variant not in GROUPS:raise ValueError("Unfrozen source provider variant")
    copy_t=copy.deepcopy(temporal)
    for code,state in copy_t["states"].items():
        allowed=code in eligible
        state["two_year_source"]["two_period_quality"]=allowed
        if not allowed:continue
        for window in WINDOWS:
            counts,nall=allstates[code][window][variant]
            orig=state["source_windows"][window]
            orig["n_all_Plantae_indexed_records"]=nall
            orig["n_US_country_GRIIS_indexed_records"]=sum(counts.values())
            orig["n_US_country_GRIIS_observed_speciesKey_categories"]=len(counts)
            orig["US_country_GRIIS_positive_speciesKey_record_counts"]=counts
            orig["source_eligible"]=True
        state["two_year_source"]["n_late_only_first_recorded_source_taxon_keys"]=len(
            set(state["source_windows"]["after_2022"]["US_country_GRIIS_positive_speciesKey_record_counts"])
           -set(state["source_windows"]["before_2022"]["US_country_GRIIS_positive_speciesKey_record_counts"]))
    copy_t["n_qualified_source_states_gate_40"]=len(eligible)>=MIN_SAME_PAIR_QUALITY_STATES
    return copy_t


def same_pair_country_results(rows):
    out={}
    for cohort_name in COHORTS:
        sub=cohort(rows,cohort_name)
        results={}
        for response in RESPONSES:
            results[response]={
                "rail_rho":partial_rank(sub,response,"RAIL") if len(sub)>=30 else None,
                "truck_rho":partial_rank(sub,response,"TRUCK") if len(sub)>=30 else None
            }
        out[cohort_name]={
            "n_source_directed_state_pairs_shared_across_all_providers":len(sub),
            "n_distinct_origin_sources":len({r["o"] for r in sub}),
            "n_distinct_dest_sources":len({r["d"] for r in sub}),
            "source_responses":results,
            "NO_pairwise_independent_replicates":True,
        }
    return out


def run(original85,inat93,fa80,flora82,prior86,census_raw):
    if (prior86.get("complete") is not True or prior86.get(
        "n_evaluable_directed_origin_destination_source_pairs")!=2255 or
        prior86.get("n_both_GBIF_source_quality_US_states")!=48 or
        prior86.get("Census_source_sha256")!=sha256(census_raw).hexdigest()):
        raise ValueError("Original v8.6 48state directional source was modified")
    geo=parse_census_centroids(census_raw.decode("utf-8"),flora82["states"])
    original_rows,original_states,missing=construct_directions(original85,fa80,geo)
    if len(original_rows)!=2255 or len(original_states)!=48:
        raise ValueError("Prior source original 48state candidate pair frame drift")
    for label in COHORTS:
        sample=cohort(original_rows,label)
        for response in RESPONSES:
            rho=partial_rank(sample,response,"RAIL")
            if rho!=prior86["all_results"][label]["results"][response]["rail_full"]:
                raise ValueError("Original original-provider FAF v8.6 rail coefficient not reproduced")
    split,eligible,fail=source_quality_states(original85,inat93)
    out={
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "original_2022_FAF6_estimated_economic_crop_freight_OD_not_real_track_traffic":True,
        "original_US_country_GRIIS_listed_not_confirmed_state_intro_taxa":True,
        "original_GBIF_source_2018_21_vs_2023_25_not_true_biological_absence":True,
        "validated_source_48_GBIF_states":48,
        "number_source_states_with_other_provider_20_GRIIS_5species_100all_records_both_windows":len(eligible),
        "source_eligible_same_provider_analysis_state_codes":eligible,
        "states_without_enough_other_provider_source_opportunity_and_reasons":fail,
        "same_source_state_selection_independent_of_plant_freight_correlation":True,
        "original_v86_48state_2255_directional_pair_coefficients_reproduced":True,
        "n_same_source_directed_pairs":None,
        "signed_common_pair_source":{},
        "allsource_covariates_held_fixed_response_only_sensitivity":{},
        "not_estimated_if_under_40_prequalified_other_provider_states":True,
        "NO_train_seed_movement_or_mowing_gust_evolutionary_causality":True,
    }
    if len(eligible)<MIN_SAME_PAIR_QUALITY_STATES:
        out["source_comparison_estimable"]=False
        out["conclusion"]="Non-iNaturalist GBIF source not enough original state representation; no signed biological effect"
        return out
    variant_maps={}
    for variant in GROUPS:
        synthetic=variant_temporal(original85,split,eligible,variant)
        rows,states,_=construct_directions(synthetic,fa80,geo)
        variant_maps[variant]={(r["o"],r["d"]):r for r in rows}
        if set(states)!=set(eligible):
            raise ValueError("State-level provider source cohort changed under speciesKey filtering")
    shared=sorted(set.intersection(*(set(v) for v in variant_maps.values())))
    out["n_same_source_directed_pairs"]=len(shared)
    if len(shared)<MIN_SHARED_DIRECTIONAL_PAIRS:
        out["source_comparison_estimable"]=False
        out["conclusion"]="Too few strictly comparable original source direction pairs across all/INAT/nonINAT"
        return out
    for variant in GROUPS:
        src=[variant_maps[variant][key] for key in shared]
        out["signed_common_pair_source"][variant]=same_pair_country_results(src)
    original_fixed=[variant_maps[GROUPS[0]][key] for key in shared]
    for variant in GROUPS:
        changed=[]
        for key,base in zip(shared,original_fixed):
            updated=dict(base)
            updated["responses"]=variant_maps[variant][key]["responses"]
            changed.append(updated)
        out["allsource_covariates_held_fixed_response_only_sensitivity"][variant]=same_pair_country_results(changed)
    out["source_comparison_estimable"]=True
    out["quality_methods"]=[
       "NON-iNaturalist source speciesKey record counts are exact per-species subtraction of independent all GBIF facets and UUID-scoped iNaturalist dataset facets with matched source filters.",
       "40 US states selected on non-iNaturalist per-year ORIGINAL source support before examining freight/plant signs; no geographic windows selected on ecological correlation.",
       "Source record candidate pools and K20 rarefaction change when dataset source is removed; the fixed-allsource-controls sensitivity holds geography/freight/record covariates fixed to isolate observational source outcome.",
       "GBIF data are not independent plot surveys, and 2022 state economic FAF OD is not actual train seed movement; outcomes cannot be interpreted as botanical arrivals."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for name in ("temporal85","inat93","fa80","flora82","prior86","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    load=lambda n:json.loads(Path(getattr(a,n)).read_text(encoding="utf-8"))
    census=requests.get(CENSUS_TIGER_URL,timeout=(20,90))
    census.raise_for_status()
    result=run(load("temporal85"),load("inat93"),load("fa80"),
               load("flora82"),load("prior86"),census.content)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    primary={g:result["allsource_covariates_held_fixed_response_only_sensitivity"].get(g,{}).get(
        COHORTS[0],{}).get("source_responses",{}).get(MAIN_RESPONSE,{})
             for g in GROUPS}
    print(json.dumps({
      "complete":result["complete"],
      "n_common_source_states":result["number_source_states_with_other_provider_20_GRIIS_5species_100all_records_both_windows"],
      "n_strict_common_pair_directions":result["n_same_source_directed_pairs"],
      "source_comparison_estimable":result.get("source_comparison_estimable"),
      "allsource_fixed_covariates_allpair_rho":primary,
      "all_source_GRIIS_recording_may_not_indicate_true_species_absence":True
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
