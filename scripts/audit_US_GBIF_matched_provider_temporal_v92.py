#!/usr/bin/env python3
"""v9.2 public GBIF datasetKey same-provider temporal audit.

Eight FROZEN U.S. Census states (v8.1 original pilot), 2018–2021
versus 2023–2025 occurrence event years, 2022 excluded. The *same*
datasetKey must be represented in BOTH periods. Identify top 3 common
dataset providers by min(earlier, later) ALL-Plantae source records,
NEVER by plant species or rail hypothesis sign. Independently source
facet each common provider by taxon speciesKey and separate its GRIIS
US-country-listed taxon categories from unknown/nonlisted records.

A late-only GBIF speciesKey is NOT an actual new plant immigrant;
matching dataset UUID does not match individual observer/field plot,
source protocol, or identify upload/backfill. No train seed carriage,
railroad mowing or wind adaptability measured.
"""
from __future__ import annotations
import argparse,json,math,time,uuid
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
from statistics import median
import requests
from harvest_US_state_GRIIS_GBIF_plant_source_v81 import (
    CENSUS_URL,PILOT_FIPS,parse_census_state_table,source_get)
from harvest_US_state_GRIIS_temporal_split_v85 import PERIODS,query_state_years
from harvest_gbif_cell_species_facets_v20 import parse_page,FACET_PAGE
from audit_gbif_detected_original_rows_v48b import parse_status_facet
from analyze_griis_alien_rail_source_v25 import valid_griis

N_PROVIDER=3
MIN_ALLPLANT_BOTH_PERIOD_PROVIDER_RECORDS=50
MIN_US_GRIIS_BOTH_PROVIDER_RECORDS=20
MAX_SAME_YEAR_DATASET_SOURCE_COUNT_DRIFT=.025
DATASET_FACET_TOP=100
MAX_TAXON_PAGES=12
SAMPLE_IS_FROZEN_ORIGINAL_EIGHT_STATES=True

def common_provider_params(state_name,window):
    out=query_state_years(state_name,window)
    out["facet"]=["datasetKey","occurrenceStatus"]
    out["facetLimit"]=DATASET_FACET_TOP
    out["datasetKey.facetLimit"]=DATASET_FACET_TOP
    out.pop("speciesKey.facetOffset",None)
    out.pop("speciesKey.facetLimit",None)
    return out

def provider_taxon_params(state_name,window,dataset_key,offset=0):
    if offset<0 or offset%FACET_PAGE or not isinstance(dataset_key,str):
        raise ValueError("Invalid provider speciesKey source parameters")
    try:uuid.UUID(dataset_key)
    except ValueError:raise ValueError("GBIF datasetKey UUID required")
    out=query_state_years(state_name,window,offset)
    out["datasetKey"]=dataset_key
    out["facet"]=["speciesKey","occurrenceStatus"]
    out["facetLimit"]=FACET_PAGE
    out["speciesKey.facetLimit"]=FACET_PAGE
    return out

def provider_facet(session,state_name,window,baseline_count):
    raw=source_get(session,common_provider_params(state_name,window))
    n,providers=parse_page(raw,"datasetKey")
    status=parse_status_facet(raw)
    if status["indexed_source_count"]!=n:
        raise ValueError("Dataset/status source count mismatch")
    if baseline_count<=0 or abs(n-baseline_count)/baseline_count>MAX_SAME_YEAR_DATASET_SOURCE_COUNT_DRIFT:
        raise ValueError("Live original GBIF source changed >2.5% since independently verified v8.5 indexed source")
    counts={}
    for key,value in providers:
        try:
            normalized=str(uuid.UUID(key))
        except ValueError:
            raise ValueError("GBIF datasetKey facet must be a UUID")
        if normalized in counts:raise ValueError("Duplicate top source dataset UUID")
        counts[normalized]=value
    if sum(counts.values())>n:
        raise ValueError("Top100 dataset counts exceed occurrence source count")
    return {"source_record_count":n,"top_dataset_UUID_counts":counts,
            "top100_source_count_coverage":round(sum(counts.values())/n,6),
            "n_top_facet_provider_categories":len(counts),
            "top100_not_the_exhaustive_provider_universe":True}

def select_prior_common_provider_pairs(early,late):
    before=early["top_dataset_UUID_counts"]
    after=late["top_dataset_UUID_counts"]
    common=set(before)&set(after)
    choices=[p for p in common
             if min(before[p],after[p])>=MIN_ALLPLANT_BOTH_PERIOD_PROVIDER_RECORDS]
    choices.sort(key=lambda k:(-min(before[k],after[k]),
                               -(before[k]+after[k]),k))
    return choices[:N_PROVIDER],{"n_top100_dataset_UUIDs_in_both_windows":len(common),
        "n_shared_providers_meeting_50_AllPlantae_records_each_window":len(choices),
        "n_selected_preoutcome_source_provider_UUIDs":min(N_PROVIDER,len(choices)),
        "provider_selected_independent_of_GRIIS_taxa_or_FAF_rail_effect":True}

def provider_GRIIS_facets(session,name,window,provider,expected_total,accepted):
    source={}
    count=None
    n_pages=0
    statuses=None
    for page in range(MAX_TAXON_PAGES):
        response=source_get(session,provider_taxon_params(name,window,provider,page*FACET_PAGE))
        actual,cats=parse_page(response,"speciesKey")
        if count is None:
            count=actual
            status=parse_status_facet(response)
            statuses=status
            if status["indexed_source_count"]!=actual:
                raise ValueError("Provider-filtered species status and count diverged")
        elif actual!=count:
            raise ValueError("Provider-filtered source count drift over taxon pagination")
        for key,num in cats:
            if key in source:raise ValueError("Duplicate speciesKey across provider source facet pages")
            if num>0:source[key]=num
        n_pages+=1
        if len(cats)<FACET_PAGE:break
    else:
        raise ValueError("Provider-filtered GBIF species facet pagination exhausted")
    if expected_total<=0 or abs(count-expected_total)/expected_total>MAX_SAME_YEAR_DATASET_SOURCE_COUNT_DRIFT:
        raise ValueError("Provider source filtered number changed versus independently obtained datasetKey facet")
    if sum(source.values())>count:
        raise ValueError("Provider taxa source count exceeds independent occurrence total")
    unknown=(statuses["unclassified_count"]+statuses["statuses"].get("NOT_DETECTED",0))
    if unknown>count:
        raise ValueError("Provider non-DETECTED/unknown count invalid")
    griis={k:num for k,num in source.items() if k in accepted}
    return {
      "original_indexed_all_Plantae_records":count,
      "n_US_country_GRIIS_indexed_records":sum(griis.values()),
      "n_US_country_GRIIS_recorded_species_key_categories":len(griis),
      "US_country_GRIIS_positive_speciesKey_record_counts":griis,
      "source_non_DETECTED_or_unknown_fraction":round(unknown/count,7) if count else None,
      "complete_faceted_taxon_source_pages":n_pages,
      "same_dataset_UUID_does_not_identify_same_observers_or_plots":True,
    }

def evaluate_state(before,after,orig_state,summaries,chosen,accepted):
    """Conserve original all-provider source species set; never compare
    raw novelty counts without annotating how much source was excluded."""
    e=orig_state["source_windows"]["before_2022"]
    l=orig_state["source_windows"]["after_2022"]
    early=set(e["US_country_GRIIS_positive_speciesKey_record_counts"])
    late=set(l["US_country_GRIIS_positive_speciesKey_record_counts"])
    original_new=late-early
    provider_early=set()
    provider_late=set()
    provider_e_records=0
    provider_l_records=0
    source_support=0
    usable=0
    for provider in chosen:
        ev=summaries[provider]["before_2022"]
        lv=summaries[provider]["after_2022"]
        if (ev["n_US_country_GRIIS_indexed_records"]>=MIN_US_GRIIS_BOTH_PROVIDER_RECORDS and
            lv["n_US_country_GRIIS_indexed_records"]>=MIN_US_GRIIS_BOTH_PROVIDER_RECORDS):
            usable+=1
            provider_e_records+=ev["original_indexed_all_Plantae_records"]
            provider_l_records+=lv["original_indexed_all_Plantae_records"]
            provider_early.update(ev["US_country_GRIIS_positive_speciesKey_record_counts"])
            provider_late.update(lv["US_country_GRIIS_positive_speciesKey_record_counts"])
        source_support+=1
    if not provider_early.issubset(early) or not provider_late.issubset(late):
        raise ValueError("Independently selected provider-specific taxa are not subset of original GBIF all-source state-period taxa")
    provider_only=provider_late-provider_early
    supported_new=original_new&provider_late
    return {
        "original_all_provider_country_GRIIS_early_source_taxa":len(early),
        "original_all_provider_country_GRIIS_late_source_taxa":len(late),
        "original_late_only_source_taxon_categories_NOT_colonizations":len(original_new),
        "n_preoutcome_same_UUID_providers_selected":source_support,
        "n_same_UUID_providers_with_both_20record_GRIIS_support":usable,
        "two_period_common_provider_source_support":usable>=1,
        "n_matched_provider_early_GRIIS_source_taxon_categories":len(provider_early),
        "n_matched_provider_late_GRIIS_source_taxon_categories":len(provider_late),
        "matched_provider_only_late_taxon_categories":len(provider_only),
        "n_matched_provider_late_only_already_recorded_early_by_OTHER_datasets":len(provider_only&early),
        "n_original_allsource_late_only_taxa_also_visible_from_retained_providers":len(supported_new),
        "fraction_original_late_only_supported_by_common_UUID_providers":
           round(len(supported_new)/len(original_new),6) if original_new else None,
        "retained_common_providers_early_allPlantae_record_coverage":
           round(provider_e_records/e["n_all_Plantae_indexed_records"],6) if usable else 0.,
        "retained_common_providers_late_allPlantae_record_coverage":
           round(provider_l_records/l["n_all_Plantae_indexed_records"],6) if usable else 0.,
        "recorded_late_only_taxa_NOT_new_biological_introductions":True,
    }

def run(session,census,listed,baseline):
    if baseline.get("complete") is not True or baseline.get("source_periods")!=PERIODS:
        raise ValueError("Original v8.5 independent before/post 2022 GBIF source not available")
    if set(PILOT_FIPS)-set(baseline["states"]):
        raise ValueError("Original source 8 Census state identities missing")
    if len(census)!=51 or len(listed)<100:
        raise ValueError("Original Census FIPS or country-GRIIS source missing")
    out={
        "complete":True,"source_utc":datetime.now(timezone.utc).isoformat(),
        "study":"REAL U.S. v8.1 fixed 8-state shared GBIF datasetKey source provider pre/post audit",
        "frozen_state_codes":list(PILOT_FIPS),"source_periods":PERIODS,
        "provider_source_rule":"top 3 UUIDs in both 2018-21 and 2023-25 all-Plantae top100 facets with at least 50 records in each, rank by minimum records; GRIIS >=20 each required for usable source",
        "no_2022_plant_record_year_in_either_source":True,
        "states":{},"source_errors":{},
        "NOT_actual_rail_seed_dispersal_or_new_invaders":True}
    for code in PILOT_FIPS:
        name=census[code]["name"]
        if baseline["states"][code]["state_name"]!=name or not baseline["states"][code][
                "two_year_source"]["two_period_quality"]:
            raise ValueError("Original v8.1 geography and independently prequalified two-window source changed")
        state={"state_name":name,"fips":code,"source_provider_facet":{},"selected_shared_provider_UUIDs":[],
               "source_provider_subset":{}}
        try:
            for w in PERIODS:
                n=baseline["states"][code]["source_windows"][w]["n_all_Plantae_indexed_records"]
                state["source_provider_facet"][w]=provider_facet(session,name,w,n)
            picked,stats=select_prior_common_provider_pairs(
                state["source_provider_facet"]["before_2022"],
                state["source_provider_facet"]["after_2022"])
            state["preoutcome_source_provider_selection"]=stats
            state["selected_shared_provider_UUIDs"]=picked
            for provider in picked:
                state["source_provider_subset"][provider]={}
                for win in PERIODS:
                    exp=state["source_provider_facet"][win]["top_dataset_UUID_counts"][provider]
                    state["source_provider_subset"][provider][win]=provider_GRIIS_facets(
                        session,name,win,provider,exp,listed)
            state["source_comparison"]=evaluate_state(
                state["source_provider_facet"]["before_2022"],
                state["source_provider_facet"]["after_2022"],baseline["states"][code],
                state["source_provider_subset"],picked,listed)
            state["complete"]=True
        except Exception as e:
            state["complete"]=False
            state["source_error"]=type(e).__name__+": "+str(e)[:300]
            out["source_errors"][code]=state["source_error"]
        out["states"][code]=state
        summary=state.get("source_comparison",{})
        print(json.dumps({
            "state":code,"complete":state["complete"],
            "selected_shared_UUID_sources":len(state["selected_shared_provider_UUIDs"]),
            "usable_GRIIS_sources":summary.get("n_same_UUID_providers_with_both_20record_GRIIS_support"),
            "matched_provider_record_coverage_late":summary.get("retained_common_providers_late_allPlantae_record_coverage"),
            "baseline_post_only_source_taxa":summary.get("original_late_only_source_taxon_categories_NOT_colonizations"),
            "retained_provider_post_only_source_taxa":summary.get("n_original_allsource_late_only_taxa_also_visible_from_retained_providers"),
            "n_errors":len(out["source_errors"])},ensure_ascii=False),flush=True)
    out["complete"]=len(out["states"])==len(PILOT_FIPS) and not out["source_errors"]
    comparisons=[r["source_comparison"] for r in out["states"].values() if r.get("complete")]
    usable=[r for r in comparisons if r["two_period_common_provider_source_support"]]
    frac=[r["fraction_original_late_only_supported_by_common_UUID_providers"]
          for r in usable if r["fraction_original_late_only_supported_by_common_UUID_providers"] is not None]
    out["summary"]={
        "n_source_completed_8_original_states":len(comparisons),
        "n_states_with_at_least_one_same_UUID_providing_20_GRIIS_indexed_records_in_each_window":len(usable),
        "median_fraction_of_baseline_late_only_source_taxa_visible_in_matched_provider_late_sources":
            round(median(frac),6) if frac else None,
        "median_matched_provider_fraction_of_late_allPlantae_indexed_records":
            round(median(r["retained_common_providers_late_allPlantae_record_coverage"] for r in usable),6) if usable else None,
        "n_states_with_any_false_provider_only_late_taxa_previously_recorded_in_other_dataset":
            sum(r["n_matched_provider_late_only_already_recorded_early_by_OTHER_datasets"]>0
                for r in usable),
    }
    out["limitations"]=[
        "GBIF datasetKey UUID identity only fixes SOURCE DATASET, not observer, revisits, survey sites, detection protocol, sampling effort or 2026 GBIF retrospective data uploads.",
        "Top100 dataset facets selected on all-Plantae volume may omit many smaller earlier/later providers; matched-provider fraction cannot be used as unweighted randomized new-species colonization probability.",
        "Matched provider-source record support may fall to zero in some states. Missing taxa are unrecorded, NOT biological absence or locally introduced species.",
        "Country-level U.S. GRIIS alien status is NOT local state-native status. The same species may be naturally resident in a state for decades.",
        "2018-21 and 2023-25 differ in length. Source-record depth/effort, geographic location and taxonomy may vary over time, even under same dataset UUID.",
        "2022 FAF6 economic OD or physical 2026 mapped track path is not a record of viable seeds, train freight car contamination, trackside flora, mowing/gust exposure or evolved adaptation."
    ]
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("griis","baseline85","output"):p.add_argument("--"+k,required=True)
    a=p.parse_args()
    load=lambda x:json.loads(Path(getattr(a,x)).read_text(encoding="utf-8"))
    census_request=requests.get(CENSUS_URL,timeout=(20,90))
    census_request.raise_for_status()
    census=parse_census_state_table(census_request.text)
    accepted=valid_griis(load("griis")).get("US")
    if not accepted:raise ValueError("U.S. GRIIS source lacks accepted taxa")
    with requests.Session() as session:
        session.headers.update({"User-Agent":"tetsu-GBIF-eightstate-matched-provider-source-v92"})
        out=run(session,census,accepted,load("baseline85"))
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],
        "source_errors":out["source_errors"],"NOT_real_plant_colonization":True},
        ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Matched dataset provider source incomplete; do not infer true new species.")

if __name__=="__main__":
    main()
