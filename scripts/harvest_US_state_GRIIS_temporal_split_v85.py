#!/usr/bin/env python3
"""v8.5: independently census-fixed 51 U.S. state GBIF GRIIS plant
SOURCE-WINDOW audit, using 2018–21 (pre-2022 FAF base year) and
2023–25 (post) with 2022 NEVER queried/overlapped.

This is NOT species colonization, establishment, true absence,
independent survey detection, state-level alien status, dated
rail freight track/seed movement or railway mowing/wind adaptation.
Raw 2026 indexed speciesKey facets contain historical event YEARs;
dataset upload/backfill and status/provider process are uncontrolled.
"""
from __future__ import annotations
import argparse,json,math,time
from collections import Counter
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import requests

from harvest_US_state_GRIIS_GBIF_plant_source_v81 import (
    CENSUS_URL,GBIF_URL,parse_census_state_table,source_get)
from harvest_gbif_cell_species_facets_v20 import parse_page,FACET_PAGE
from audit_gbif_detected_original_rows_v48b import parse_status_facet
from analyze_griis_alien_rail_source_v25 import valid_griis

PERIODS={"before_2022":"2018,2021","after_2022":"2023,2025"}
SOURCE_QUALITY_MAX_UNKNOWN_OR_NOT_DETECTED=.02
MIN_TOTAL=100
MIN_DATASETS=3
MIN_GRIIS_TAXA=5
MIN_GRIIS_RECORDS=20
MAX_PAGES=12
MIN_BOTH_PERIOD_QUALIFIED_STATE_SOURCES=40

def query_state_years(state_name,window,offset=0):
    if not state_name or window not in PERIODS or type(offset)!=int or offset<0 or offset%FACET_PAGE!=0:
        raise ValueError("FROZEN 2018-21 and 2023-25 GBIF source windows required")
    return {
        "country":"US","stateProvince":state_name,
        "kingdomKey":6,"taxonRank":"SPECIES",
        "year":PERIODS[window],
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "limit":0,"facet":["speciesKey","occurrenceStatus","datasetKey"],
        "facetLimit":FACET_PAGE,
        "speciesKey.facetLimit":FACET_PAGE,"speciesKey.facetOffset":offset,
        "occurrenceStatus.facetLimit":8,"datasetKey.facetLimit":5,
    }

def get_window(session,state_name,window,accepted):
    counts={}
    indexed=None
    source_datasets=None
    non_status=None
    pages=0
    for page in range(MAX_PAGES):
        raw=source_get(session,query_state_years(state_name,window,page*FACET_PAGE))
        count,categories=parse_page(raw,"speciesKey")
        if indexed is None:
            indexed=count
            _,publishers=parse_page(raw,"datasetKey")
            source_datasets=len(publishers)
            status=parse_status_facet(raw)
            if status["indexed_source_count"]!=count:raise ValueError("GBIF status/data-source count discrepancy")
            non_status=status["statuses"].get("NOT_DETECTED",0)+status["unclassified_count"]
        elif indexed!=count:
            raise ValueError("GBIF 2026 search source changed during speciesKey pagination")
        for key,n in categories:
            if key in counts:raise ValueError("Duplicate GBIF speciesKey across windows/paged categories")
            if type(n)!=int or n<0:raise ValueError("Bad original GBIF speciesKey count")
            if n>0:counts[key]=n
        pages+=1
        if len(categories)<FACET_PAGE:break
        time.sleep(.05)
    else:
        raise ValueError("GBIF state/source faceted species list exceeds fixed 60000 speciesKey budget")
    if indexed is None or sum(counts.values())>indexed:
        raise ValueError("GBIF state-period indexed source species records exceed total")
    non_rate=non_status/indexed if indexed>0 else None
    taxa={k:n for k,n in counts.items() if k in accepted}
    reasons=[]
    if indexed<MIN_TOTAL:reasons.append("indexed_all_Plantae_records_below_100")
    if source_datasets<MIN_DATASETS:reasons.append("independently_published_dataset_UUID_categories_below_3")
    if non_rate is None or non_rate>SOURCE_QUALITY_MAX_UNKNOWN_OR_NOT_DETECTED:
        reasons.append("source_NOT_DETECTED_or_unknown_status_above_2pct")
    if len(taxa)<MIN_GRIIS_TAXA:reasons.append("GRIIS_US_country_list_observed_species_keys_below_5")
    if sum(taxa.values())<MIN_GRIIS_RECORDS:reasons.append("GRIIS_US_country_list_indexed_records_below_20")
    return {
        "complete":True,
        "source_eligible":not reasons,
        "n_all_Plantae_indexed_records":indexed,
        "n_all_Plantae_positive_species_key_categories":len(counts),
        "n_dataset_UUID_categories_at_least":source_datasets,
        "source_NOT_DETECTED_or_unknown_record_fraction":
            round(non_rate,7) if non_rate is not None else None,
        "n_US_country_GRIIS_observed_speciesKey_categories":len(taxa),
        "n_US_country_GRIIS_indexed_records":sum(taxa.values()),
        "US_country_GRIIS_positive_speciesKey_record_counts":taxa,
        "source_speciesKey_facet_pages_complete":pages,
        "source_eligibility_failures":reasons,
        "NOT_species_absence_nor_establishment_nor_state_introduced_status":True,
    }

def audit_temporal_window(source_early,source_late):
    """No new-colonization inferred from an indexed taxon recorded in only
    one window; treat it ONLY as first record in this two-bin dataset.
    """
    if not all(x.get("complete") is True for x in (source_early,source_late)):
        return {"two_period_quality":False,"n_early_only_source_taxon_keys":None,
                "n_late_only_first_recorded_source_taxon_keys":None}
    e=source_early["US_country_GRIIS_positive_speciesKey_record_counts"]
    l=source_late["US_country_GRIIS_positive_speciesKey_record_counts"]
    return {
      "two_period_quality":source_early["source_eligible"] and source_late["source_eligible"],
      "n_early_only_source_taxon_keys":len(set(e)-set(l)),
      "n_late_only_first_recorded_source_taxon_keys":len(set(l)-set(e)),
      "n_taxon_keys_recorded_both_year_bins":len(set(e)&set(l)),
      "n_late_only_species_recorded_is_NOT_biological_introduction":True,
    }

def run(griis,states,sess,fa79,fa82):
    if (fa79.get("complete") is not True or fa79.get("FAF_release")!="FAF6.0 / 2022 BASE YEAR / STATE OD CSV"
        or fa82.get("complete") is not True or fa82.get("n_requested_state_fips_codes")!=51
        or set(states)!=set(fa82["states"]) or len(states)!=51):
        raise ValueError("Independent 2022 FAF / 2026 GBIF 51-state Census FIPS source universe changed")
    accepted=valid_griis(griis).get("US")
    if accepted is None or len(accepted)<100:raise ValueError("Original verified U.S. GRIIS backbone incomplete")
    out={
      "complete":True,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "geographical_study":"Census FIPS 51 US states/DC predeclared, with 2022 FAF agricultural freight source",
      "source_periods":PERIODS,
      "source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows":True,
      "counted_entities":"2026 GBIF indexed recorded source speciesKey categories, NOT independently surveyed plant colonization",
      "US_country_GRIIS_accepted_taxon_key_count":len(accepted),
      "states":{},"errors":{},
      "not_seed_delivery_causation_nor_mowing_wind_adaptation":True,
    }
    for code,state in sorted(states.items()):
        periods={}
        for win in PERIODS:
            try:
                periods[win]=get_window(sess,state["name"],win,accepted)
            except Exception as exc:
                periods[win]={"source_error":True}
                out["errors"][f"{code}_{win}"]=type(exc).__name__+": "+str(exc)[:160]
            time.sleep(.12)
        count=audit_temporal_window(periods["before_2022"],periods["after_2022"])
        out["states"][code]={
          "state_name":state["name"],"state_abbreviation":state["abbreviation"],
          "source_windows":periods,"two_year_source":count}
        print(json.dumps({
          "FIPS":code,"early_ok":periods["before_2022"].get("source_eligible"),
          "late_ok":periods["after_2022"].get("source_eligible"),
          "post_only_recorded_taxa":count["n_late_only_first_recorded_source_taxon_keys"],
          "source_errors":len(out["errors"])},ensure_ascii=False),flush=True)
    all_complete=len(out["states"])==51 and not out["errors"]
    qualified={code for code,v in out["states"].items()
               if v["two_year_source"]["two_period_quality"] is True}
    out["complete"]=all_complete
    out["summary"]={
       "n_US_Census_FIPS_states_source_requested":51,
       "n_early_window_GRIIS_and_GBIF_source_eligible":sum(
         x["source_windows"]["before_2022"].get("source_eligible") is True
         for x in out["states"].values()),
       "n_late_window_GRIIS_and_GBIF_source_eligible":sum(
         x["source_windows"]["after_2022"].get("source_eligible") is True
         for x in out["states"].values()),
       "n_both_window_GRIIS_US_recording_source_quality_states":len(qualified),
       "n_source_error_windows":len(out["errors"]),
       "n_country_US_GRIIS_taxon_keys_first_recorded_in_any_qualified_late_state":
           len({taxon for c in qualified
                for taxon in (set(out["states"][c]["source_windows"]["after_2022"][
                    "US_country_GRIIS_positive_speciesKey_record_counts"])
                    -set(out["states"][c]["source_windows"]["before_2022"][
                    "US_country_GRIIS_positive_speciesKey_record_counts"]))}),
    }
    out["n_qualified_source_states_gate_40"]=len(qualified)>=MIN_BOTH_PERIOD_QUALIFIED_STATE_SOURCES
    out["limits"]=[
      "GBIF indexed record YEAR is not publication/accession date; retrospective additions, new observers and source providers can generate fake new taxa without colonization.",
      "Country-wide USA GRIIS listed introduced species may be indigenous within individual states; do not call state-specific alien or new invasion.",
      "Missing source taxa in 2018–21 does not establish biological absence; first positive GBIF key in 2023–25 is NOT colonization or seed arrival.",
      "One 2022 FAF6 estimated economic OD matrix, not observed rails or historical pre/post train operating routes, causes no quasi-experiment.",
      "Different length of 2018–21 (four years) versus 2023–25 (three years) and GBIF dataset/effort drift bias first-record fractions; no inferences until matched record coverage or provider controls.",
      "2022 is excluded from BOTH flora bins, which only establishes ordering of recorded YEAR windows around a freight base year, not true shipment-caused colonization.",
      "No measured train-delivered propagules, rail-versus-road germination, mowing frequency, gust exposure, plasticity or evolutionary adaptation.",
    ]
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("griis","fa79","flora82","output"):p.add_argument("--"+k,required=True)
    args=p.parse_args()
    load=lambda key:json.loads(Path(getattr(args,key)).read_text(encoding="utf-8"))
    census=requests.get(CENSUS_URL,timeout=(20,90))
    census.raise_for_status()
    states=parse_census_state_table(census.text)
    with requests.Session() as session:
        session.headers.update({"User-Agent":"tetsu-v85-2026-GBIF-2022before-after-state-source"})
        results=run(load("griis"),states,session,load("fa79"),load("flora82"))
    results["source_Census_FIPS_state_txt_SHA256"]=sha256(census.content).hexdigest()
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":results["complete"],"summary":results["summary"],
          "source_support_40":results["n_qualified_source_states_gate_40"],
          "not_first_appearance_or_evolution":True},ensure_ascii=False),flush=True)
    if not results["complete"]:
        raise SystemExit("51-state split-year GBIF source incomplete; NO taxa absence may be inferred")
if __name__=="__main__":main()
