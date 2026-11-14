#!/usr/bin/env python3
"""v8.1: public Census FIPS state identities and indexed GBIF stateProvince
speciesKey source for eight preselected U.S. states. No freight/plants
association, no verified introduced status at the state level, no seed
arrival or railway filtering. This is only a feasibility/quality pilot.

Source vintage risk: old GRIIS country-introduced GBIF Backbone keys
from v2.4 versus GBIF interpreted 2026 speciesKey facets. Unknowns never
mean species absence or native status. All original source rows remain
in external sources; only per-state source counts/quality are exported.
"""
from __future__ import annotations
import argparse,csv,io,json,math,time
from pathlib import Path
from datetime import datetime,timezone
from hashlib import sha256
from collections import Counter
import requests

from analyze_griis_alien_rail_source_v25 import valid_griis
from harvest_gbif_cell_species_facets_v20 import parse_page,FACET_PAGE
from audit_gbif_detected_original_rows_v48b import parse_status_facet

CENSUS_URL="https://www2.census.gov/geo/docs/reference/state.txt"
GBIF_URL="https://api.gbif.org/v1/occurrence/search"
PILOT_FIPS=("06","12","17","19","36","37","48","53")
SOURCE_PERIODS=((2010,2017),(2018,2025))
MAX_FACET_PAGES=12
SLEEP_BETWEEN_REQUESTS=.25
MAX_NON_DETECTED_OR_UNKNOWN=.02
MIN_DATASET_UUID_CATEGORIES=3
MIN_ALL_PLANT_INDEXED_RECORDS=100
MIN_GRIIS_COUNTRY_TAXA=5
MIN_GRIIS_COUNTRY_INDEXED_RECORDS=20
FIPS_TERRITORY_EXCLUSIONS=frozenset(("60","66","69","72","74","78"))


def parse_census_state_table(text):
    reader=csv.DictReader(io.StringIO(text),delimiter="|")
    required={"STATE","STUSAB","STATE_NAME","STATENS"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError("US Census official FIPS state code source table changed")
    states={}
    for row in reader:
        code=(row.get("STATE") or "").strip()
        name=(row.get("STATE_NAME") or "").strip()
        abbr=(row.get("STUSAB") or "").strip()
        if not code or not name or not abbr or code in states:
            raise ValueError("Invalid or duplicate official FIPS state row")
        if code not in FIPS_TERRITORY_EXCLUSIONS:
            states[code]={"name":name,"abbreviation":abbr}
    if len(states)!=51 or any(c not in states for c in PILOT_FIPS):
        raise ValueError("No exact 50-state plus DC matching census identity")
    return states


def query_params(name,period,offset=0):
    if period not in SOURCE_PERIODS or not name.strip() or (
            offset%FACET_PAGE!=0 or offset<0):
        raise ValueError("GBIF state source query outside frozen year/page bounds")
    return {
        "country":"US","stateProvince":name,
        "kingdomKey":6,"taxonRank":"SPECIES",
        "year":f"{period[0]},{period[1]}",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "limit":0,"facet":["speciesKey","occurrenceStatus","datasetKey"],
        "facetLimit":FACET_PAGE,
        "speciesKey.facetLimit":FACET_PAGE,
        "speciesKey.facetOffset":offset,
        "occurrenceStatus.facetLimit":8,
        "datasetKey.facetLimit":5,
    }


def source_get(session,params):
    for i in range(4):
        try:
            r=session.get(GBIF_URL,params=params,timeout=(25,90))
            r.raise_for_status()
            data=r.json()
            if not isinstance(data,dict):
                raise ValueError("GBIF API JSON object missing")
            return data
        except (requests.RequestException,ValueError):
            if i==3:raise
            time.sleep(i+1)


def get_state_period(session,state_name,period,keys):
    records={}
    total=None
    datasets=None
    status=None
    pages=0
    for i in range(MAX_FACET_PAGES):
        raw=source_get(session,query_params(state_name,period,i*FACET_PAGE))
        indexed,categories=parse_page(raw,"speciesKey")
        if total is None:
            total=indexed
            _,data= parse_page(raw,"datasetKey")
            datasets=len(data)
            status=parse_status_facet(raw)
            if status["indexed_source_count"]!=total:
                raise ValueError("Source indexed state status and species key facet totals disagree")
        elif indexed!=total:
            raise ValueError("Source indexed GBIF state record counts drift during pagination")
        for key,n in categories:
            if key in records:
                raise ValueError("SpeciesKey facet duplicate across pages")
            if n>0:records[key]=n
        pages+=1
        if len(categories)<FACET_PAGE:
            break
        time.sleep(SLEEP_BETWEEN_REQUESTS)
    else:
        raise ValueError("US state speciesKey source facet pagination exceeds 60000 categories")
    if total is None or sum(records.values())>total:
        raise ValueError("Source all-Plantae state species facet sum > GBIF total")
    if total==0:
        return {"source_complete":True,"source_eligible":False,
                "source_all_Plantae_indexed_records":0,
                "n_source_dataset_UUID_categories_at_least":0,
                "n_source_GRIIS_country_introduced_species_categories":0,
                "source_GRIIS_speciesKey_indexed_counts":{},
                "reason_not_qualified":"no_indexed_state_source_not_species_absence"}
    non=status["statuses"].get("NOT_DETECTED",0)+status["unclassified_count"]
    non_rate=non/total
    alien={k:n for k,n in records.items() if k in keys}
    reasons=[]
    if total<MIN_ALL_PLANT_INDEXED_RECORDS:reasons.append("below_100_record_source_gate")
    if datasets<MIN_DATASET_UUID_CATEGORIES:reasons.append("fewer_than_3_publisher_datasets")
    if non_rate>MAX_NON_DETECTED_OR_UNKNOWN:reasons.append("too_many_NOT_DETECTED_or_UNKNOWN_status")
    if len(alien)<MIN_GRIIS_COUNTRY_TAXA:reasons.append("fewer_than_5_country_GRIIS_keys")
    if sum(alien.values())<MIN_GRIIS_COUNTRY_INDEXED_RECORDS:reasons.append("fewer_than_20_country_GRIIS_indexed_records")
    return {
        "source_complete":True,"source_eligible":not reasons,
        "source_all_Plantae_indexed_records":total,
        "source_full_species_key_facet_categories":len(records),
        "n_source_dataset_UUID_categories_at_least":datasets,
        "source_NOT_DETECTED_or_unknown_fraction":round(non_rate,7),
        "n_source_GRIIS_country_introduced_species_categories":len(alien),
        "n_source_GRIIS_country_introduced_species_records":sum(alien.values()),
        "source_GRIIS_speciesKey_indexed_counts":alien,
        "source_facets_fully_paginated_n_pages":pages,
        "reason_not_qualified":reasons,
        "GRIIS_country_level_introduced_DOES_NOT_imply_state_level_introduced":True,
    }


def run(checklist,census_text,session):
    country=valid_griis(checklist).get("US")
    if not country or len(country)<100:
        raise ValueError("Original GRIIS source USA country introduced species list not verified")
    census=parse_census_state_table(census_text)
    out={
        "complete":True,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "study":"2026 GBIF indexed source US stateProvince pilot for GRIIS USA-country introduced taxa and FAF6 2022 state FIPS geography",
        "pilot_FIPS_state_codes_chosen_without_any_plant_or_freight_association_sign":list(PILOT_FIPS),
        "pilot_geography":"8 distinct preselected states spanning US regions, Census FIPS code state name map",
        "source_year_periods":[list(p) for p in SOURCE_PERIODS],
        "source_census_FIPS_SHA256":sha256(census_text.encode()).hexdigest(),
        "n_US_states_plus_DC_in_official_Census_source":len(census),
        "n_valid_GRIIS_country_US_plant_keys":len(country),
        "source_GBIF_status_filter_legacy_PRESENT_NOT_used":True,
        "states":{},"errors":{},
        "NO_state_level_nativity_or_causal_railway_plant_inference":True
    }
    for code in PILOT_FIPS:
        state=census[code]
        result={"name":state["name"],"abbreviation":state["abbreviation"],"periods":{}}
        for period in SOURCE_PERIODS:
            key=f"{period[0]}_{period[1]}"
            try:
                result["periods"][key]=get_state_period(session,state["name"],period,country)
            except Exception as e:
                out["errors"][code+"_"+key]=str(type(e).__name__)+": "+str(e)[:190]
                result["periods"][key]={"source_error":True}
            time.sleep(SLEEP_BETWEEN_REQUESTS)
        result["source_qualified_in_BOTH_year_windows"]=all(
            result["periods"][f"{a}_{b}"].get("source_eligible") is True
            for a,b in SOURCE_PERIODS)
        out["states"][code]=result
        print(json.dumps({
            "FIPS":code,"state":state["name"],
            "eligible_both_periods":result["source_qualified_in_BOTH_year_windows"],
            "late_GRIIS_taxa":result["periods"]["2018_2025"].get(
                "n_source_GRIIS_country_introduced_species_categories"),
            "source_errors":len(out["errors"])},ensure_ascii=False),flush=True)
    out["complete"]=not out["errors"] and len(out["states"])==len(PILOT_FIPS)
    out["summary"]={
        "n_states_source_complete":sum(all(x.get("source_complete") is True for x in v["periods"].values())
                                          for v in out["states"].values()),
        "n_US_states_two_period_source_qualified":sum(v["source_qualified_in_BOTH_year_windows"]
                                                       for v in out["states"].values()),
        "n_US_states_with_recent_at_least_5_country_GRIIS_taxa":
            sum(v["periods"]["2018_2025"].get("n_source_GRIIS_country_introduced_species_categories",0)>=5
                for v in out["states"].values())}
    out["limitations"]=[
        "stateProvince is a GBIF interpreted source string, NOT full Census-polygon occurrence inclusion; state-name/abbreviation inconsistencies and missing source locations can bias the results.",
        "GRIIS USA-country naturalized/taxonomic status does not verify INTRODUCED status in each STATE.",
        "Legacy GBIF accepted speciesKey versus modern 2026 interpreted species taxonomy may drift; missing keys never mean plant absences.",
        "Indexed GBIF positive speciesKey facets are not original deduplicated observer/plot data, and indexed source DETECTED status is not separable species-by-species.",
        "2022 FAF6 economic OD cannot confirm colonization after freight or active rail service; freight and plant data have different time/spatial grains.",
        "This is an 8-state SOURCE FEASIBILITY pilot; no chosen ecological association has been fitted or tested."
    ]
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--griis",required=True);ap.add_argument("--output",required=True)
    a=ap.parse_args()
    aliases=json.loads(Path(a.griis).read_text(encoding="utf-8"))
    r=requests.get(CENSUS_URL,timeout=(20,60),headers={"User-Agent":"tetsu-census-US-FIPS-state-GBIF-plant-v81"})
    r.raise_for_status()
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-Census-FAF6-GBIF-state-source-v81"})
    out=run(aliases,r.text,session)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],"errors":out["errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Not all preselected US states returned valid source data; DO NOT infer biological absence")

if __name__=="__main__":main()
