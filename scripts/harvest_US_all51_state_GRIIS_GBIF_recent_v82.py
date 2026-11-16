#!/usr/bin/env python3
"""v8.2: all 51 Census US states plus DC, 2018–25 GBIF speciesKey
positive-record source, fixed independently of 2022 freight or plant
association signs. Join Census FIPS state names, country-level USA
GRIIS accepted species keys and current GBIF stateProvince facets.

This is an indexed-source quality matrix only, not native/introduced
status validated at each state, original deduplicated specimens, true
species absence or freight/plant ecological effects.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import requests

from analyze_griis_alien_rail_source_v25 import valid_griis
from harvest_US_state_GRIIS_GBIF_plant_source_v81 import (
    CENSUS_URL,parse_census_state_table,get_state_period)
from audit_BTS_FAF6_2022_agricultural_rail_freight_v79 import VERSION

PERIOD=(2018,2025)
MIN_SOURCE_STATES=40
SOURCE_STAGE="Census FIPS US 50 states plus DC recent GBIF country-GRIIS source QA"


def run(griis,states,session,previous_faf):
    if previous_faf.get("complete") is not True or previous_faf.get("FAF_release")!=VERSION or (
         previous_faf.get("data",{}).get("n_distinct_state_code_categories_in_source")!=51):
        raise ValueError("Original 2022 FAF6 freight state source not independently verified")
    accepted=valid_griis(griis).get("US")
    if not accepted or len(accepted)<100:
        raise ValueError("United States GRIIS source country accepted plant taxonomy invalid")
    if len(states)!=51:
        raise ValueError("Original Census 50-state + DC source frame missing")
    out={
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "study":SOURCE_STAGE,
        "country":"US","date_period":"2018_2025",
        "FAF6_original_2022_state_OD_code_country_frame_verified":True,
        "US_Census_FIPS_51_states_selected_before_species_or_FAF_effect_sign":True,
        "n_requested_state_fips_codes":len(states),
        "n_GRIIS_accepted_US_country_speciesKey_categories":len(accepted),
        "states":{},"source_errors":{},
        "no_colonization_causality_native_status_or_seed_transport_assigned":True}
    for code,st in sorted(states.items()):
        try:
            result=get_state_period(session,st["name"],PERIOD,accepted)
            if not result.get("source_complete"):
                raise ValueError("No complete state source")
        except Exception as exc:
            result={"source_error":True}
            out["source_errors"][code]=type(exc).__name__+": "+str(exc)[:180]
        out["states"][code]={
            "official_Census_state_name":st["name"],
            "official_Census_state_abbreviation":st["abbreviation"],
            "GBIF_2018_2025_GRIIS_country_taxa_source":result,
        }
        print(json.dumps({"FIPS":code,"source_complete":result.get("source_complete") is True,
            "quality_eligible":result.get("source_eligible") is True,
            "n_country_GRIIS_taxa":result.get("n_source_GRIIS_country_introduced_species_categories"),
            "n_source_errors":len(out["source_errors"])},ensure_ascii=False),flush=True)
        time.sleep(.15)
    ncomplete=sum(x["GBIF_2018_2025_GRIIS_country_taxa_source"].get("source_complete") is True
                  for x in out["states"].values())
    eligible=sum(x["GBIF_2018_2025_GRIIS_country_taxa_source"].get("source_eligible") is True
                 for x in out["states"].values())
    out["complete"]=ncomplete==51 and not out["source_errors"]
    out["summary"]={
        "n_Census_FIPS_51_state_source_completed":ncomplete,
        "n_state_FIPS_with_adequate_2018_2025_GRIIS_plant_records":eligible,
        "n_state_FIPS_with_inadequate_GBIF_source":ncomplete-eligible,
        "n_state_source_errors":len(out["source_errors"]),
        "n_distinct_US_GRIIS_country_list_keys_observed_2018_2025":len({
            k for x in out["states"].values()
            for k in x["GBIF_2018_2025_GRIIS_country_taxa_source"].get(
                "source_GRIIS_speciesKey_indexed_counts",{})}),
    }
    out["source_quality_minimum_for_freight_network_comparison"]=eligible>=MIN_SOURCE_STATES
    out["limits"]=[
        "stateProvince is an interpreted GBIF text field, not Census polygon geometry; geographic completeness, synonym and abbreviation mismatches remain unverified.",
        "GRIIS USA COUNTRY-level introduced status does not establish introduction status or colonization in an individual US state; do not call non-listed taxa native.",
        "2026 GBIF indexed source speciesKey count facets do not verify independent field survey effort or species-level DETECTED; missing taxa are unrecorded, not confirmed absent.",
        "This source is 2018–25 botanical RECORDS, with one 2022 FAF6 shipment estimate. Historical sequence, freights' delivery of actual seeds, road versus rail establishment and evolution are NOT identified.",
        "FAF6 states are official economic FIPS locations, but line routing, crop contamination, species transport and disturbance exposure are all unmeasured.",
        "Source quality and country state cohort were defined independent of economic freight/botanical association coefficients.",
    ]
    return out


def main():
    ap=argparse.ArgumentParser()
    for x in ("griis","faf-v79","output"):
        ap.add_argument("--"+x,required=True)
    a=ap.parse_args()
    r=requests.get(CENSUS_URL,timeout=(20,90),
        headers={"User-Agent":"tetsu-state-flora-Census-FAF6-all51-v82"})
    r.raise_for_status()
    states=parse_census_state_table(r.text)
    aliases=json.loads(Path(a.griis).read_text(encoding="utf-8"))
    previous=json.loads(Path(a.faf_v79).read_text(encoding="utf-8"))
    with requests.Session() as session:
        session.headers.update({"User-Agent":"tetsu-US-51state-GBIF-GRIIS-indexed-source-v82"})
        out=run(aliases,states,session,previous)
    out["census_source_url"]=CENSUS_URL
    out["census_source_sha256"]=sha256(r.content).hexdigest()
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],
        "source_quality_minimum_for_freight_network_comparison":
            out["source_quality_minimum_for_freight_network_comparison"],
        "n_errors":len(out["source_errors"])},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("GBIF all-state source incomplete; do NOT manufacture plant absences")

if __name__=="__main__":
    main()
