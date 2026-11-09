#!/usr/bin/env python3
"""v1.8 reproducible GBIF macro-country x period vascular-plant coverage audit.

This uses GBIF SEARCH counts only, not downloaded occurrences. Count
queries do not hit the 100k records-retrieval ceiling, but no species,
coordinate, sampling effort or network effect is estimated here.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import time
import requests

COUNTRIES=("JP","KR","TW","US","CA","GB","DE","FR","AU","NZ")
PERIODS=((2010,2017),(2018,2025))
BASIS=(None,"HUMAN_OBSERVATION","PRESERVED_SPECIMEN")
BASE="https://api.gbif.org/v1/occurrence/search"
TAXON_SOURCE="https://api.gbif.org/v1/species/6"
USER_AGENT="tetsu-macroecology-gbif-coverage-v1.8 (country-source-counts; ecological effect not estimated)"
HARD_MAX_RETRIES=5


def validate_frame(countries=COUNTRIES,periods=PERIODS):
    if len(countries)!=len(set(countries)) or len(countries)<5:
        raise ValueError("Bad country study frame")
    if any(len(x)!=2 or len(x[0])!=2 or not x.isalpha() for x in countries):
        raise ValueError("ISO alpha-2 country frame required")
    if len(periods)!=2 or periods[0][1]+1!=periods[1][0]:
        raise ValueError("Time windows must be disjoint and contiguous")
    if periods[0][0]!=2010 or periods[1][1]!=2025:
        raise ValueError("Study window bounds changed")
    if any(a>b for a,b in periods):
        raise ValueError("Bad year interval")


def query_params(country,period,basis):
    if country not in COUNTRIES or period not in PERIODS or basis not in BASIS:
        raise ValueError("Only frozen GBIF source opportunities may be queried")
    params={
        "country":country, "kingdomKey":6, "taxonRank":"SPECIES",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "year":f"{period[0]},{period[1]}","limit":0,
    }
    if basis:
        params["basisOfRecord"]=basis
    return params


def parse_count(payload):
    if not isinstance(payload,dict):
        raise ValueError("Unexpected GBIF response type")
    n=payload.get("count")
    if isinstance(n,bool) or not isinstance(n,int) or n<0:
        raise ValueError("Missing or invalid GBIF result count")
    if payload.get("offset") not in (None,0):
        raise ValueError("Unexpected GBIF paging offset")
    if payload.get("limit") not in (None,0):
        raise ValueError("The source endpoint did not honor limit=0")
    return n


def get_json(session,url,params=None):
    for attempt in range(HARD_MAX_RETRIES):
        try:
            response=session.get(url,params=params,timeout=(20,70))
            if response.status_code in (429,500,502,503,504):
                raise requests.HTTPError(f"GBIF response {response.status_code}")
            response.raise_for_status()
            return response.json()
        except (requests.RequestException,ValueError):
            if attempt+1>=HARD_MAX_RETRIES:
                raise
            time.sleep(min(8,2**attempt))
    raise RuntimeError("GBIF request exhausted retry budget")


def audit(session,request_func=get_json,sleep_func=time.sleep):
    validate_frame()
    # Numeric backbone kingdom 6 is intentionally verified instead of
    # assuming that an updated service still maps it to Plantae.
    taxon=request_func(session,TAXON_SOURCE)
    if not isinstance(taxon,dict) or str(taxon.get("canonicalName","")).lower()!="plantae":
        raise ValueError("GBIF backbone numeric key 6 is not explicitly verified as Plantae")
    summary={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "stage":"source-country-period coverage only, no taxon incidence or railway effect",
        "gbif_taxon_key_6_verified_canonical_name":"Plantae",
        "taxon_filter":"kingdomKey=6 and species-rank (includes nonvascular plants until future Tracheophyta filter)",
        "coordinate_filter":"hasCoordinate=true and hasGeospatialIssue=false; no numeric uncertainty check yet",
        "country_frame_frozen_list":list(COUNTRIES),
        "year_windows_frozen":[list(x) for x in PERIODS],
        "source_endpoint":BASE,
        "species_level_taxon_status_not_audited":True,
        "gbif_public_search_100k_download_limit_NOT_applied_to_counts":True,
        "counts":{},
        "errors":{},
    }
    for country in COUNTRIES:
        periods=[]
        for period in PERIODS:
            row={"start_year":period[0],"end_year":period[1],"by_basis":{}}
            try:
                for basis in BASIS:
                    label=basis if basis else "ALL_BASIS"
                    params=query_params(country,period,basis)
                    row["by_basis"][label]=parse_count(
                        request_func(session,BASE,params)
                    )
                    sleep_func(.2)
                if (row["by_basis"]["HUMAN_OBSERVATION"]+
                    row["by_basis"]["PRESERVED_SPECIMEN"] >
                    row["by_basis"]["ALL_BASIS"]):
                    raise ValueError("Basis-class counts exceed all-basis records")
                periods.append(row)
                print(json.dumps({"country":country,"period":period,
                                  "counts":row["by_basis"]}),flush=True)
            except Exception as exc:
                summary["errors"][f"{country}_{period[0]}_{period[1]}"]=(
                    f"{type(exc).__name__}: {str(exc)[:180]}")
                print(json.dumps({"source_error":country,
                    "period":period,"type":type(exc).__name__}),flush=True)
        summary["counts"][country]=periods
    summary["n_successful_country_periods"]=sum(
        len(items) for items in summary["counts"].values()
    )
    summary["full_coverage"]=(
        summary["n_successful_country_periods"]==len(COUNTRIES)*len(PERIODS)
        and not summary["errors"]
    )
    summary["not_valid_for_biological_inference"]=[
        "Counts are indexed occurrence records, not deduplicated plant individuals or independently sampled survey plots.",
        "Kingdom Plantae includes non-vascular taxa; vascular plant restriction and accepted-species identification are not yet joined.",
        "Coordinate availability does not establish coordinate precision sufficient for a near-track claim.",
        "No observation-to-active-rail alignment, non-native country checklist, road-network covariates or actual establishment response is present.",
        "Differences in human observations reflect changing global citizen-science participation and cannot be treated as plant population expansion.",
        "For >100k actual source records, request an authenticated DOI-cited GBIF download; do not scrape or paginate search results past 100k.",
    ]
    return summary


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/gbif_10countries_2periods_plantae_source_counts_v18.json")
    args=p.parse_args()
    client=requests.Session()
    client.headers.update({"User-Agent":USER_AGENT,"Accept":"application/json"})
    out=audit(client)
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "full_coverage":out["full_coverage"],
        "success_country_periods":out["n_successful_country_periods"],
        "errors":out["errors"]
    },ensure_ascii=False),flush=True)
    if not out["full_coverage"]:
        raise SystemExit("Incomplete country-period source coverage; do not treat missing query as zero")


if __name__=="__main__":
    main()
