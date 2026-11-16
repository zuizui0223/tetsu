#!/usr/bin/env python3
"""v4.8 bounded GBIF ORIGINAL OCCURRENCE-ROW metadata quality audit.

Query 6 original country macro source blocks × EXACT 2 pre-fixed 50km
cell positions x1y1 and x2y2 × BOTH periods (2010-17, 2018-25) from
GBIF occurrence/search API. Retrieve only FIRST 300 indexed occurrence
rows per bin: bounded, NONRANDOM source retrieval. It is NOT a
DOI-backed full occurrence download or representative audit of all
plant record rows. The official service has 300 records/page and a
100000-record search offset+limit ceiling.

Aggregate and store NO raw occurrence IDs, precise coordinates or
individual recordedBy names. Calculate only real fetched-row coverage:
dataset concentration, missing recorder, missing/large coordinates
uncertainty, missing event date, same rounded point/species/date proxy
collisions and country/50km query spatial consistency.

Do not convert first-page selection into survey effort correction,
species absence, seed occurrence at a railway, or national population
rates. Requires authentic API; fake tests are not biological evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time

import requests
from shapely.geometry import Point
from shapely import wkt

from audit_gbif_equalarea_macroblocks_v29 import make_grid, BLOCKS
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX, PERIODS

API="https://api.gbif.org/v1/occurrence/search"
FIXED_OFFSETS=((1,1),(2,2))
PAGE_LIMIT=300
MIN_GENUINE_SOURCE_RECORDS_PER_BIN=0
MAX_QUERY_RETRIES=4


def frozen_cells(block):
    if block not in BLOCKS_SIX:
        raise ValueError("Source geographic blocks have changed")
    original={r["cell_id"]:r for r in make_grid(block)}
    names=[f"{block}_x{x}_y{y}" for x,y in FIXED_OFFSETS]
    if any(n not in original for n in names):
        raise ValueError("Missing fixed 50km target source cells")
    return [original[name] for name in names]


def source_params(cell,period,limit=PAGE_LIMIT):
    if period not in PERIODS:
        raise ValueError("Unfrozen historical source period")
    if isinstance(limit,bool) or not isinstance(limit,int) or not(1<=limit<=300):
        raise ValueError("GBIF endpoint does not allow >300 source rows per page")
    return {
        "country":cell["country"],
        "geometry":cell["gbif_geometry_wgs84"],
        "kingdomKey":6,
        "taxonRank":"SPECIES",
        "occurrenceStatus":"PRESENT",
        "hasCoordinate":"true",
        "hasGeospatialIssue":"false",
        "year":period.replace("_",","),
        "offset":0,
        "limit":limit,
    }


def source_page(session,params):
    for retry in range(MAX_QUERY_RETRIES):
        try:
            response=session.get(API,params=params,timeout=(20,90))
            response.raise_for_status()
            data=response.json()
            if not isinstance(data,dict):
                raise ValueError("Non-JSON-object GBIF source occurrence page")
            return data
        except (ValueError,requests.RequestException):
            if retry+1==MAX_QUERY_RETRIES:
                raise
            time.sleep(2**retry)
    raise AssertionError("Unreachable request retry branch")


def is_numeric_positive(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) and x>=0


def audit_source_page(cell,period,payload,prior_flora=None):
    if period not in PERIODS:
        raise ValueError("Unexpected source time period")
    total=payload.get("count")
    records=payload.get("results")
    if isinstance(total,bool) or not isinstance(total,int) or total<0:
        raise ValueError("Invalid full-source count")
    if not isinstance(records,list) or len(records)>PAGE_LIMIT:
        raise ValueError("Invalid GBIF 300-record source page")
    if len(records)>total:
        raise ValueError("Source page exceeds indexed count")
    expected=min(total,PAGE_LIMIT)
    if len(records)!=expected:
        raise ValueError(f"Unexpected incomplete source first page {len(records)} vs {expected}")
    if total<=PAGE_LIMIT and payload.get("endOfRecords") is not True:
        raise ValueError("Source claims last partial page but not end of records")
    source_poly=wkt.loads(cell["gbif_geometry_wgs84"])
    if not source_poly.is_valid or source_poly.is_empty:
        raise ValueError("Fixed GBIF source macro cell polygon invalid")
    dataset=Counter()
    basis=Counter()
    issues=Counter()
    identifier_keys=set()
    duplicate_gbif_ids=0
    missing_dataset=0
    missing_recorded_by=0
    missing_event_date=0
    missing_occurrence_id=0
    missing_event_id=0
    missing_coordinate_uncertainty=0
    source_uncertainty_gt1000m=0
    source_uncertainty_gt25000m=0
    source_uncertainty_gt50000m=0
    bad_coordinate=0
    coordinate_outside_query_polygon=0
    invalid_source_species_rank=0
    invalid_source_year=0
    missing_species_key=0
    repeated_dataset_occurrence_keys=0
    repeated_round_coord_species_date=0
    occurrence_keys=set()
    coarse_keys=set()
    for row in records:
        if not isinstance(row,dict):
            raise ValueError("Non-dictionary actual source GBIF record")
        g=row.get("key")
        if isinstance(g,bool) or not isinstance(g,int):
            raise ValueError("GBIF occurrence row lacks indexed numeric key")
        if g in identifier_keys:
            duplicate_gbif_ids+=1
        identifier_keys.add(g)
        data_key=str(row.get("datasetKey") or "").strip()
        if not data_key:
            missing_dataset+=1
            data_key="__MISSING_DATASET__"
        dataset[data_key]+=1
        basis[str(row.get("basisOfRecord") or "UNKNOWN")]+=1
        for issue in row.get("issues",[]) or []:
            issues[str(issue)]+=1
        recorded_by=row.get("recordedBy")
        if not isinstance(recorded_by,str) or not recorded_by.strip():
            missing_recorded_by+=1
        date=str(row.get("eventDate") or "").strip()
        if not date:
            missing_event_date+=1
        if not str(row.get("eventID") or "").strip():
            missing_event_id+=1
        source_occid=str(row.get("occurrenceID") or "").strip()
        if not source_occid:
            missing_occurrence_id+=1
        else:
            key=(data_key,source_occid)
            if key in occurrence_keys:
                repeated_dataset_occurrence_keys+=1
            occurrence_keys.add(key)
        uncertainty=row.get("coordinateUncertaintyInMeters")
        if not is_numeric_positive(uncertainty):
            missing_coordinate_uncertainty+=1
        else:
            if uncertainty>1000:source_uncertainty_gt1000m+=1
            if uncertainty>25000:source_uncertainty_gt25000m+=1
            if uncertainty>50000:source_uncertainty_gt50000m+=1
        lon=row.get("decimalLongitude")
        lat=row.get("decimalLatitude")
        if not (is_numeric_positive(abs(lon)) if isinstance(lon,(int,float)) and not isinstance(lon,bool) else False) or not (
            is_numeric_positive(abs(lat)) if isinstance(lat,(int,float)) and not isinstance(lat,bool) else False):
            bad_coordinate+=1
            continue
        if not (-180<=lon<=180 and -90<=lat<=90):
            bad_coordinate+=1
            continue
        if not source_poly.covers(Point(lon,lat)):
            coordinate_outside_query_polygon+=1
        species=row.get("speciesKey")
        if isinstance(species,bool) or not isinstance(species,int):
            missing_species_key+=1
        if str(row.get("rank") or row.get("taxonRank") or "").upper()!="SPECIES":
            invalid_source_species_rank+=1
        yr=row.get("year")
        a,b=(int(x) for x in period.split("_"))
        if isinstance(yr,bool) or not isinstance(yr,int) or not a<=yr<=b:
            invalid_source_year+=1
        if species is not None and date:
            coarse=(species,date,round(lon,4),round(lat,4))
            if coarse in coarse_keys:
                repeated_round_coord_species_date+=1
            coarse_keys.add(coarse)
    n=len(records)
    hhi=sum((v/n)**2 for v in dataset.values()) if n else None
    results={
        "source_country":cell["country"],
        "source_cell_id":cell["cell_id"],
        "period":period,
        "indexed_source_search_total":total,
        "actual_original_record_rows_examined":n,
        "actual_search_page_offset":0,
        "max_GBIF_search_page_size":PAGE_LIMIT,
        "records_exhaustively_scanned":total<=PAGE_LIMIT,
        "first_page_is_nonrandom_if_more_records":total>PAGE_LIMIT,
        "fraction_of_indexed_records_returned":round(n/total,7) if total else None,
        "has_valid_GBIF_indexed_keys":duplicate_gbif_ids==0,
        "unique_indexed_GBIF_keys":len(identifier_keys),
        "n_missing_dataset_key":missing_dataset,
        "n_distinct_source_datasets_in_first_page":len(dataset),
        "top_dataset_record_share_in_first_page":round(max(dataset.values())/n,6) if n else None,
        "first_page_source_dataset_HHI":round(hhi,6) if hhi is not None else None,
        "basis_of_record_first_page":dict(basis),
        "n_missing_recordedBy_first_page":missing_recorded_by,
        "n_missing_eventDate_first_page":missing_event_date,
        "n_missing_eventID_first_page":missing_event_id,
        "n_missing_occurrenceID_first_page":missing_occurrence_id,
        "n_missing_coordinateUncertaintyInMeters_first_page":missing_coordinate_uncertainty,
        "n_coordinate_uncertainty_gt1000m_first_page":source_uncertainty_gt1000m,
        "n_coordinate_uncertainty_gt25000m_first_page":source_uncertainty_gt25000m,
        "n_coordinate_uncertainty_gt50000m_first_page":source_uncertainty_gt50000m,
        "n_invalid_or_missing_coordinates_first_page":bad_coordinate,
        "n_coords_outside_original_query_polygon_first_page":coordinate_outside_query_polygon,
        "n_wrong_rank_first_page":invalid_source_species_rank,
        "n_wrong_historical_year_first_page":invalid_source_year,
        "n_missing_speciesKey_first_page":missing_species_key,
        "n_duplicate_dataset_occurrenceID_first_page":repeated_dataset_occurrence_keys,
        "n_same_rounded_point_species_eventdate_collisions_first_page":repeated_round_coord_species_date,
        "top_GBIF_index_issues_first_page":issues.most_common(10),
        "raw_row_ID_coordinates_and_human_recorder_names_NOT_saved":True,
        "not_an_independent_random_survey_sample":True,
    }
    if isinstance(prior_flora,dict):
        c=prior_flora.get("periods",{}).get(period,{})
        results["prior_source_GRIIS_qualified_for_ecological_model"]=bool(
            c.get("alien_observed_minimum_5_species_and_20_records") and
            prior_flora.get("source_pre_outcome_both_period_recording_eligible")
        )
        prior=c.get("indexed_all_plantae_source_records")
        results["prior_indexed_source_record_count"]=prior
        results["observed_source_count_index_drift"]=total-prior if isinstance(prior,int) else None
    return results


def run(flora,query=source_page,sleep=time.sleep):
    if not flora.get("complete") or flora.get("source_errors"):
        raise ValueError("Original GRIIS 50km source was incomplete")
    if not set(BLOCKS_SIX).issubset(flora.get("regions",{})):
        raise ValueError("Six geographic region source frame missing")
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-macro-GBIF-original-row-bounded-metadata-audit-v4.8"})
    output={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_api":API,
        "scope":"first 300 GBIF original occurrence rows for 2 fixed cells per country/period, not a full dataset download",
        "original_country_blocks":list(BLOCKS_SIX),
        "frozen_cell_offsets":[list(t) for t in FIXED_OFFSETS],
        "regions":{},"errors":{},
        "not_DOI_citable_bulk_occurrence_download":True,
        "not_a_random_sample_of_original_source_records":True,
    }
    for block in BLOCKS_SIX:
        orig=flora["regions"][block]
        if orig.get("country")!=BLOCKS[block][0]:
            raise ValueError("GBIF source country block mismatches")
        rows={}
        for cell in frozen_cells(block):
            cid=cell["cell_id"]
            results={}
            for p in PERIODS:
                try:
                    query_response=query(session,source_params(cell,p))
                    results[p]=audit_source_page(cell,p,query_response,
                                                 prior_flora=orig["cells"].get(cid))
                except Exception as exc:
                    output["errors"][cid+"_"+p]=f"{type(exc).__name__}: {str(exc)[:240]}"
                    results[p]={"source_error":True}
                sleep(.1)
            rows[cid]=results
        output["regions"][block]={"country":BLOCKS[block][0],
            "n_fixed_50km_cell_origins":len(rows),
            "n_cell_period_pages_successful":sum(
                "actual_original_record_rows_examined" in p for row in rows.values() for p in row.values()),
            "cells":rows}
        totals=[x for row in rows.values() for x in row.values() if "actual_original_record_rows_examined" in x]
        print(json.dumps({"block":block,"country":BLOCKS[block][0],
            "n_successful_original_row_bins":len(totals),
            "n_actual_original_rows_inspected":sum(x["actual_original_record_rows_examined"] for x in totals),
            "missing_coord_uncertainty":sum(x["n_missing_coordinateUncertaintyInMeters_first_page"] for x in totals),
            "missing_recorded_by":sum(x["n_missing_recordedBy_first_page"] for x in totals),
            "bins_nonrandom_first_page":sum(x["first_page_is_nonrandom_if_more_records"] for x in totals),
            "source_errors":len(output["errors"])
        },ensure_ascii=False),flush=True)
    output["complete"]=len(output["regions"])==len(BLOCKS_SIX) and not output["errors"]
    output["summary"]={
        "n_country_blocks":len(output["regions"]),
        "n_original_record_source_bins":len(BLOCKS_SIX)*len(FIXED_OFFSETS)*len(PERIODS),
        "n_bins_successful":sum(x["n_cell_period_pages_successful"] for x in output["regions"].values()),
        "n_original_row_records_inspected":sum(
            p.get("actual_original_record_rows_examined",0) for b in output["regions"].values()
            for c in b["cells"].values() for p in c.values()),
        "n_nonrandom_first_pages":sum(
            p.get("first_page_is_nonrandom_if_more_records",False) for b in output["regions"].values()
            for c in b["cells"].values() for p in c.values()),
    }
    output["limits"]=[
        "Only FIRST page (<=300) of GBIF search results in each cell and period, relevance/order nonrandom; no population proportions or unbiased quality estimates possible where more than 300 records exist.",
        "Record counts are standard kingdom Plantae and species rank, NOT restricted to country-GRIIS listed introduced flora.",
        "No full DOI-citable GBIF download request was made; authenticated bulk user credential is needed for event-level complete data and citations.",
        "Missing coordinateUncertaintyInMeters does not prove the true position is inaccurate; issue-free GBIF coords do not mean coordinate uncertainty is small.",
        "Matching eventDate/species/rounded coordinates can include legitimately multiple records; collision is NOT verified duplicate.",
        "Original human recorder names, exact lat/lon and record ID keys are not saved, published or committed.",
        "No actual rail-adjacent specimen or freight/seed transfer outcome is audited by this source QA.",
    ]
    return output


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--griis-flora-source",required=True)
    ap.add_argument("--output",default="results/audit/GBIF_original_firstpage_6country_metadata_v48.json")
    a=ap.parse_args()
    flora=json.loads(Path(a.griis_flora_source).read_text(encoding="utf-8"))
    out=run(flora)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],
                      "errors":out["errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("GBIF occurrence original-row diagnostic incomplete; cannot claim global event quality")


if __name__=="__main__":
    main()
