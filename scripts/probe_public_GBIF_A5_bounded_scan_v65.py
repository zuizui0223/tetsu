#!/usr/bin/env python3
"""v6.5 independent, bounded US/FR GBIF public August-2026 original A5 file audit.

Why: v6.4 France failed because the first FIVE LEXICOGRAPHIC Parquet
objects contained zero French-Burgundy source rows; it was incorrect to
interpret this as absence of botanical original events. BEFORE reading
country-GRIIS plant outcome signs, re-query a prespecified bounded
set of up to 48 actual source Parquet objects in the preexisting
France-Burgundy A5 centroid partitions.

Unlike v6.0, do not cherry-pick the FIRST positive object: scan ALL
selected source objects up to an explicit resource budget and report
zero and positive source object counts, row totals, and completeness
of this selected prefix. This remains a DEVELOPMENT PREVIEW derivative
of GBIF's August-2026 original archive DOI 10.15468/dl.8yrbe7,
NOT a newly authenticated DOI download and NOT a complete spatial
sample. A5 source centroid partitions can miss boundary overlap.

Only aggregate source-row counts from real GBIF interpreted events,
never raw GBIF record IDs or coordinate rows. No seed transport effect.
"""
from __future__ import annotations

import argparse,json,re
from datetime import datetime,timezone
from pathlib import Path

from probe_public_GBIF_A5_plant_events_v60 import (
    BUCKET,ROOT,UNDERLYING_GBIF_DOI,MAX_ONE_PARQUET_OBJECT_BYTES,
    FROZEN_BLOCKS,listing,source_a5_cells,
    find_small_geographic_plant_files,read_true_public_parquet_source
)

MAX_BOUNDED_SOURCE_OBJECTS=48
MIN_BOUNDED_SOURCE_OBJECTS=6
RESOURCE_LIMITATION_NOT_FREQUENTIST_SAMPLING=True


def scan_bounded_source(db,block,candidates,query=read_true_public_parquet_source,
                        scan_limit=MAX_BOUNDED_SOURCE_OBJECTS):
    if block not in FROZEN_BLOCKS:
        raise ValueError("Block not in original frozen geography")
    if not isinstance(candidates,list) or len(candidates)<MIN_BOUNDED_SOURCE_OBJECTS:
        raise ValueError("Insufficient public Parquet candidate source objects for v6.4 repair")
    if isinstance(scan_limit,bool) or not isinstance(scan_limit,int) or scan_limit<MIN_BOUNDED_SOURCE_OBJECTS or scan_limit>MAX_BOUNDED_SOURCE_OBJECTS:
        raise ValueError("Bad source-bounded pre-outcome file selection count")
    filenames=[p.get("key") for p in candidates]
    if any(not isinstance(key,str) or not key.startswith(ROOT) for key in filenames) or len(set(filenames))!=len(filenames):
        raise ValueError("Source candidate object key prefix/uniqueness wrong")
    # Natural-number ordering is frozen independently of botanical
    # outcomes: lexical sort data_1, data_11, ... before data_2
    # silently oversamples source files with two-digit suffixes.
    def source_sort_key(item):
        key=item["key"]
        leaf=key.rsplit("/",1)[-1]
        m=re.fullmatch(r"data_(\\d+)\\.parquet",leaf)
        if m:
            return (key.rsplit("/",1)[0],0,int(m.group(1)),key)
        return (key.rsplit("/",1)[0],1,0,key)
    selected=sorted(candidates,key=source_sort_key)[:scan_limit]
    observed=[]
    failed=[]
    for candidate in selected:
        try:
            row=query(db,block,candidate)
            n=row.get("n_actual_original_interpreted_rows_in_pilot_bbox")
            if not isinstance(n,int) or isinstance(n,bool) or n<0:
                raise ValueError("Actual original source row count missing")
            observed.append({
                "source_key":row["selected_public_parquet_key"],
                "source_size_bytes":row["selected_public_parquet_size_bytes"],
                "original_interpreted_source_rows_in_200km_bounding_envelope":n,
                "n_source_Present_status":row["n_Present_status_source_rows"],
                "n_source_Absent_status":row["n_Absent_status_source_rows"],
                "n_source_missing_coord_uncertainty":row["n_missing_coordinate_uncertainty_source_rows"],
                "n_source_datasets_within_this_file":row["n_distinct_GBIF_dataset_UUIDs_in_sample"],
                "actual_source_rows_read":True,
            })
            print(json.dumps({"block":block,
                "bounded_source_object_file":candidate["key"].split("/")[-1],
                "source_rows_in_region":n,
                "source_present":row["n_Present_status_source_rows"]
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            failed.append({"source_key":candidate["key"],
                           "source_error":f"{type(exc).__name__}: {str(exc)[:180]}"})
            print(json.dumps({"block":block,"source_query_failure":candidate["key"],
                              "error":failed[-1]["source_error"]},ensure_ascii=False),flush=True)
    total=sum(x["original_interpreted_source_rows_in_200km_bounding_envelope"] for x in observed)
    stats={
        "block":block,
        "source_candidate_objects_in_centroid_A5_bound":len(candidates),
        "n_source_objects_selected_by_natural_numeric_filename_BEFORE_reading_species_outcomes":len(selected),
        "n_source_objects_actually_scanned":len(observed),
        "n_source_objects_failed_API_or_IO":len(failed),
        "n_source_objects_with_nonzero_actual_source_occurrences":sum(x["original_interpreted_source_rows_in_200km_bounding_envelope"]>0 for x in observed),
        "n_source_objects_with_zero_actual_source_occurrences":sum(x["original_interpreted_source_rows_in_200km_bounding_envelope"]==0 for x in observed),
        "n_actual_interpreted_event_source_rows_in_selected_files_and_bbox":total,
        "n_Present_records_in_selected_files":sum(x["n_source_Present_status"] for x in observed),
        "n_Absent_records_in_selected_files":sum(x["n_source_Absent_status"] for x in observed),
        "n_missing_coord_uncertainty_in_selected_source_files":sum(x["n_source_missing_coord_uncertainty"] for x in observed),
        "full_centroid_partition_source_object_list_covered":len(candidates)<=scan_limit and not failed,
        "all_200km_original_16_cell_geography_fully_observed":False,
        "source_object_results":observed,
        "source_errors":failed,
        "genuine_GBIF_original_interpreted_rows_verified":total>0 and not failed,
        "NO_country_GRIIS_species_composition_or_railway_effect_estimated":True,
    }
    return stats


def run(blocks,query=read_true_public_parquet_source,scan_limit=MAX_BOUNDED_SOURCE_OBJECTS):
    import boto3,duckdb
    from botocore import UNSIGNED
    from botocore.config import Config
    if not blocks or len(set(blocks))!=len(blocks) or any(x not in FROZEN_BLOCKS for x in blocks):
        raise ValueError("Source geography must be original US Midwest and/or France Burgundy")
    s3=boto3.client("s3",region_name="us-east-1",config=Config(signature_version=UNSIGNED))
    db=duckdb.connect()
    db.execute("SET memory_limit='3GB'; SET threads=2")
    db.execute("INSTALL httpfs; LOAD httpfs;")
    db.execute("SET s3_region='us-east-1';")
    db.execute("INSTALL a5 FROM community; LOAD a5;")
    result={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source":"GBIF experimental public August2026 original interpreted Parquet A5",
        "underlying_original_GBIF_DWCA_DOI":UNDERLYING_GBIF_DOI,
        "geography_frozen_before_any_new_flora_outcome":list(blocks),
        "max_source_objects_per_block_predeclared":scan_limit,
        "regions":{},"errors":{},
    }
    try:
        for block in blocks:
            try:
                original_centroids=source_a5_cells(db,block)
                groups,candidates=find_small_geographic_plant_files(s3,original_centroids)
                if len(candidates)<MIN_BOUNDED_SOURCE_OBJECTS:
                    raise ValueError("No adequate original public source object listing for expanded bounded scan")
                d=scan_bounded_source(db,block,candidates,query,scan_limit)
                d["a5_r2_source_centroid_partitions"]=original_centroids
                d["plant_taxon_partition_prefixes"]=groups
                result["regions"][block]=d
                print(json.dumps({"block":block,
                    "candidate_source_files":d["source_candidate_objects_in_centroid_A5_bound"],
                    "scanned_source_files":d["n_source_objects_actually_scanned"],
                    "positive_source_files":d["n_source_objects_with_nonzero_actual_source_occurrences"],
                    "actual_source_rows_in_selected_files":d["n_actual_interpreted_event_source_rows_in_selected_files_and_bbox"],
                    "source_preview_not_exhaustive":True
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                result["errors"][block]=f"{type(exc).__name__}: {str(exc)[:260]}"
                print(json.dumps({"source_error":block,"detail":result["errors"][block]},
                                 ensure_ascii=False),flush=True)
    finally:
        db.close()
    result["complete"]=not result["errors"] and len(result["regions"])==len(blocks) and all(
        x["genuine_GBIF_original_interpreted_rows_verified"] for x in result["regions"].values())
    result["limits"]=[
        "The original 50km A5 source partitions are identified using cell CENTROIDS, not complete coverage of cell boundary A5 tiles.",
        "The source A5 snapshot is an August2026 EXPERIMENTAL preview; it is not a new DOI authenticated download or assuredly complete country population.",
        "Candidate objects are only <=400MB due to original source listing cap; larger files and unselected A5 partitions may contain relevant events.",
        "Every scanned Parquet is reported, including zeros and failed files; selection by prespecified natural numeric file order and a source resource budget is not random representative sampling.",
        "The source aggregation is a COUNTRY 200km BOUNDING BOX, not necessarily exact original sixteen EPSG6933 geographic cell membership.",
        "No country GRIIS alien species status, original row dedup, railway active freight service, or seed transport effects are calculated in this acquisition QA.",
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--blocks",nargs="+",choices=FROZEN_BLOCKS,default=["fr_burgundy"])
    p.add_argument("--scan-limit",type=int,default=MAX_BOUNDED_SOURCE_OBJECTS)
    p.add_argument("--output",default="results/audit/GBIF_public_A5_original_bounded_scan_v65.json")
    a=p.parse_args()
    r=run(a.blocks,scan_limit=a.scan_limit)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":r["complete"],"source_block_count":len(r["regions"]),
                      "errors":r["errors"]},ensure_ascii=False),flush=True)
    if not r["complete"]:
        raise SystemExit("Bounded original GBIF experimental source does not yield validated original event rows")


if __name__=="__main__":
    main()
