#!/usr/bin/env python3
"""v6.8 source-only audit of 33 previously UNREAD GBIF public Parquet objects.

Official GBIF August2026 experimental original A5 (underlying source
DWCA DOI 10.15468/dl.8yrbe7) has more France Burgundy source objects
than v6.5's preselected first 48 files. v6.7 compared original 16
50km cell polygon A5 vs centroid cell IDs: their union lists 81 public
Plantae Parquet objects, 33 unseen and 10 >400MB.

BEFORE any botanical GRIIS taxon or railway effect inspection, query
ONLY these 33 unscanned objects for true ALL PLANTAE original indexed
row counts in the identical source France Burgundy 200km bounding
envelope and 2018–2025 time period. Compare to the verified 671,528
rows from previously scanned 48 files. All error/missing large source
objects are reported, never converted into botanical absence or zero.

No individual GBIF IDs, original observation coordinates, observer
names or species-specific GRIIS responses exported. No new DOI, full
population, active freight or train-carried seed effect inferred.
"""
from __future__ import annotations

import argparse,json,re
from datetime import datetime,timezone
from pathlib import Path

from audit_public_GBIF_FR_A5_polygon_coverage_v67 import (
    BLOCK,ROOT,BUCKET,polygon_cells,list_all_source_parquets
)
from probe_public_GBIF_A5_plant_events_v60 import (
    listing,read_true_public_parquet_source,UNDERLYING_GBIF_DOI
)

SOURCE_ORIGINAL_FILES=48
EXPECTED_ADDITIONAL_FILES=33
SOURCE_TOTAL_FILES=81
MAX_SINGLE_PARQUET_BYTES=1_000_000_000


def classify_files(previous,geometry,all_files):
    if previous.get("complete") is not True or set(previous.get("regions",{}))!={BLOCK}:
        raise ValueError("Original genuine v6.5 file-level source audit missing")
    prior=previous["regions"][BLOCK]
    if prior.get("n_source_objects_actually_scanned")!=SOURCE_ORIGINAL_FILES or (
            prior.get("n_actual_interpreted_event_source_rows_in_selected_files_and_bbox")!=671528):
        raise ValueError("Precommitted original GBIF France 48-file + denominator contract changed")
    if geometry.get("original_region")!=BLOCK or geometry.get("n_total_source_parquet_files_in_all_true_A5_polygon_intersecting_partitions")!=SOURCE_TOTAL_FILES:
        raise ValueError("Independent v6.7 A5 polygon+centroid source census not 81 files")
    recorded={x.get("source_key") for x in prior.get("source_object_results",[])}
    bykey={x.get("key"):x for x in all_files}
    if len(recorded)!=SOURCE_ORIGINAL_FILES or len(bykey)!=len(all_files) or len(bykey)!=SOURCE_TOTAL_FILES:
        raise ValueError("Incomplete source Parquet object keys or duplicate S3 data")
    if not recorded.issubset(bykey):
        raise ValueError("48 original GBIF Parquet objects missing from same-source expanded A5 listing")
    skipped=[bykey[k] for k in sorted(bykey.keys()-recorded)]
    if len(skipped)!=EXPECTED_ADDITIONAL_FILES:
        raise ValueError("Missing file count changes source completion design")
    for v in skipped:
        if v["size_bytes"]>MAX_SINGLE_PARQUET_BYTES or v["size_bytes"]<=0:
            raise ValueError("GBIF source Parquet exceeded explicit 1GB object handling budget; no absent inference")
        match=re.search(r"/a5_r2=(\d+)/",v["key"])
        if not match:
            raise ValueError("A5 source object lacks source partition ID")
        v["a5_r2"]=int(match.group(1))
        v["taxon_partition"]="taxonpartition=Plantae"
    return skipped


def scan_unread(db,previous,geommeta,sourcefiles,query=read_true_public_parquet_source):
    missing=classify_files(previous,geommeta,sourcefiles)
    done=[];errors=[]
    for item in missing:
        try:
            source=query(db,BLOCK,item)
            count=source.get("n_actual_original_interpreted_rows_in_pilot_bbox")
            if isinstance(count,bool) or not isinstance(count,int) or count<0:
                raise ValueError("Missing actual GBIF original Parquet country-year source count")
            done.append({
                "source_Parquet_key":item["key"],
                "source_compressed_bytes":item["size_bytes"],
                "n_original_Plantae_rows_in_2018to2025_FR_200km_bbox":count,
                "n_rowgrain_PRESENT_status":source["n_Present_status_source_rows"],
                "n_rowgrain_ABSENT_status":source["n_Absent_status_source_rows"],
                "n_missing_stated_coordinate_uncertainty":source["n_missing_coordinate_uncertainty_source_rows"]
            })
            print(json.dumps({"block":BLOCK,"source_object_file":item["key"].rsplit("/",1)[-1],
                "original_new_source_rows":count,
                "source_anomaly_note":"This file was not in original 48-file source sample",
                "not_GRIIS_alien_count":True
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            errors.append({"source_Parquet_key":item["key"],
                           "error":f"{type(exc).__name__}: {str(exc)[:240]}"})
            print(json.dumps({"block":BLOCK,"unread_Parquet_source_error":errors[-1]},ensure_ascii=False),flush=True)
    total=sum(r["n_original_Plantae_rows_in_2018to2025_FR_200km_bbox"] for r in done)
    previous_rows=previous["regions"][BLOCK]["n_actual_interpreted_event_source_rows_in_selected_files_and_bbox"]
    return {
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_country":"FR",
        "block":BLOCK,
        "original_GBIF_August2026_source_DOI":UNDERLYING_GBIF_DOI,
        "source_EXPERIMENTAL_Parquet_preview":True,
        "source_Parquet_metadata_total_in_A5_polygon_centroid_UNION":SOURCE_TOTAL_FILES,
        "previously_scanned_source_Parquet_objects":SOURCE_ORIGINAL_FILES,
        "new_additional_source_Parquet_objects_planned":EXPECTED_ADDITIONAL_FILES,
        "n_new_additional_source_Parquet_objects_successful":len(done),
        "n_new_additional_source_Parquet_objects_with_at_least_one_source_event":sum(
            r["n_original_Plantae_rows_in_2018to2025_FR_200km_bbox"]>0 for r in done),
        "2018to2025_raw_all_Plantae_original_rows_from_old_48_files":previous_rows,
        "2018to2025_raw_all_Plantae_original_rows_from_new_33_files":total,
        "2018to2025_combined_81_file_raw_all_Plantae_source_rows_in_same_bbox":
            previous_rows+total if not errors else None,
        "fraction_of_combined_bounded_81_file_rows_covered_by_old_48":
            round(previous_rows/(previous_rows+total),7) if not errors and previous_rows+total>0 else None,
        "n_source_files_unread_or_invalid":len(errors),
        "source_additional_object_count_diagnostics":done,
        "source_query_errors":errors,
        "complete":not errors and len(done)==EXPECTED_ADDITIONAL_FILES,
        "all_source_country_GRIIS_taxon_outcomes_UNOBSERVED":True,
        "not_a_globally_COMPLETE_GBIF_census_or_new_DOI":True,
        "original_source_50km_true_A5_polygon_cover_DISAGREED_with_centroid":geommeta.get(
            "source_50km_polygon_plus_centroid_A5_grid_audit",{}).get("polygon_cover_disagreed_with_centroid_lookup"),
        "limitations":[
            "Only original 16-cell France 200km BOUNDING RECTANGLE, not true equal-area within-cell event incidence; no GRIIS taxon grouping is done in this source metadata/QC step.",
            "Original 48+33 S3 object-union covers only identified A5 candidate IDs, and A5 polygon function disagreed with some centroid partition IDs; complete regional event population is not proven.",
            "A5 source candidates include other geography and taxa, and source snapshot is a DEVELOPMENT preview, not a fresh GBIF authenticated download or issued DOI.",
            "S3 source rows may be indexed original interpreted events from non-independent observers; counts are not field sampling effort, ecological absence or individual plants.",
            "2026 rail network does not identify actual historical freight operations or seed dispersal.",
        ]
    }


def run(previous,geography):
    import boto3,duckdb
    from botocore import UNSIGNED
    from botocore.config import Config
    s3=boto3.client("s3",region_name="us-east-1",config=Config(signature_version=UNSIGNED))
    db=duckdb.connect()
    db.execute("SET memory_limit='3GB'; SET threads=2; SET enable_progress_bar=false")
    db.execute("INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1'")
    db.execute("INSTALL spatial; LOAD spatial; INSTALL a5 FROM community; LOAD a5")
    try:
        grid,meta=polygon_cells(db,BLOCK)
        plants,_,truncated=listing(s3,ROOT)
        if truncated:raise ValueError("Public GBIF source partition taxon listing truncated")
        taxon_groups=[v for v in plants if v[len(ROOT):].startswith("taxonpartition=Plantae")]
        if not taxon_groups:raise ValueError("No true public source Plantae taxon partitions")
        allfiles=[]
        for group in taxon_groups:
            for id in grid:
                files,_=list_all_source_parquets(s3,group+"a5_r2="+str(id)+"/")
                allfiles+=files
        result=scan_unread(db,previous,geography,allfiles)
        result["source_A5_conservative_partition_IDs"]=grid
        print(json.dumps({"complete":result["complete"],
            "n_new_public_source_objects":result["n_new_additional_source_Parquet_objects_successful"],
            "n_new_source_files_with_positive_rows":result["n_new_additional_source_Parquet_objects_with_at_least_one_source_event"],
            "n_old_original_all_plants":result["2018to2025_raw_all_Plantae_original_rows_from_old_48_files"],
            "n_new_original_all_plants":result["2018to2025_raw_all_Plantae_original_rows_from_new_33_files"],
            "fraction_of_combined_raw_source_already_in_old_48":result["fraction_of_combined_bounded_81_file_rows_covered_by_old_48"],
            "errors":len(result["source_query_errors"])
        },ensure_ascii=False),flush=True)
        return result
    finally:
        db.close()


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--old-source-v65",required=True)
    p.add_argument("--full-metadata-v67",required=True)
    p.add_argument("--output",default="results/audit/GBIF_original_FR_added_33_source_file_counts_v68.json")
    a=p.parse_args()
    old=json.loads(Path(a.old_source_v65).read_text(encoding="utf-8"))
    meta=json.loads(Path(a.full_metadata_v67).read_text(encoding="utf-8"))
    r=run(old,meta)
    file=Path(a.output);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding="utf-8")
    if not r["complete"]:
        raise SystemExit("Some original public GBIF S3 source files could not be counted; no full-source inference")


if __name__=="__main__":
    main()
