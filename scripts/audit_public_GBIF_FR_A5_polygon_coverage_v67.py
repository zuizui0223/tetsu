#!/usr/bin/env python3
"""v6.7 GBIF public A5 complete source-POLYGON spatial partition metadata audit.

v6.0-v6.6 chose A5 r2 partitions containing original 50km cell
CENTROIDS, which is NOT guaranteed to cover cells crossing A5 grid
boundaries. Do not interpret original plant zero observations until the
official A5_geometry_to_cells(polygon,2) coverage is compared.

Retrieve FULL A5 partition coverage using the published GBIF DuckDB A5
community extension on all 16 frozen 50km source polygons. Then list
ALL public GBIF experimental August2026 Plantae Parquet object metadata
for those A5 partitions through anonymous S3 (no occurrence records,
no credentials, no GBIF DOI issued). Compare old v6.5 48 scanned file
objects with this complete spatial candidate list, count unscanned or
oversized Parquet object coverage and true resource budget.

This does not assert that even the experimental snapshot itself is
globally complete, identical to current GBIF indexed data or biologically
unbiased. No rail/seed ecology or species outcomes are used.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from audit_gbif_equalarea_macroblocks_v29 import BLOCKS,make_grid
from probe_public_GBIF_A5_plant_events_v60 import (
    BUCKET,ROOT,UNDERLYING_GBIF_DOI,
    source_a5_cells,listing,MAX_ONE_PARQUET_OBJECT_BYTES
)

BLOCK="fr_burgundy"
MAX_A5_SOURCE_IDS=20
MAX_SOURCE_FILE_METADATA=2500
MAX_S3_LIST_PAGE=250


def polygon_cells(db,block=BLOCK):
    if block!=BLOCK:raise ValueError("Original source geography not frozen France Burgundy")
    all_cells=set()
    audit={}
    rows=make_grid(block)
    if len(rows)!=16:
        raise ValueError("Expected original 4x4 equal-area geography")
    for r in rows:
        polygon=r["gbif_geometry_wgs84"]
        got=db.execute(
            "SELECT unnest(a5_geometry_to_cells(ST_GeomFromText(?), 2))",
            [polygon]).fetchall()
        ids={int(v[0]) for v in got}
        if not ids or any(k<=0 for k in ids):
            raise ValueError("Official A5 polygon coverage not found")
        audit[r["cell_id"]]=sorted(ids)
        all_cells|=ids
    if len(all_cells)>MAX_A5_SOURCE_IDS:
        raise ValueError("True source polygons cover too many A5 grid partitions for bounded audit")
    center=set(source_a5_cells(db,block))
    # Official A5 polygon cell cover is NOT always consistent with
    # independently calculated A5 lon/lat centroid membership for
    # small source geometries. Do not conceal the disagreement or
    # call one an exhaustive spatial cover: combine both as a
    # conservative candidate set, and publish the discrepancy.
    polygon_only=set(all_cells)
    not_in_polygon=sorted(center-polygon_only)
    additional=sorted(polygon_only-center)
    all_cells.update(center)
    if len(all_cells)>MAX_A5_SOURCE_IDS:
        raise ValueError("Conservative A5 source ID union exceeds bounded source metadata budget")
    return sorted(all_cells),{
        "n_original_50km_source_cells":len(rows),
        "n_A5_polygon_cover_partition_ids":len(polygon_only),
        "n_A5_conservative_polygon_plus_centroid_partition_ids":len(all_cells),
        "n_A5_centroid_partition_ids":len(center),
        "extra_A5_partition_ids_missed_by_centroids":additional,
        "centroid_A5_IDs_missing_from_official_polygon_function":not_in_polygon,
        "polygon_cover_disagreed_with_centroid_lookup":bool(not_in_polygon),
        "true_complete_polygon_intersection_cover_NOT_verified_if_disagreement":bool(not_in_polygon),
        "source_A5_r2_cell_cover_by_original_50km_cells":audit,
        "original_centroid_A5_partition_ids":sorted(center),
    }


def list_all_source_parquets(client,root):
    if not root.startswith(ROOT+"taxonpartition=Plantae") or "a5_r2=" not in root or not root.endswith("/"):
        raise ValueError("Unapproved S3 GBIF A5 Plantae source metadata prefix")
    files=[]
    continuation=None
    for i in range(10):
        args={"Bucket":BUCKET,"Prefix":root,"MaxKeys":MAX_S3_LIST_PAGE}
        if continuation:args["ContinuationToken"]=continuation
        r=client.list_objects_v2(**args)
        for obj in r.get("Contents",[]) or []:
            key=obj.get("Key");size=obj.get("Size")
            if not isinstance(key,str) or not key.startswith(root) or not key.endswith(".parquet") or not isinstance(size,int) or size<=0:
                continue
            files.append({"key":key,"size_bytes":size})
        if len(files)>MAX_SOURCE_FILE_METADATA:
            raise ValueError("Full A5 file metadata source list exceeds safe cap")
        if not r.get("IsTruncated"):
            if len({x["key"] for x in files})!=len(files):
                raise ValueError("Duplicate Parquet file source objects returned across S3 listing")
            return files,i+1
        continuation=r.get("NextContinuationToken")
        if not continuation:
            raise ValueError("Truncated source S3 page with no continuation token")
    raise ValueError("GBIF A5 object listing still truncated after bounded pagination; never infer complete object population")


def analyze(old,db,client):
    if old.get("complete") is not True or set(old.get("regions",{}))!={BLOCK}:
        raise ValueError("Prior v6.5 source file scan success required")
    oldregion=old["regions"][BLOCK]
    if oldregion.get("n_source_objects_actually_scanned")!=48:
        raise ValueError("Original file selection not the frozen 48 source objects")
    oldkeys={r.get("source_key") for r in oldregion["source_object_results"]}
    if len(oldkeys)!=48 or any(not isinstance(v,str) for v in oldkeys):
        raise ValueError("Missing actual v6.5 source selected file keys")
    a5ids,geog=polygon_cells(db,BLOCK)
    roots,_,truncated=listing(client,ROOT)
    if truncated:
        raise ValueError("Global GBIF Plantae source partition prefix listing truncated")
    plants=sorted(v for v in roots if v[len(ROOT):].startswith("taxonpartition=Plantae"))
    if not plants:
        raise ValueError("No official experimental Plantae taxonomic source partitions")
    group_results={}
    all_files=[]
    for parent in plants:
        for a5id in a5ids:
            prefix=parent+"a5_r2="+str(a5id)+"/"
            files,npages=list_all_source_parquets(client,prefix)
            all_files.extend(files)
            group_results[prefix]={
                "n_source_Parquet_objects":len(files),
                "n_S3_metadata_listing_pages":npages,
                "n_objects_larger_than_old_400MB_read_limit":sum(
                    f["size_bytes"]>MAX_ONE_PARQUET_OBJECT_BYTES for f in files),
                "total_source_compressed_Parquet_bytes":sum(f["size_bytes"] for f in files),
            }
    keys={f["key"] for f in all_files}
    if len(keys)!=len(all_files):
        raise ValueError("Duplicate GBIF Parquet objects across A5 partitions")
    if not oldkeys.issubset(keys):
        raise ValueError("Original 48 v6.5 source objects missing from full polygon A5 source list")
    output={
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_region":BLOCK,
        "country":"FR",
        "original_country_area_source":"16 equal-area 2500km2 EPSG6933 cells",
        "underlying_source_GBIF_DWCA_DOI":UNDERLYING_GBIF_DOI,
        "actual_source_original_event_rows_NOT_downloaded_in_this_metadata_audit":True,
        "source_50km_polygon_plus_centroid_A5_grid_audit":geog,
        "source_public_Plant_parquet_prefixes":[v[len(ROOT):] for v in plants],
        "n_total_source_parquet_files_in_all_true_A5_polygon_intersecting_partitions":len(all_files),
        "n_original_lexically_fixed_48_Parquet_source_files":len(oldkeys),
        "n_additional_source_files_not_covered_by_old_v65_query":len(keys-oldkeys),
        "n_large_files_omitted_by_old_400MB_filter":sum(
            f["size_bytes"]>MAX_ONE_PARQUET_OBJECT_BYTES for f in all_files),
        "source_compressed_all_candidate_size_bytes":sum(x["size_bytes"] for x in all_files),
        "fraction_files_scanned_by_v65":round(len(oldkeys)/len(keys),6) if keys else None,
        "source_file_metadata_by_A5_partition":group_results,
        "complete_source_object_metadata_listing_WITHIN_THE_COMPUTED_UNION_A5_IDs":True,
        "full_original_50km_polygon_spatial_A5_coverage_guaranteed":not geog["polygon_cover_disagreed_with_centroid_lookup"],
        "not_a_full_GBIF_original_event_population":True,
        "not_causal_train_seed_dispersal":True,
        "limits":[
            "This is a metadata-only full Parquet object listing within the conservative UNION of official polygon-cover and centroid-index A5 IDs; the two functions disagreed on some IDs, so full polygon spatial coverage is NOT independently guaranteed.",
            "GBIF preview is developmental derivative of August2026 official DWCA DOI, NOT a new DOI or user-authenticated bulk download.",
            "Number of source objects or compressed bytes does NOT equal number of original country GRIIS plant records.",
            "Total source objects include many taxonomic/other geography records in A5; row-level country/year/true 50km-cell spatial filters remain essential.",
            "Old v6.5 source selection used only A5 cell CENTROIDS and <=400MB file size, max48 object budget; a missed source object does not imply missing botanical taxa or absence.",
            "No original GBIF human names or observation locations, no railway freight line operations, no seed propagules or biological effect inferred.",
        ]
    }
    print(json.dumps({
        "block":BLOCK,
        "n_original_50km_cells":geog["n_original_50km_source_cells"],
        "n_official_polygon_A5_source_partition_ids":geog["n_A5_polygon_cover_partition_ids"],
        "n_conservative_union_A5_source_partition_ids":geog["n_A5_conservative_polygon_plus_centroid_partition_ids"],
        "n_centroid_A5_IDs_not_in_polygon_cover":len(geog["centroid_A5_IDs_missing_from_official_polygon_function"]),
        "n_centroid_A5_source_partition_ids":geog["n_A5_centroid_partition_ids"],
        "n_centroid_missed_A5_source_partition_ids":len(geog["extra_A5_partition_ids_missed_by_centroids"]),
        "n_all_source_parquet_objects":len(all_files),
        "n_source_files_v65_scanned":len(oldkeys),
        "n_missing_additional_parquet_objects":len(keys-oldkeys),
        "n_above_original_400MB_cap":output["n_large_files_omitted_by_old_400MB_filter"],
        "complete_metadata_only":True
    },ensure_ascii=False),flush=True)
    return output


def main():
    import boto3,duckdb
    from botocore import UNSIGNED
    from botocore.config import Config
    p=argparse.ArgumentParser()
    p.add_argument("--original-v65-source-json",required=True)
    p.add_argument("--output",default="results/audit/GBIF_FR_full_polygon_A5_source_metadata_v67.json")
    a=p.parse_args()
    old=json.loads(Path(a.original_v65_source_json).read_text(encoding="utf-8"))
    s3=boto3.client("s3",region_name="us-east-1",
                    config=Config(signature_version=UNSIGNED))
    db=duckdb.connect()
    db.execute("INSTALL spatial; LOAD spatial")
    db.execute("INSTALL a5 FROM community; LOAD a5")
    try:result=analyze(old,db,s3)
    finally:db.close()
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,"metadata_source_only":True}),flush=True)


if __name__=="__main__":
    main()
