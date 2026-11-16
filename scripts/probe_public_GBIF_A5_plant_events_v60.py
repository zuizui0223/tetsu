#!/usr/bin/env python3
"""v6.0: REAL original GBIF event rows from public A5 geographic parquet PREVIEW.

GBIF's public *experimental* August-2026 DWCA-derived A5 parquet
snapshot is an independently documented derivative of official
10.15468/dl.8yrbe7. No credentials are required. This pilot reads
only preselected original 200km 50km-grid US Midwest and FR Burgundy
geographic sources, from a SMALL, explicitly non-exhaustive subset of
public Plantae taxonomic × A5-r2 geographic partition files.

IMPORTANT: A single selected Parquet file is NOT the entire GBIF
original archive, or a source-complete country flora sample. Check
actual interpreted row status PRESENT, coordinate uncertainty,
event identifiers, taxon identifiers; output only counts and source
file metadata, NEVER GBIF individual record IDs, human identities,
or precise occurrence lat/lon.

Preflight source risk: gbif/occurrence #609 reports ~5% missing files
in regular monthly snapshots 2025-07 onward; we are using the separate
AUGUST 2026 DWCA development preview only and must separately audit
completeness before ecological analysis. DOI underlying published
original snapshot is NOT newly issued by tetsu.

If geographic pruning or remote query fails, fail closed; no botanical
effect or seed dispersal is inferred.
"""
from __future__ import annotations
import argparse,json,math
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import quote

from shapely import wkt
from shapely.ops import unary_union
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS,make_grid

BUCKET="gbif-public-data"
ROOT="development/2026-08-01-DWCA/occurrence/p_taxon_a5/"
UNDERLYING_GBIF_DOI="10.15468/dl.8yrbe7"
FROZEN_BLOCKS=("us_midwest","fr_burgundy")
MAX_TAXON_PREFIX_GROUPS=60
MAX_SOURCE_FILES_TO_TRY=5
MAX_ONE_PARQUET_OBJECT_BYTES=400_000_000
MAX_REGION_EVENT_ROW_COUNT_PER_SELECTED_FILE=1000000


def equal_area_source_cells(block):
    if block not in FROZEN_BLOCKS:raise ValueError("Nonfrozen original 50km source block")
    rows=make_grid(block)
    if len(rows)!=16 or any(r["country"]!=BLOCKS[block][0] for r in rows):
        raise ValueError("Original 50km source area changed")
    return rows


def source_a5_cells(db,block):
    """Call the authoritative DuckDB A5 community extension on 16 source centroids.

    This is a limited source pilot: a5_r2 centroid covers for the
    original source region, NOT an exhaustive all-overlapping-cells
    listing needed for full data completeness.
    """
    cells=equal_area_source_cells(block)
    values=[]
    for cell in cells:
        p=wkt.loads(cell["gbif_geometry_wgs84"]).centroid
        cid=db.execute("SELECT a5_lonlat_to_cell(?,?,2)",[p.x,p.y]).fetchone()[0]
        if isinstance(cid,bool) or not isinstance(cid,int) or cid<=0:
            raise ValueError("Actual A5 source r2 encoded ID invalid")
        values.append(int(cid))
    return sorted(set(values))


def listing(client,prefix,delimiter="/"):
    page=client.list_objects_v2(Bucket=BUCKET,Prefix=prefix,
                                Delimiter=delimiter,MaxKeys=250)
    groups=[a.get("Prefix") for a in page.get("CommonPrefixes",[])
            if isinstance(a,dict) and a.get("Prefix")]
    files=[{"key":r["Key"],"size_bytes":int(r["Size"])}
           for r in page.get("Contents",[])
           if isinstance(r,dict) and isinstance(r.get("Key"),str)
           and isinstance(r.get("Size"),int) and r["Key"].endswith(".parquet")]
    return groups,files,bool(page.get("IsTruncated"))


def find_small_geographic_plant_files(s3,source_a5):
    roots,_,truncated=listing(s3,ROOT)
    plants=sorted(x for x in roots if x[len(ROOT):].lower().startswith("taxonpartition=plantae"))
    if not plants:
        raise ValueError("No anonymous public GBIF Plantae partition prefix")
    if truncated:
        raise ValueError("GBIF source taxon partition list truncated; cannot safely select")
    candidates=[]
    geos=sorted(set(source_a5))
    for parent in plants:
        for cell in geos:
            prefix=parent+"a5_r2="+str(cell)+"/"
            _,objects,partial=listing(s3,prefix)
            if partial:raise ValueError("A5 source object listing truncated")
            for obj in objects:
                if obj["size_bytes"]<=MAX_ONE_PARQUET_OBJECT_BYTES and obj["size_bytes"]>100:
                    candidates.append({**obj,"a5_r2":cell,
                                       "taxon_partition":parent[len(ROOT):].rstrip("/")})
    candidates.sort(key=lambda r:(r["key"],r["size_bytes"]))
    return plants,candidates


def read_true_public_parquet_source(db,block,candidate):
    key=candidate["key"]
    if not key.startswith(ROOT) or ".." in key:
        raise ValueError("Cannot query nonapproved external GBIF preview object")
    country=BLOCKS[block][0]
    file=f"s3://{BUCKET}/{key}"
    # No SQL injection: key from the approved public S3 bucket listing,
    # exact approved prefix. Single-quote source names are rejected.
    if "'" in file:raise ValueError("Unsafe GBIF source key")
    source=f"read_parquet('{file}',hive_partitioning=true,hive_types={{'a5_r2':UBIGINT}})"
    schema=db.execute(f"DESCRIBE SELECT * FROM {source}").fetchall()
    types={str(row[0]).lower():str(row[0]) for row in schema}
    mandatory=("gbifid","datasetkey","specieskey","countrycode","year",
               "decimallatitude","decimallongitude","occurrencestatus",
               "coordinateuncertaintyinmeters")
    if any(x not in types for x in mandatory):
        raise ValueError("Public original GBIF Parquet lacks expected interpreted event schema")
    def q(n):
        return '"'+types[n].replace('"','""')+'"'
    rows=equal_area_source_cells(block)
    footprint=unary_union([wkt.loads(c["gbif_geometry_wgs84"]) for c in rows])
    west,south,east,north=footprint.bounds
    bbox=(
        f"{q('decimallongitude')} BETWEEN {west:.9f} AND {east:.9f} "
        f"AND {q('decimallatitude')} BETWEEN {south:.9f} AND {north:.9f}"
    )
    where=(f"{q('countrycode')}='{country}' AND {q('year')} BETWEEN 2018 AND 2025 "
           f"AND {bbox}")
    # Count sample only; exact original 50km projected cell membership
    # and original country GRIIS status would be analyzed downstream,
    # not inferred from this smaller country bbox source sample.
    query=f"""SELECT
       count(*) AS source_rows_within_bounding_envelope,
       count(*) FILTER (WHERE UPPER({q('occurrencestatus')})='PRESENT') AS indexed_present_rows,
       count(*) FILTER (WHERE UPPER({q('occurrencestatus')})='ABSENT') AS indexed_absent_rows,
       count(*) FILTER (WHERE {q('coordinateuncertaintyinmeters')} IS NULL) AS missing_coordinate_uncertainty,
       count(DISTINCT {q('datasetkey')}) AS n_source_datasets,
       count(DISTINCT {q('specieskey')}) AS n_source_specieskeys
       FROM {source}
       WHERE {where}"""
    found=db.execute(query).fetchone()
    if not found or len(found)!=6 or found[0]<0 or found[0]>MAX_REGION_EVENT_ROW_COUNT_PER_SELECTED_FILE:
        raise ValueError("Remote source returned invalid/unexpectedly high pilot record count")
    return {
        "selected_public_parquet_key":key,
        "selected_public_parquet_size_bytes":candidate["size_bytes"],
        "selected_a5_r2":candidate["a5_r2"],
        "selected_source_taxon_partition":candidate["taxon_partition"],
        "source_query_country":country,
        "source_query_earliest_year":2018,
        "source_query_latest_year":2025,
        "source_query_only_WGS84_bounding_envelope_NOT_original_16_cell_exact_membership":True,
        "n_actual_original_interpreted_rows_in_pilot_bbox":found[0],
        "n_Present_status_source_rows":found[1],
        "n_Absent_status_source_rows":found[2],
        "n_missing_coordinate_uncertainty_source_rows":found[3],
        "n_distinct_GBIF_dataset_UUIDs_in_sample":found[4],
        "n_distinct_species_keys_in_sample":found[5],
        "real_public_original_Parquet_source_rows_READ":True,
        "the_sample_is_NOT_full_GBIF_source_or_GRIIS_intro_only":True,
    }


def run(blocks):
    import boto3,duckdb
    from botocore import UNSIGNED
    from botocore.config import Config
    if not blocks or len(blocks)!=len(set(blocks)) or any(b not in FROZEN_BLOCKS for b in blocks):
        raise ValueError("Source pilot only in predefined US/FR 50km macro blocks")
    s3=boto3.client("s3",region_name="us-east-1",
                   config=Config(signature_version=UNSIGNED,retries={"max_attempts":2}))
    db=duckdb.connect()
    db.execute("SET threads=2")
    db.execute("SET memory_limit='3GB'")
    db.execute("SET enable_progress_bar=false")
    db.execute("INSTALL httpfs; LOAD httpfs;")
    db.execute("SET s3_region='us-east-1';")
    db.execute("INSTALL a5 FROM community; LOAD a5;")
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_public_preview_root":f"s3://{BUCKET}/{ROOT}",
        "underlying_GBIF_August2026_original_download_DOI":UNDERLYING_GBIF_DOI,
        "not_a_NEW_GBIF_DWCA_download_or_DOI":True,
        "preview_is_experimental_and_NOT_completeness_validated":True,
        "regions":{},"errors":{},
    }
    for block in blocks:
        try:
            ids=source_a5_cells(db,block)
            plants,candidates=find_small_geographic_plant_files(s3,ids)
            if not candidates:
                raise ValueError("No bounded public Plantae source Parquet objects at original 50km centroids")
            sampled=[]
            failures=[]
            for candidate in candidates[:MAX_SOURCE_FILES_TO_TRY]:
                try:
                    row=read_true_public_parquet_source(db,block,candidate)
                    sampled.append(row)
                    print(json.dumps({"block":block,"public_original_event_source":row},
                                    ensure_ascii=False),flush=True)
                    if row["n_actual_original_interpreted_rows_in_pilot_bbox"]>0:break
                except Exception as exc:
                    failures.append(f"{type(exc).__name__}: {str(exc)[:180]}")
            if not sampled or all(x["n_actual_original_interpreted_rows_in_pilot_bbox"]==0 for x in sampled):
                raise ValueError("No actual original GBIF source rows could be verified in bounded plant preview; "+str(failures[:2]))
            out["regions"][block]={
                "country":BLOCKS[block][0],
                "n_real_a5_r2_source_centroid_partitions":len(ids),
                "source_plant_taxon_partition_prefixes":[p[len(ROOT):] for p in plants],
                "n_actual_candidate_GBIF_public_Plant_parquet_objects_in_source_a5_centroids":len(candidates),
                "n_pilot_public_parquet_objects_scanned":len(sampled),
                "source_query_failures_before_success":failures,
                "pilot_interpreted_occurrence_files":sampled,
                "not_full_geographically_partitioned_source":True,
            }
        except Exception as exc:
            out["errors"][block]=f"{type(exc).__name__}: {str(exc)[:250]}"
            print(json.dumps({"source_error":block,"error":out["errors"][block]},
                             ensure_ascii=False),flush=True)
    db.close()
    out["complete"]=not out["errors"] and len(out["regions"])==len(blocks)
    out["limits"]=[
        "Only one or few bounded actual Parquet objects selected by lexical source key, not all plant taxon A5 partitions or 50km polygons.",
        "A5 resolution 2 partitions were selected by original grid centroids, not full source region boundary overlap; border partitions can be missing.",
        "Source query uses country+period+BBOX, not exact equal-area cell membership; do not claim original regional study counts or ecological inference.",
        "Preview is experimental August 2026 DOI-derived original GBIF data; full original snapshot extraction/completeness not verified.",
        "GBIF record IDs, human names and actual individual source coordinates were not exported.",
        "No 2026 physical railway freight or plant dispersal, observer-quality correction or country GRIIS taxa analysis in this source-pilot.",
    ]
    print(json.dumps({"complete":out["complete"],
                      "n_blocks_with_true_public_original_occurrence_rows":len(out["regions"]),
                      "errors":out["errors"]},ensure_ascii=False),flush=True)
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--blocks",nargs="+",choices=FROZEN_BLOCKS,default=["us_midwest"])
    p.add_argument("--output",default="results/audit/GBIF_public_A5_original_Plant_event_source_v60.json")
    a=p.parse_args()
    result=run(a.blocks)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    if not result["complete"]:
        raise SystemExit("No real GBIF public original source observed in one-file A5 geographic pilot")


if __name__=="__main__":
    main()
