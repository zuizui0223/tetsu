#!/usr/bin/env python3
"""v5.9: Probe actual ANONYMOUS GBIF original-occurrence Parquet partitions.

A public, historical GBIF-derived DWCA snapshot from 2026-08-01 is
documented by GBIF on 2026-10-19:
 s3://gbif-public-data/development/2026-08-01-DWCA/
 occurrence/p_taxon_a5/ and occurrence/p_taxon_h3/
 DOI of UNDERLYING published August full download:
 https://doi.org/10.15468/dl.8yrbe7

This is NOT a claim the development preview is itself a separate DOI
download, 100% complete, identical to the November indexed GBIF
facets, or a production-grade archive. ONLY public S3 object metadata
is read and reported. No credentials, observer IDs, exact coordinates,
occurrence events or private identifiers are stored.

Additional safety gate: check a known missing object in the older
monthly GBIF Open Data October 2026 snapshot (gbif/occurrence #609).
Anonymous S3 HEAD 000019 should work while 000020 is reported missing,
despite 2026-10 snapshot citation full-download count. If verified,
DO NOT treat the ordinary snapshot as a complete population census.

This source preflight does not execute ecological statistics or download
a 3-billion-row data corpus.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

SOURCE_BUCKET="gbif-public-data"
SOURCE_DATE="2026-08-01"
PREVIEW_PREFIX="development/2026-08-01-DWCA"
PREVIEW_PATHS=(
    "occurrence/p_taxon_a5/",
    "occurrence/p_taxon_h3/",
)
SOURCE_PUBLISHED_DWCA_DOI="10.15468/dl.8yrbe7"
AWS_PUBLIC_SNAPSHOT_BUCKET="gbif-open-data-us-east-1"
SUSPECT_MONTH="2026-10-01"
LEGACY_CHUNKS=("000019","000020")
MAX_PREFIX_SOURCES=100
MAX_INSPECTED_FILES=3


def _is_missing_s3_error(exc):
    response=getattr(exc,"response",{})
    code=str((response.get("Error") or {}).get("Code",""))
    return code in ("404","NoSuchKey","NotFound")


def list_public_source_prefix(client,bucket,prefix):
    """One bounded, nonrecursive list of source taxon/geo partitions."""
    result=client.list_objects_v2(
        Bucket=bucket,Prefix=prefix,Delimiter="/",MaxKeys=MAX_PREFIX_SOURCES
    )
    folders=[x.get("Prefix") for x in result.get("CommonPrefixes",[])
             if isinstance(x,dict) and isinstance(x.get("Prefix"),str)]
    direct=[{"filename":Path(x["Key"]).name,
             "size_bytes":x["Size"],
             "last_modified_utc":x["LastModified"].isoformat() if hasattr(x.get("LastModified"),"isoformat") else None
            } for x in result.get("Contents",[])
            if isinstance(x,dict) and x.get("Key") and isinstance(x.get("Size"),int)]
    if not folders and not direct:
        raise ValueError("No actual GBIF public preview objects in selected prefix")
    return {
        "source_s3_bucket":bucket,
        "source_prefix":prefix,
        "n_partition_prefixes_in_bounded_listing":len(folders),
        "n_direct_objects_in_bounded_listing":len(direct),
        "preview_prefixes":folders[:MAX_PREFIX_SOURCES],
        "direct_object_metadata":direct[:MAX_INSPECTED_FILES],
        "S3_list_is_truncated":bool(result.get("IsTruncated")),
        "does_NOT_certify_whole_snapshot_completeness":True,
    }


def validate_object_sample(client,bucket,key,expected_present):
    try:
        h=client.head_object(Bucket=bucket,Key=key)
        return {
            "key":key,
            "present":True,
            "size_bytes":int(h.get("ContentLength",0)),
            "reported_content_type":str(h.get("ContentType") or ""),
            "matches_expected_presence":expected_present,
        }
    except Exception as exc:
        if _is_missing_s3_error(exc):
            return {"key":key,"present":False,
                    "matches_expected_presence":not expected_present}
        raise


def run(client_preview,client_legacy):
    audit={
        "utc":datetime.now(timezone.utc).isoformat(),
        "preview_bucket":SOURCE_BUCKET,
        "preview_date":SOURCE_DATE,
        "preview_is_EXPERIMENTAL":True,
        "preview_underlying_source_GBIF_DWCA_DOI":SOURCE_PUBLISHED_DWCA_DOI,
        "GBIF_preview_documentation":"https://data-blog.gbif.org/post/extending-parquet-snapshots/",
        "previous_monthly_snapshot_reported_incomplete":"https://github.com/gbif/occurrence/issues/609",
        "public_preview_source_prefixes":{},
        "legacy_suspect_snapshot_examples":{},
        "errors":{},
    }
    for path in PREVIEW_PATHS:
        prefix=PREVIEW_PREFIX+"/"+path
        try:
            audit["public_preview_source_prefixes"][path]=list_public_source_prefix(
                client_preview,SOURCE_BUCKET,prefix)
            print(json.dumps({"source":path,
                "prefix_groups":audit["public_preview_source_prefixes"][path]["n_partition_prefixes_in_bounded_listing"],
                "direct_files":audit["public_preview_source_prefixes"][path]["n_direct_objects_in_bounded_listing"],
                "first_groups":audit["public_preview_source_prefixes"][path]["preview_prefixes"][:8]
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            audit["errors"]["preview_"+path]=f"{type(exc).__name__}: {str(exc)[:180]}"
    for filename in LEGACY_CHUNKS:
        key=f"occurrence/{SUSPECT_MONTH}/occurrence.parquet/{filename}"
        try:
            audit["legacy_suspect_snapshot_examples"][filename]=validate_object_sample(
                client_legacy,AWS_PUBLIC_SNAPSHOT_BUCKET,key,expected_present=filename=="000019")
        except Exception as exc:
            audit["errors"]["legacy_"+filename]=f"{type(exc).__name__}: {str(exc)[:180]}"
    audit["preview_list_source_accessible"]=len(audit["public_preview_source_prefixes"])==len(PREVIEW_PATHS)
    audit["legacy_reproduces_published_missing_chunk_issue"]=(
        audit["legacy_suspect_snapshot_examples"].get("000019",{}).get("present") is True
        and audit["legacy_suspect_snapshot_examples"].get("000020",{}).get("present") is False
    )
    audit["complete_probe"]=not audit["errors"] and audit["preview_list_source_accessible"]
    audit["no_original_occurrence_rows_claimed_or_downloaded"]=True
    audit["no_complete_GBIF_event_population_claimed"]=True
    audit["interpretation_limits"]=[
        "Anonymous S3 listing/HEAD certifies only specific bucket prefixes and objects, not data content, source completeness, taxonomic interpretation or country-cell event quality.",
        "The preview is a DEVELOPMENT layout using August 2026 GBIF source, with its underlying original download DOI; it is not a fresh DOI created by this audit.",
        "Old public monthly snapshot has an independently reported approximately 4.9% missing-object problem from July 2025 onward; never infer full population from an incomplete directory listing.",
        "GBIF's October 19 2026 blog provides DuckDB A5/H3 geographic partition queries; further schema and locality tests are necessary before remote botanical source sampling.",
        "No user GBIF credentials were read, no new download request submitted, and no railway seed-propagule result was calculated."
    ]
    print(json.dumps({
        "preview_list_source_accessible":audit["preview_list_source_accessible"],
        "legacy_oct2026_missing_20th_chunk_reproduced":audit["legacy_reproduces_published_missing_chunk_issue"],
        "complete_probe":audit["complete_probe"],
        "errors":audit["errors"]
    },ensure_ascii=False),flush=True)
    return audit


def main():
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/GBIF_public_A5_H3_source_preflight_v59.json")
    args=p.parse_args()
    client=boto3.client("s3",region_name="us-east-1",
                        config=Config(signature_version=UNSIGNED,
                                      retries={"max_attempts":3,"mode":"standard"}))
    audit=run(client,client)
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding="utf-8")
    if not audit["complete_probe"]:
        raise SystemExit("Public GBIF Parquet source not yet accessible/verified; no inference")


if __name__=="__main__":
    main()
