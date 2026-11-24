#!/usr/bin/env python3
"""v368a: inspect AUTHENTIC worldwide iNaturalist public monthly S3 metadata.

Official source observations.csv.gz and taxa.csv.gz are MONTHLY,
photo-linked observation data, not full standardized vegetation plot
survey. Important missing schema: geoprivacy/coordinates_obscured.
Never pretend exact rail80m inference without external privacy QA.

Only HEAD source compressed sizes/ETags/times and a tiny gzip prefix.
No original photo, observer, original GPS, or taxon-by-observation rows
persisted/rebroadcast.
"""
from __future__ import annotations
import argparse,hashlib,json,zlib,datetime
from pathlib import Path

BUCKET="inaturalist-open-data"
FILES=("observations.csv.gz","taxa.csv.gz")
REQ_OBS={"observation_uuid","observer_id","latitude","longitude",
         "positional_accuracy","taxon_id","quality_grade","observed_on"}
REQ_TAX={"taxon_id","ancestry","rank","name","active"}
PRIVACY_REQUIRED_BUT_ABSENT={"geoprivacy","taxon_geoprivacy",
                              "coordinates_obscured"}
REGION="us-east-1"
RANGE_BYTES=1_048_575

def source_s3():
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    return boto3.client("s3",region_name=REGION,
           config=Config(signature_version=UNSIGNED,region_name=REGION,
                         connect_timeout=30,read_timeout=90,
                         retries={"max_attempts":3}))

def read_metadata(client,key):
    info=client.head_object(Bucket=BUCKET,Key=key)
    size=int(info["ContentLength"])
    if size<10_000:raise ValueError("Global original source metadata suspiciously small")
    with client.get_object(Bucket=BUCKET,Key=key,
         Range=f"bytes=0-{min(size-1,RANGE_BYTES)}")["Body"] as f:
        fragment=f.read()
    if not fragment.startswith(b"\x1f\x8b"):
        raise ValueError("Original AWS global iNaturalist metadata not GZIP")
    decompressor=zlib.decompressobj(16+zlib.MAX_WBITS)
    data=decompressor.decompress(fragment,3_000_000)
    lines=data.split(b"\n")
    if len(lines)<2:raise ValueError("Original global metadata header sample not parsed")
    cols=lines[0].decode("utf-8-sig").strip().split("\t")
    if len(cols)!=len(set(cols)) or len(cols)<5:
        raise ValueError("Original S3 TSV header invalid")
    expected=REQ_OBS if key.startswith("observations") else REQ_TAX
    if not expected.issubset(cols):
        raise ValueError("Original global source schema changed: "+str(expected-set(cols)))
    date=info.get("LastModified")
    return {
      "source_s3_key":f"s3://{BUCKET}/{key}",
      "source_bucket_region":REGION,
      "source_original_compressed_bytes":size,
      "original_source_ETag":info.get("ETag","").strip('"'),
      "original_source_LastModified_UTC":date.isoformat() if date else None,
      "original_source_version_id":info.get("VersionId"),
      "original_source_content_encoding":info.get("ContentEncoding"),
      "original_source_SHA256_from_s3_header":info.get("ChecksumSHA256"),
      "source_gzip_prefix_header_sha256":hashlib.sha256(lines[0]).hexdigest(),
      "original_TSV_source_schema_fields":cols,
      "source_audited_compressed_prefix_bytes":len(fragment),
      "source_first_uncompressed_line_values_not_redistributed":True}

def run():
    client=source_s3()
    data={k:read_metadata(client,k) for k in FILES}
    observation_cols=set(data["observations.csv.gz"]["original_TSV_source_schema_fields"])
    taxon_cols=set(data["taxa.csv.gz"]["original_TSV_source_schema_fields"])
    if PRIVACY_REQUIRED_BUT_ABSENT & observation_cols:
        msg="Original source includes privacy metadata, inspect source version and update frozen design"
    else:msg="Original observational metadata lacks geoprivacy/coordinates_obscured; fine-scale railway point inference is NOT authorized without independent original API privacy verification"
    return {"version":"v368a-real-monthly-global-iNaturalist-observation-and-taxon-metadata-source",
     "source_fetched_at_model_date":"2026-11-24",
     "original_rail_source_2026_11_18_already_globally_complete":True,
     "original_global_iNaturalist_open_data_monthly_source":data,
     "original_Plant_ancestry_taxon_id_candidate":47126,
     "original_public_metadata_global_taxonomy_species_rank_MUST_be_checked_before_plant_incidence":True,
     "original_source_observation_privacy_schema_gate":msg,
     "must_NOT_claim_all_iNaturalist_observations_in_photos_linked_source":True,
     "must_NOT_claim_railway_80m_ecological_plant_effect":True,
     "original_observation_geoprivacy_not_verified":True,
     "raw_sourced_photo_GPS_observer_IDs_and_taxon_incidence_not_exported":True}

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True)
    a=p.parse_args()
    obj=run()
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(obj,ensure_ascii=False),flush=True)
