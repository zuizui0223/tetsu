#!/usr/bin/env python3
"""v356: list EVERY official Overture 2026-11-18.0 transportation segment Parquet
source object anonymously and inspect actual RAIL connector/name/OSM-source schema.

First-stage ONLY, no biological warmth/similarity result and no global
extraction claim before ALL source objects are scanned. Source S3 object
listing must be complete, with continuation tokens fully consumed.
"""
from __future__ import annotations
import argparse,hashlib,importlib.metadata,json,re,time
from collections import Counter
from pathlib import Path

RELEASE="2026-11-18.0"
REGION="us-west-2"
BUCKET="overturemaps-us-west-2"
PREFIX=f"release/{RELEASE}/theme=transportation/type=segment/"
PREDECLARED_BBOX=(14.420,50.060,14.470,50.120) # Prague, source schema-only first-party pilot
EXPECTED_PRIMARY_SCHEMA="v2.0.0"
MAX_KEYS=30000

def list_original_parquet_objects():
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    config=Config(region_name=REGION,signature_version=UNSIGNED,
        connect_timeout=20,read_timeout=65,retries={"max_attempts":3})
    s3=boto3.client("s3",region_name=REGION,config=config)
    keys=[];marker=None;pages=0;allobjects=0;source_sums=0
    while True:
        kwargs={"Bucket":BUCKET,"Prefix":PREFIX,"MaxKeys":1000}
        if marker:kwargs["ContinuationToken"]=marker
        response=s3.list_objects_v2(**kwargs)
        pages+=1
        if pages>MAX_KEYS:raise ValueError("Original S3 listing repeated pagination: unexpected >30k pages")
        for obj in response.get("Contents",[]):
            allobjects+=1
            key=obj["Key"]
            if not key.endswith(".parquet"):continue
            size=int(obj["Size"])
            if size<=0:raise ValueError("Original Overture source object size <=0")
            source_sums+=size
            keys.append({"key":key,"bytes":size,"s3_etag":obj.get("ETag","").strip('"')})
            if len(keys)>MAX_KEYS:raise ValueError("Full source registry exceeds declared 30k entries; stop")
        if not response.get("IsTruncated",False):break
        marker=response.get("NextContinuationToken")
        if not marker:raise ValueError("Original ListObjects continuation token missing; census incomplete")
    keys.sort(key=lambda x:x["key"])
    if not keys or len({x["key"] for x in keys})!=len(keys):
        raise ValueError("S3 original complete source listing empty/duplicated")
    checked=f"{BUCKET}\n{PREFIX}\n"+''.join(x["key"]+"\t"+str(x["bytes"])+"\t"+x["s3_etag"]+"\n" for x in keys)
    hash=hashlib.sha256(checked.encode()).hexdigest()
    return {"original_s3_2026_transportation_segment_directory":f"s3://{BUCKET}/{PREFIX}",
            "source_listing_all_pages_consumed":True,
            "original_listing_n_pages":pages,
            "source_total_objects_in_directory":allobjects,
            "source_parquet_object_count":len(keys),
            "source_parquet_compressed_bytes":source_sums,
            "full_object_inventory_sha256":hash,
            "all_source_parquet_objects":keys,
            "all_2026_rail_segments_NOT_EXTRACTED_BY_INVENTORY_ALONE":True}


def inspect_real_prague_rail():
    import pyarrow.compute as pc
    from overturemaps import record_batch_reader
    out={"source_probe_bbox_WGS84":PREDECLARED_BBOX,
         "source_Overture_release_requested":RELEASE}
    for attempt in (1,2):
        try:
            reader=record_batch_reader("segment",bbox=PREDECLARED_BBOX,
                                       release=RELEASE,stac=True)
            if reader is None:raise ValueError("Source rail schema sampling no original reader")
            n=Counter();names=Counter();connectors=Counter();source_osm=Counter()
            example_fields=None
            for batch in reader:
                n["source_original_transport_rows_scanned"]+=batch.num_rows
                if n["source_original_transport_rows_scanned"]>350000:
                    raise ValueError("Unexpected >350k rail preview bbox source lines")
                if example_fields is None:example_fields=batch.schema.names
                subtype_idx=batch.schema.get_field_index("subtype")
                if subtype_idx<0:raise ValueError("Source subtype required but missing")
                rail=batch.filter(pc.equal(batch.column(subtype_idx),"rail"))
                n["original_rail_source_segments"]+=rail.num_rows
                for row in rail.to_pylist():
                    if row.get("id") is None:raise ValueError("Source GERS id missing")
                    n["rail_named"]+=bool(row.get("names") and row["names"].get("primary"))
                    n["rail_route_rules_present"]+=bool(row.get("routes"))
                    n["rail_connectors_nonempty"]+=bool(row.get("connectors"))
                    for conn in row.get("connectors") or []:
                        if conn.get("connector_id") is not None:connectors[str(conn["connector_id"])]+=1
                    classes=str(row.get("class") or "unknown")
                    names[classes]+=1
                    for src in row.get("sources") or []:
                        if str(src.get("dataset")).lower()=="openstreetmap":
                            text=str(src.get("record_id") or "")
                            source_osm["OSM_way_like"]+=bool(re.fullmatch(r"w\d+@\d+",text))
                            source_osm["OSM_relation_like"]+=bool(re.fullmatch(r"r\d+@\d+",text))
                            source_osm["OSM_missing_record_id"]+=(not bool(text))
                            source_osm["OSM_sources_total"]+=1
            if n["original_rail_source_segments"]<20:
                raise ValueError("No meaningful genuine rail pilot source schema")
            return {**out,"source_success":True,
                "source_original_rail_class_counts":dict(names),
                "source_rail_named_segments":n["rail_named"],
                "source_rail_route_rules_present":n["rail_route_rules_present"],
                "source_rail_segments_with_connectors":n["rail_connectors_nonempty"],
                "source_unique_connector_IDs":len(connectors),
                "source_physical_connector_shared_by_more_than_one_segment":sum(x>1 for x in connectors.values()),
                "source_openstreetmap_record_source_types":dict(source_osm),
                "original_schema_field_names":example_fields,
                "observed_source_connectors_have_authentic_IDs":True,
                "not_equivalent_to_named_route_membership":True,
                "no_original_rail_geometries_names_or_IDs_in_output":True}
        except Exception as err:
            out[f"source_probe_attempt_{attempt}_error"]={"type":type(err).__name__,"message":str(err)[:190]}
            if attempt<2:time.sleep(8)
    return {**out,"source_success":False,"do_not_infer_no_rail_because_probe_failed":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True)
    args=p.parse_args()
    record={"version":"v356-all-official-Overture-rail-source-inventory",
            "requested_source_release":RELEASE,
            "source_full_railway_count_is_NOT_extracted_until_every_source_file_is_scanned":True,
            "source_object_inventory":None,
            "original_small_bbox_actual_rail_schema":None}
    try:
        inv=list_original_parquet_objects()
        record["source_object_inventory"]=inv
    except Exception as err:
        record["source_listing_unavailable"]={"error_type":type(err).__name__,"short_reason":str(err)[:230]}
    record["original_small_bbox_actual_rail_schema"]=inspect_real_prague_rail()
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding="utf8")
    short={k:v for k,v in record.items() if k!="source_object_inventory"}
    if record.get("source_object_inventory"):
        short["complete_original_object_inventory_aggregate"]={k:v for k,v in record[
           "source_object_inventory"].items() if k!="all_source_parquet_objects"}
    print(json.dumps(short,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
