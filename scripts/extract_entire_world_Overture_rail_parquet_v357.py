#!/usr/bin/env python3
"""v357: full original Overture global rail-subtype extraction from 128 source parquet files.

Pilot mode reads only a bounded original single-file source prefix and audits
Arrow public anonymous S3 parquet. Full partition mode scans EVERY source
file assigned to this shard, emits rail-only GeoParquet in artifact, and
aggregates class/flags/name/connectors/geodesic lengths. NO cherry-picking of
source region, names, route or railway species. A successfully completed shard
is NOT yet a complete global network until all 16 source shards are verified.
"""
from __future__ import annotations
import argparse,hashlib,json,math,os,re,time
from collections import Counter
from pathlib import Path

from audit_global_Overture_all_rail_parquet_inventory_v356 import (
    list_original_parquet_objects,BUCKET,PREFIX,RELEASE
)

SHARDS=16
PILOT_LIMIT_SOURCE_ROWS=65_536
BATCH=32_768
EXPECTED_SOURCE_PARQUET_OBJECTS=128
EXPECTED_SOURCE_BYTES=72_770_647_616
EXPECTED_SOURCE_INVENTORY_SHA="6237da737bb783f23f352638cd5aef287bbb3d7aa8f0fe0bed8c1cc4558f6b34"
COLUMNS=("id","subtype","class","connectors","names","sources","rail_flags","geometry")
RAIL_CLASSES=("standard_gauge","narrow_gauge","unknown","subway","tram","light_rail","funicular","monorail")
MAX_RAIL_RECORD_BYTES=1_200_000_000

def inventory():
    x=list_original_parquet_objects()
    if (x["source_parquet_object_count"]!=EXPECTED_SOURCE_PARQUET_OBJECTS
        or x["source_parquet_compressed_bytes"]!=EXPECTED_SOURCE_BYTES
        or x["full_object_inventory_sha256"]!=EXPECTED_SOURCE_INVENTORY_SHA):
        raise ValueError("Official global original Overture source file inventory changed/incomplete")
    return x


def assigned_files(source,idx):
    if idx<0 or idx>=SHARDS:raise ValueError("Invalid fixed global source shard index")
    orig=source["all_source_parquet_objects"]
    rows=[r for i,r in enumerate(orig) if i%SHARDS==idx]
    if len(rows)!=EXPECTED_SOURCE_PARQUET_OBJECTS//SHARDS:
        raise ValueError("Fixed 128 global Parquet source not divided into exact 8 per shard")
    return rows


def read_batches(key,limit=None):
    import pyarrow.parquet as pq
    import pyarrow.fs as fs
    import pyarrow.compute as pc
    # Public S3 anonymous, no credentials or user geography.
    s3=fs.S3FileSystem(anonymous=True,region="us-west-2")
    name=BUCKET+"/"+key
    file=s3.open_input_file(name)
    pf=pq.ParquetFile(file)
    cols=set(pf.schema_arrow.names)
    if not set(COLUMNS)<=cols:
        file.close()
        raise ValueError("Official Overture original parquet missing one required rail graph column")
    source_census_rows=pf.metadata.num_rows
    source_groups=pf.metadata.num_row_groups
    def generator():
        scanned=0
        try:
            for batch in pf.iter_batches(batch_size=BATCH,columns=COLUMNS,use_threads=True):
                scanned+=batch.num_rows
                selected=batch.filter(pc.equal(batch.column(batch.schema.get_field_index("subtype")),"rail"))
                yield scanned,batch.num_rows,selected,source_census_rows,source_groups
                if limit is not None and scanned>=limit:break
        finally:file.close()
    return generator()


def classify_and_measure(batch,counts,geod,read_all):
    from shapely import wkb
    for row in batch.to_pylist():
        rid=row.get("id")
        if not isinstance(rid,str) or not rid:
            raise ValueError("Physical source rail feature lacks unique GERS string ID")
        counts["rail_total_segments"]+=1
        cls=row.get("class") or "source_class_null"
        counts["class:"+str(cls)]+=1
        if row.get("names") and isinstance(row["names"],dict) and row["names"].get("primary"):
            counts["source_named_segments"]+=1
        if row.get("connectors"):
            counts["segments_with_connector_info"]+=1
            counts["original_connector_incidence"]+=len(row["connectors"])
            if not all(c.get("connector_id") for c in row["connectors"]):
                counts["missing_connector_id_entries"]+=1
        for src in row.get("sources") or ():
            if str(src.get("dataset")).lower()=="openstreetmap":
                orig=str(src.get("record_id") or "")
                counts["osm_way_reference_entries"]+=bool(re.fullmatch(r"w\d+@\d+",orig))
                counts["osm_other_reference_entries"]+=not bool(re.fullmatch(r"w\d+@\d+",orig))
        flag_values=set()
        for f in row.get("rail_flags") or []:
            if isinstance(f,dict):
                flag_values.update(f.get("values") or ())
        for flag in flag_values:counts["flag:"+str(flag)]+=1
        if read_all:
            geom=row.get("geometry")
            if not isinstance(geom,(bytes,bytearray,memoryview)):raise ValueError("Missing real rail centerline WKB")
            shape=wkb.loads(bytes(geom))
            if shape.geom_type!="LineString":raise ValueError("Unexpected source rail geometry not LINESTRING")
            line_m=abs(float(geod.geometry_length(shape)))
            if not math.isfinite(line_m) or line_m<0:
                raise ValueError("Original geodesic rail centerline length invalid")
            counts["source_geodesic_rail_metres"]+=line_m


def scan(mode,idx,out):
    import pyarrow.parquet as pq
    from pyproj import Geod
    source=inventory()
    paths=assigned_files(source,idx)
    if mode=="pilot":
        paths=paths[:1]
    counts=Counter()
    perfile=[]
    writer=None
    shard_path=Path(out).with_suffix(".rail.parquet")
    last_report_time=time.monotonic()
    try:
        for i,x in enumerate(paths):
            seen=0
            source_total=None
            source_groups=None
            for scanned,n,b,rows,groups in read_batches(x["key"],limit=PILOT_LIMIT_SOURCE_ROWS if mode=="pilot" else None):
                seen=scanned
                source_total=rows
                source_groups=groups
                counts["source_original_transport_rows_scanned"]+=n
                counts["source_rail_arrow_batches"]+=bool(b.num_rows)
                if b.num_rows:
                    classify_and_measure(b,counts,Geod(ellps="WGS84"),read_all=(mode=="full"))
                    if mode=="full":
                        if writer is None:
                            shard_path.parent.mkdir(parents=True,exist_ok=True)
                            writer=pq.ParquetWriter(shard_path,b.schema,
                              compression="zstd",compression_level=5,write_statistics=True)
                        writer.write_batch(b)
                if (time.monotonic()-last_report_time)>40:
                    print(json.dumps({"source_global_shard":idx,
                        "files_finished":i,"current_source_file":i+1,
                        "source_rows_scanned":counts["source_original_transport_rows_scanned"],
                        "rail_segments_scanned":counts["rail_total_segments"]}),flush=True)
                    last_report_time=time.monotonic()
            if source_total is None:raise ValueError("No original source Parquet batch")
            if mode=="full" and seen!=source_total:
                raise ValueError(f"Partial original source file {x['key']} {seen}!={source_total}")
            perfile.append({"source_key":x["key"],"source_parquet_bytes":x["bytes"],
                            "source_original_rows_reported_in_parquet":source_total,
                            "source_original_row_groups":source_groups,
                            "actual_source_rows_scanned":seen,
                            "source_file_fully_scanned":mode=="full"})
            if mode=="full" and counts["source_geodesic_rail_metres"]<=0:
                # A shard can in principle contain no rail rows due source
                # partition geography. Do not imply no-world rail.
                pass
        if writer is not None:writer.close()
        if mode=="full" and writer is None:
            raise ValueError("Full shard had no railway rows; inspect original partitioning")
        checksum=None
        bytes_written=0
        if writer is not None:
            bytes_written=shard_path.stat().st_size
            if bytes_written>MAX_RAIL_RECORD_BYTES:
                raise ValueError("Derived rail shard exceeds 1.2GB GitHub artifact bound")
            sha=hashlib.sha256()
            with shard_path.open("rb") as f:
                for chunk in iter(lambda:f.read(2**20),b""):sha.update(chunk)
            checksum=sha.hexdigest()
        stats={"version":"v357-original-global-rail-shard"+("-PILOT" if mode=="pilot" else ""),
           "source_Overture_release":RELEASE,
           "source_inventory_SHA256":source["full_object_inventory_sha256"],
           "source_all_global_parquet_objects_count":source["source_parquet_object_count"],
           "shard_index":idx,"shards_declared":SHARDS,
           "source_shard_expected_files":len(paths),
           "all_assigned_source_files_scanned":(mode=="full" and all(r["source_file_fully_scanned"] for r in perfile)),
           "original_source_files":perfile,
           "source_counts":dict(counts),
           "source_geodesic_rail_length_km":round(counts["source_geodesic_rail_metres"]/1000,3) if mode=="full" else None,
           "derived_source_rail_only_GeoParquet_file":str(shard_path.name) if writer else None,
           "derived_source_rail_only_GeoParquet_bytes":bytes_written,
           "derived_source_rail_only_GeoParquet_SHA256":checksum,
           "derived_rail_GeoParquet_contains_original_OdbL_source_data":writer is not None,
           "ODbL_attribution_source":"Overture Maps 2026-11-18.0 / OpenStreetMap contributors (ODbL)",
           "source_GERS_connectors_not_named_railway_lines":True,
           "the_shard_ALONE_is_not_global_census":True,
           "no_plant_community_similarity_or_warmth_calculated":True}
        Path(out).parent.mkdir(parents=True,exist_ok=True)
        Path(out).write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding="utf8")
        print(json.dumps({k:v for k,v in stats.items() if k!="original_source_files"},ensure_ascii=False),flush=True)
        return stats
    finally:
        if writer is not None:
            writer.close()


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--mode",choices=["pilot","full"],required=True)
    parser.add_argument("--shard",type=int,default=0)
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    scan(args.mode,args.shard,args.output)
