#!/usr/bin/env python3
"""v357: HARD fail-closed aggregate of ALL original Overture railway Parquet scan shards.

A single skipped original file or duplicate Parquet source key forbids
a "worldwide railway all extracted" statement. Aggregates original source
segments and all classes; actual named OSM route relations are NOT included
in Overture and require independent source and coverage audit.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path
from audit_global_Overture_all_rail_parquet_inventory_v356 import list_original_parquet_objects,RELEASE
from extract_entire_world_Overture_rail_parquet_v357 import (
    SHARDS,EXPECTED_SOURCE_PARQUET_OBJECTS,EXPECTED_SOURCE_BYTES,
    EXPECTED_SOURCE_INVENTORY_SHA
)

PUBLISHER_RAIL_LENGTH_KM=2_892_198

def combine(parts,source_manifest):
    if len(parts)!=SHARDS:
        raise ValueError("Original full world railway extraction missing one of 16 shards")
    if (source_manifest.get("source_parquet_object_count")!=EXPECTED_SOURCE_PARQUET_OBJECTS
        or source_manifest.get("source_parquet_compressed_bytes")!=EXPECTED_SOURCE_BYTES
        or source_manifest.get("full_object_inventory_sha256")!=EXPECTED_SOURCE_INVENTORY_SHA
        or not source_manifest.get("source_listing_all_pages_consumed")):
        raise ValueError("Original world transportation file source inventory incomplete/drifted")
    expected={x["key"]:x for x in source_manifest["all_source_parquet_objects"]}
    seen={}
    total=Counter()
    shards=[]
    for shard_index in range(SHARDS):
        p=parts[shard_index]
        if p.get("version")!="v357-original-global-rail-shard" or p.get("shard_index")!=shard_index:
            raise ValueError("Global source shard role/order not verified")
        if (p.get("shards_declared")!=SHARDS or p.get("source_inventory_SHA256")!=EXPECTED_SOURCE_INVENTORY_SHA
            or p.get("source_Overture_release")!=RELEASE
            or not p.get("all_assigned_source_files_scanned")
            or len(p.get("original_source_files",[]))!=8):
            raise ValueError("Original worldwide rail source shard incomplete")
        for x in p["original_source_files"]:
            k=x["source_key"]
            if k not in expected or k in seen or (x["source_parquet_bytes"]!=expected[k]["bytes"]
                or not x["source_file_fully_scanned"]
                or x["actual_source_rows_scanned"]!=x["source_original_rows_reported_in_parquet"]):
                raise ValueError("Original transport source file duplicate/missing/partial scan/size mismatch")
            seen[k]=shard_index
        counts=p["source_counts"]
        if counts["source_original_transport_rows_scanned"]!=sum(
             x["source_original_rows_reported_in_parquet"] for x in p["original_source_files"]):
            raise ValueError("Rail shard original batch source row census mismatch")
        n=counts.get("rail_total_segments",0)
        if n>0 and (not p.get("derived_source_rail_only_GeoParquet_SHA256")
                    or p.get("derived_source_rail_only_GeoParquet_bytes",0)<=0):
            raise ValueError("Original source rail rows extracted but rail-only Parquet file not preserved")
        for key,value in counts.items():
            if isinstance(value,(int,float)) and math.isfinite(value):total[key]+=value
        shards.append({"shard_index":shard_index,
          "original_transport_source_files_n":len(p["original_source_files"]),
          "original_transport_source_rows":counts["source_original_transport_rows_scanned"],
          "original_rail_segments":n,
          "original_rail_geodesic_km":p["source_geodesic_rail_length_km"],
          "source_rail_only_Parquet_artifact_filename":p["derived_source_rail_only_GeoParquet_file"],
          "source_rail_Parquet_sha256":p["derived_source_rail_only_GeoParquet_SHA256"],
          "source_rail_Parquet_bytes":p["derived_source_rail_only_GeoParquet_bytes"]})
    if set(seen)!=set(expected):
        raise ValueError(f"Worldwide rail original source listing completeness failure: missing {len(set(expected)-set(seen))}")
    total_rail=total.get("rail_total_segments",0)
    km=total.get("source_geodesic_rail_metres",0.0)/1000
    if not 1_000_000<=total_rail<=100_000_000:
        raise ValueError("Global rail source census implausible original rail segment count")
    if not 2_000_000<=km<=4_000_000:
        raise ValueError("Global rail source geodesic length implausible; cannot declare all railway network")
    rel=(km-PUBLISHER_RAIL_LENGTH_KM)/PUBLISHER_RAIL_LENGTH_KM
    if abs(rel)>0.08:
        raise ValueError("Global official 2026-11-18 full rail length reference differs by >8%")
    byclass={key.removeprefix("class:"):val for key,val in total.items() if key.startswith("class:")}
    if sum(byclass.values())!=total_rail:
        raise ValueError("Original all rail class grouping does not sum to all original rail segments")
    return {"version":"v357-original-entire-world-Overture-2026-all-rail-scan-complete",
      "source_release":RELEASE,
      "source_original_rail_features_entire_release_n":int(total_rail),
      "source_real_geodesic_all_rail_km":round(km,3),
      "source_reported_publisher_rail_subtype_length_km":PUBLISHER_RAIL_LENGTH_KM,
      "source_geodesic_vs_publisher_relative_delta":round(rel,8),
      "source_all_segment_source_objects_complete":len(seen),
      "source_all_source_transportation_rows_scanned":int(total["source_original_transport_rows_scanned"]),
      "original_S3_source_Parquet_inventory_sha256":EXPECTED_SOURCE_INVENTORY_SHA,
      "original_rail_class_feature_count":dict(sorted(byclass.items())),
      "source_all_rail_named_feature_n":int(total["source_named_segments"]),
      "source_rail_features_with_physical_connector_info":int(total["segments_with_connector_info"]),
      "source_total_rail_connector_incidence":int(total["original_connector_incidence"]),
      "original_osm_way_source_ref_incidence":int(total["osm_way_reference_entries"]),
      "original_osm_non_way_source_ref_incidence":int(total["osm_other_reference_entries"]),
      "rail_flag_feature_n":{key.removeprefix("flag:"):val for key,val in total.items() if key.startswith("flag:")},
      "source_shards":shards,
      "entire_original_rail_subtype_includes_subway_tram":True,
      "physical_graph_components_NOT_equivalent_to_named_railway_lines":True,
      "raw_rail_only_Parquet_data_available_as_GitHub_Actions_artifacts_not_committed":True,
      "no_plant_community_similarity_estimated":True,
      "no_worldwide_specific_active_train_lines_or_population_effects_proved":True,
      "source_license":"Overture Maps, OSM contributors; ODbL"}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    import boto3
    source=list_original_parquet_objects()
    parts={}
    for i in range(SHARDS):
        path=Path(args.source_dir,f"v357_shard_{i}.json")
        if not path.is_file():
            raise FileNotFoundError("Missing complete original Overture rail source shard "+str(i))
        parts[i]=json.loads(path.read_text(encoding="utf8"))
    obj=combine(parts,source)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({k:v for k,v in obj.items() if k!="source_shards"},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
