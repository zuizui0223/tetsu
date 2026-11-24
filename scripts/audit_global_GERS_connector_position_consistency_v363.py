#!/usr/bin/env python3
"""v363b true geodesic connector UUID positional source-integrity audit.

Source globally complete original Overture 2026-11-18, 128 Parquet source
files summarized as 16 original rail-only SHA-pinned GeoParquet shards
from v357; 2,054,855 rail GERS segments, not train operated lines.

Sample shared connector UUIDs independently of geometric distances by
SHA256(connector_id) first 4 bytes % 256==0; compare all occurrences of
that true UUID after interpolation to original source LineString.at.
Output no UUID, original railway geometry or source botanical locations.
"""
from __future__ import annotations
import argparse,hashlib,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path
from build_global_Overture_physical_connector_components_v361 import (
    validate_files,SOURCE_RELEASE,EXPECTED_INVENTORY_SHA
)
from audit_global_Overture_all_rail_parquet_inventory_v356 import RELEASE

MODULUS=256
EXPECTED_SOURCE_ALL_RAIL=2054855
MIN_SELECTED_CONNECTORS=1000
MAX_SELECTED_CONNECTORS=50000
SOURCE_ALLOWED_AT_EPS=1e-12
THRESHOLDS=(5,25,100,1000)

def pick_uuid(uuid):
    if not isinstance(uuid,str) or not uuid:return False
    first4=hashlib.sha256(uuid.encode("utf8")).digest()[:4]
    return int.from_bytes(first4,"big")%MODULUS==0

def true_WGS84_geodesic_at(line,at,geod):
    """Overture source at = travelled WGS84 ellipsoid length / full length.

    DO NOT use shapely.LineString.interpolate(normalized=True), because
    planar longitude/latitude coordinate lengths disagree with the
    geodesic distance fraction published by Overture, especially for
    curved track segments and high latitudes.
    """
    coords=list(line.coords)
    if len(coords)<2:
        raise ValueError("Source rail line lacks two original WGS84 vertices")
    if not (0.0<=at<=1.0):
        raise ValueError("Original at location is outside normalized source range")
    if at<=1e-12:return float(coords[0][0]),float(coords[0][1])
    if at>=1-1e-12:return float(coords[-1][0]),float(coords[-1][1])
    dist=[]
    total=0.
    for (x0,y0),(x1,y1) in zip(coords,coords[1:]):
        az,_,length=geod.inv(x0,y0,x1,y1)
        if not math.isfinite(length) or length<0:
            raise ValueError("Original source rail line WGS84 segment length invalid")
        dist.append((x0,y0,az,length))
        total+=length
    if total<=0:raise ValueError("Original source rail line has no true geodesic extent")
    goal=total*at
    previous=0.
    for x,y,az,length in dist:
        if goal<=previous+length:
            xnew,ynew,_=geod.fwd(x,y,az,min(length,max(0.0,goal-previous)))
            return float(xnew),float(ynew)
        previous+=length
    return float(coords[-1][0]),float(coords[-1][1])


def sampled_connector_positions(shard_paths):
    import pyarrow.parquet as pq
    from shapely import wkb
    from pyproj import Geod
    geod=Geod(ellps="WGS84")
    sampled={}
    counts=Counter()
    for path in shard_paths:
        source=pq.ParquetFile(path)
        for batch in source.iter_batches(batch_size=20000,columns=["id","connectors","geometry"],use_threads=True):
            for row in batch.to_pylist():
                counts["source_original_rail_segment_rows_scanned"]+=1
                conns=row.get("connectors") or []
                candidates=[]
                seen=set()
                for c in conns:
                    if not isinstance(c,dict):
                        counts["source_malformed_connector_incidence"]+=1
                        continue
                    cid=str(c.get("connector_id") or "")
                    if not cid:
                        counts["source_missing_connector_id"]+=1
                        continue
                    counts["source_original_total_connector_incidence"]+=1
                    if not pick_uuid(cid):continue
                    counts["source_sampled_connector_occurrences"]+=1
                    try:at=float(c["at"])
                    except (ValueError,TypeError,KeyError):
                        counts["source_sampled_invalid_at"]+=1;continue
                    if not math.isfinite(at) or not -SOURCE_ALLOWED_AT_EPS<=at<=1+SOURCE_ALLOWED_AT_EPS:
                        counts["source_sampled_invalid_at"]+=1;continue
                    at=min(1.,max(0.,at))
                    if cid in seen:
                        counts["source_sampled_duplicate_uuid_within_same_segment"]+=1
                        continue
                    seen.add(cid)
                    candidates.append((cid,at))
                if not candidates:continue
                shape=wkb.loads(bytes(row["geometry"]))
                if shape.geom_type!="LineString" or shape.is_empty:
                    raise ValueError("Original source GERS rail line geometry is malformed")
                for cid,at in candidates:
                    lon,lat=true_WGS84_geodesic_at(shape,at,geod)
                    if not all(math.isfinite(z) for z in (lon,lat)):
                        raise ValueError("Invalid original WGS84 source connector position")
                    if not (-180<=lon<=180 and -90<=lat<=90):
                        raise ValueError("Original coordinate out of WGS84 domain")
                    if cid not in sampled:
                        sampled[cid]=[lon,lat,1,0.0]
                    else:
                        previous=sampled[cid]
                        metres=abs(geod.inv(lon,lat,previous[0],previous[1])[2])
                        if not math.isfinite(metres):raise ValueError("Nonfinite source connector displacement")
                        previous[2]+=1
                        previous[3]=max(metres,previous[3])
    if counts["source_original_rail_segment_rows_scanned"]!=EXPECTED_SOURCE_ALL_RAIL:
        raise ValueError("All 2,054,855 worldwide original source rails were not read")
    if not MIN_SELECTED_CONNECTORS<=len(sampled)<=MAX_SELECTED_CONNECTORS:
        raise ValueError("Stable source connector sampling hash population unexpected")
    multi=[r for r in sampled.values() if r[2]>=2]
    if len(multi)<100:
        raise ValueError("Too few original repeated source connector IDs to QA global rail graph")
    bins={f"max_source_UUID_displacement_above_{threshold}_m":sum(r[3]>threshold for r in multi)
          for threshold in THRESHOLDS}
    degrees=Counter(r[2] for r in sampled.values())
    return {"original_source_all_rail_segments_audited":counts["source_original_rail_segment_rows_scanned"],
      "source_original_connector_incidence":counts["source_original_total_connector_incidence"],
      "source_selected_uuid_sha256_modulus":MODULUS,
      "original_source_connector_uuid_selected":len(sampled),
      "original_source_selected_connector_incidence":counts["source_sampled_connector_occurrences"],
      "selected_uuid_at_least_two_original_rail_segments":len(multi),
      "selected_uuid_single_original_segment":len(sampled)-len(multi),
      "source_selected_invalid_connector_at":counts["source_sampled_invalid_at"],
      "source_selected_duplicate_UUID_within_same_segment":counts["source_sampled_duplicate_uuid_within_same_segment"],
      "source_repeated_connector_UUID_maximum_geodesic_separation_m":
         round(max(r[3] for r in multi),3),
      "source_repeated_connector_geodesic_displacement_median_m":
         round(statistics.median(r[3] for r in multi),3),
      "source_uuid_mismatch_counts_exceeding_original_25_100_and_1000m":bins,
      "sampled_UUID_incidence_degree_quantiles_NOT_route_degrees":dict(sorted(
         (str(k),v) for k,v in degrees.items() if k<=8)),
      "source_original_UUIDs_source_segment_names_GPS_not_exported":True,
      "physical_network_connector_consistency_not_train_route_operation_proof":True}

def run(source_dir,census_path):
    original=json.loads(Path(census_path).read_text(encoding="utf8"))
    paths=validate_files(source_dir,original)
    result=sampled_connector_positions(paths)
    return {"version":"v363b-true-WGS84-geodesic-global-GERS-connector-UUID-positional-QA",
      "original_global_rail_source_release":SOURCE_RELEASE,
      "original_genuine_full_128_Overture_source_inventory_SHA256":EXPECTED_INVENTORY_SHA,
      "original_global_full_16_rail_only_GeoParquet_SHA_verified":True,
      "stable_selection":"SHA256(connector_id)[:4] interpreted as 32-bit int mod 256 == 0; not source data outcome or botanist site selection",
      "correct_original_Overture_at_interpretation":"WGS84 geodetic distance fraction along original line segments, NOT planar lonlat shapely normalized interpolate",
      "audit":result,
      "NO_worldwide_botanical_Jaccard_or_train_traffic_result":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",required=True)
    p.add_argument("--census",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    obj=run(a.source_dir,a.census)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(obj,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
