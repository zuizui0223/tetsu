#!/usr/bin/env python3
"""Bounded OSM road HISTORY pilot (2017 vs 2024) near Utsunomiya LRT.

Uses open HeiGIT ohsome API /v1/elements/geometry for road ways mapped at
specific historical dates. Road mapped in OSM in 2017 is *positive*
evidence of an OSM mapped road; missing OSM way does NOT establish
physical road absence or its construction date.

Official MLIT 2024 N13 road centres are used as a contemporaneous mapping
completeness reference. MLIT N05 noncommercial railway geometry, OSM
feature geometries, locations and identifiers are transient, not exported.

No seed, greenness/plant-cover, mowing, wind or railway causal effect.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import tempfile

import requests
from shapely.geometry import Point,shape
from shapely.ops import transform as transform_geom
from shapely.strtree import STRtree
from pyproj import Transformer

from audit_n05_historical_railway import bounded_zip_download
from measure_utsunomiya_s2_ndvi_pixels import (
    get_lrt_lines,build_fixed_pairs
)
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_verge_pixel_eligibility import near_points_from_140m
from tetsu.rail_road_landbalance import build_n13_road_index

OHSOME_BASE="https://api.ohsome.org/v1"
OHSOME_ROADS=OHSOME_BASE+"/elements/geometry"
OHSOME_META=OHSOME_BASE+"/metadata"
OSM_QUERY_FILTER=(
    "type:way and highway in (motorway,motorway_link,trunk,trunk_link,"
    "primary,primary_link,secondary,secondary_link,tertiary,tertiary_link,"
    "residential,service,unclassified,living_street,road)"
)
DATES={"2017":"2017-12-31","2024":"2024-09-30"}
ANCHORS=(4,12,20,28,36,44)
LAT_HALF=0.020
LON_HALF=0.028
MAX_RESPONSE_BYTES=18_000_000
NEAR_RAIL_OFFSETS=(20,35)
DISTANCE_THRESHOLD_M=25.


def date_aoi(lon,lat):
    return (
        round(lon-LON_HALF,6),round(lat-LAT_HALF,6),
        round(lon+LON_HALF,6),round(lat+LAT_HALF,6)
    )


def year_source_geojson(client,bbox,year):
    if year not in DATES:
        raise ValueError("Unapproved OSM snapshot year")
    if not (139<=bbox[0]<bbox[2]<=141 and 35<=bbox[1]<bbox[3]<=38):
        raise ValueError("Out-of-study-area historical road request")
    args={
        "bboxes":",".join(f"{v:.6f}" for v in bbox),
        "time":DATES[year],
        "filter":OSM_QUERY_FILTER,
        "clipGeometry":"true",
    }
    with client.post(
        OHSOME_ROADS,data=args,stream=True,timeout=(20,125)
    ) as resp:
        resp.raise_for_status()
        body=bytearray()
        for chunk in resp.iter_content(chunk_size=1024*512):
            body.extend(chunk)
            if len(body)>MAX_RESPONSE_BYTES:
                raise ValueError("Historical OSM single response too large")
    try:
        doc=json.loads(bytes(body))
    except ValueError as exc:
        raise ValueError("Historical OSM returned invalid GeoJSON") from exc
    if not isinstance(doc,dict) or doc.get("type")!="FeatureCollection":
        raise ValueError("Ohsome failed to return OSM FeatureCollection")
    features=doc.get("features")
    if not isinstance(features,list):
        raise ValueError("Ohsome FeatureCollection features missing")
    return features,hashlib.sha256(bytes(body)).hexdigest(),len(body)


def lines_from_osm(features,utm):
    road=[]
    counts=Counter()
    for feature in features:
        geom=feature.get("geometry")
        if not geom:
            counts["no_geometry"]+=1
            continue
        try:
            parsed=shape(geom)
        except (ValueError,TypeError) as exc:
            raise ValueError("Invalid historical road feature geometry") from exc
        if parsed.geom_type=="LineString":
            options=[parsed]
        elif parsed.geom_type=="MultiLineString":
            options=list(parsed.geoms)
        else:
            counts["non_road_line_geometry"]+=1
            continue
        for line in options:
            if line.is_empty or len(line.coords)<2:
                continue
            metro=transform_geom(utm.transform,line)
            if metro.length>=1:
                road.append(metro)
        counts["OSM_feature_rows"]+=1
    return road,dict(counts)


def nearest_metric(points,road_lines):
    if not road_lines:
        # No 2017 OSM road ways recorded near the site means UNKNOWN
        # physical road exposure, not a 0 count of real roads.
        return {k:None for k in points}
    index=STRtree(road_lines)
    return {
        key:float(Point(*loc).distance(road_lines[int(index.nearest(Point(*loc)))]))
        for key,loc in points.items()
    }


def pick_study_sites():
    n05=bounded_zip_download()
    lines,_=get_lrt_lines(n05)
    points,pairs,_=build_fixed_pairs(lines,48)
    selected=[]
    for anchor in ANCHORS:
        candidates=[(near,far) for near,far in pairs
                    if int(near.split("_",1)[0][1:])==anchor]
        if candidates:
            selected.append(sorted(candidates)[0])
    if len(selected)<4:
        raise ValueError("Fewer than 4 selected history pilot anchors")
    sites=[]
    ll=Transformer.from_crs(32654,4326,always_xy=True)
    for oldnear,far in selected:
        study_pts={}
        for offset in NEAR_RAIL_OFFSETS:
            pp,pair=near_points_from_140m(points,[(oldnear,far)],offset)
            near_id=pair[0][0]
            study_pts[near_id]=pp[near_id]
        study_pts[far]=points[far]
        lon,lat=ll.transform(*points[oldnear])
        sites.append({
            "anchor_index":int(oldnear.split("_",1)[0][1:]),
            "bbox":date_aoi(lon,lat),
            "points":study_pts
        })
    return sites,hashlib.sha256(n05).hexdigest()


def run():
    samples,n05_sha=pick_study_sites()
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-LRT-2017-openstreetmap-road-history-v1.4"})
    metadata=client.get(OHSOME_META,timeout=(15,40))
    metadata.raise_for_status()
    info=metadata.json()
    temporal=info.get("extractRegion",{}).get("temporalExtent",{})
    if not str(temporal.get("fromTimestamp","")).startswith("2007"):
        raise ValueError("Ohsome backend 2017 temporal coverage unverified")
    utm=Transformer.from_crs(4326,32654,always_xy=True)
    result={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "OSM_api":OHSOME_ROADS,
        "OSM_snapshot_dates":DATES,
        "OSM_way_filter":OSM_QUERY_FILTER,
        "OSM_attribution":"© OpenStreetMap contributors / ODbL, via ohsome API, HeiGIT",
        "OSM_source_api_version":info.get("apiVersion"),
        "OSM_source_temporal_extent":temporal,
        "N05_source_SHA256":n05_sha,
        "N05_noncommercial_source_geometry_not_redistributed":True,
        "n_predeclared_rail_anchors":len(samples),
        "bbox_area":"~5km × 4km around each selected rail anchor, actual OSM geometries clipped",
        "N13_2024_matched_reference":{},
        "OSM_road_snapshot_audits":{},
        "compare_snapshots_and_current_N13":{},
        "failures":{},
        "limits":[
            "2017 OpenStreetMap geometry describes ROAD WAYS MAPPED BY 2017-12-31, not all roads actually physically existing then.",
            "Missing roads in 2017 OSM but mapped in 2024 might reflect mapping edits, imagery tracing or changed contributor coverage, not railway construction.",
            "OSM service and local current MLIT N13 road geometry use different classification completeness, road types and survey protocols.",
            "2017 preconstruction roads are not directly operator-certified and a mapped highway may omit carriageway width and tram sharing.",
            "Urban light rail is largely co-located with roads; 10m Sentinel pixel cannot separate tram ballast vs street tree/gardens.",
            "No causal effect, train-carried seeds, vegetation management, adaptation, species arrival, wind or greenness inferred.",
            "All source geometries, OSM road IDs, track coordinates and individual observation points are omitted from public artifacts."
        ]
    }
    with tempfile.TemporaryDirectory() as root:
        nearest_n13,n13_src=build_n13_road_index(download_dir=Path(root))
        result["N13_2024_matched_reference"]=n13_src
        source_blocked_403=False
        for sample in samples:
            aid=sample["anchor_index"]
            pointset=sample["points"]
            contemporary={key:nearest_n13(coord) for key,coord in pointset.items()}
            successful={}
            for year in (() if source_blocked_403 else DATES):
                try:
                    features,sha,nbytes=year_source_geojson(client,sample["bbox"],year)
                    lines,counts=lines_from_osm(features,utm)
                    nearest=nearest_metric(pointset,lines)
                    successful[year]=nearest
                    result["OSM_road_snapshot_audits"][f"anchor{aid:02d}_{year}"]={
                        "n_source_features":len(features),
                        "n_projected_road_lines":len(lines),
                        "source_sha256":sha,"source_bytes":nbytes,
                        "quality":counts,
                        "nonzero_ways_in_query":bool(lines),
                    }
                except requests.HTTPError as exc:
                    response=getattr(exc,"response",None)
                    if response is not None and response.status_code==403:
                        source_blocked_403=True
                        result["failures"][f"anchor{aid:02d}_{year}"]=(
                            "403 source extraction unavailable or API-key required"
                        )
                        break  # Do not issue 11 more known-blocked requests.
                    result["failures"][f"anchor{aid:02d}_{year}"]=(
                        f"{type(exc).__name__}: {str(exc)[:200]}"
                    )
                except Exception as exc:
                    result["failures"][f"anchor{aid:02d}_{year}"]=(
                        f"{type(exc).__name__}: {str(exc)[:200]}"
                    )
            for key in pointset:
                group=("rail20" if key.endswith("_near20m") else
                       "rail35" if key.endswith("_near35m") else "far1600")
                x2017=successful.get("2017",{}).get(key)
                x2024=successful.get("2024",{}).get(key)
                n13=contemporary[key]
                result["compare_snapshots_and_current_N13"].setdefault(group,[]).append({
                    "ohsome_2017_road_le25m":x2017 is not None and x2017<=DISTANCE_THRESHOLD_M,
                    "ohsome_2024_road_le25m":x2024 is not None and x2024<=DISTANCE_THRESHOLD_M,
                    "n13_2024_road_le25m":n13<=DISTANCE_THRESHOLD_M,
                    "ohsome_2017_response_success":key in successful.get("2017",{}),
                    "ohsome_2024_response_success":key in successful.get("2024",{}),
                    "ohsome_2017_nearest_road_known":x2017 is not None,
                    "ohsome_2024_nearest_road_known":x2024 is not None,
                })
            print(json.dumps({
                "anchor_number":aid,
                "2017_osm_road_features":result["OSM_road_snapshot_audits"].get(f"anchor{aid:02d}_2017",{}).get("n_source_features"),
                "2024_osm_road_features":result["OSM_road_snapshot_audits"].get(f"anchor{aid:02d}_2024",{}).get("n_source_features"),
                "failed_queries":sum(key.startswith(f"anchor{aid:02d}") for key in result["failures"])
            },ensure_ascii=False),flush=True)
    # Remove individual location-level classification; only counts.
    group_audit={}
    for group,entries in result.pop("compare_snapshots_and_current_N13").items():
        group_audit[group]={
            "n_sample_sites":len(entries),
            "2017_road_way_mapped_within25m":(
                sum(x["ohsome_2017_road_le25m"] for x in entries)
                if any(x["ohsome_2017_response_success"] for x in entries)
                else None
            ),
            "2024_road_way_mapped_within25m":(
                sum(x["ohsome_2024_road_le25m"] for x in entries)
                if any(x["ohsome_2024_response_success"] for x in entries)
                else None
            ),
            "2024_N13_road_within25m":sum(x["n13_2024_road_le25m"] for x in entries),
            "2017_OSM_queries_success":sum(x["ohsome_2017_response_success"] for x in entries),
            "2024_OSM_queries_success":sum(x["ohsome_2024_response_success"] for x in entries),
            "2017_OSM_nearest_line_known":sum(x["ohsome_2017_nearest_road_known"] for x in entries),
            "2024_OSM_nearest_line_known":sum(x["ohsome_2024_nearest_road_known"] for x in entries),
            "missing_from_OSM_not_physical_absence":True,
        }
    result["grouped_road_history_comparison"]=group_audit
    result["historical_OSM_source_status"]=(
        "blocked_403_no_2017_road_evidence" if source_blocked_403
        else "some_historical_query_responses" if
        result["OSM_road_snapshot_audits"] else
        "no_historical_query_responses"
    )
    result["403_known_provider_context"]=(
        "HeiGIT has restricted legacy ohsome v1 /elements/geometry; "
        "new v2 extraction requires authenticated API key. Do not retry "
        "anonymous extraction or treat unavailable data as mapped-road absence."
    )
    print(json.dumps({"history_summary":group_audit,"errors":result["failures"]},
                     ensure_ascii=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/lrt_2017_osm_road_history_vs_2024_N13_v14.json")
    a=p.parse_args()
    report=run()
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    if report["historical_OSM_source_status"]=="blocked_403_no_2017_road_evidence":
        print("Historical OSM extraction BLOCKED by provider 403; 2017 road exposure is UNKNOWN, not zero.")
        return
    if not any(v.get("nonzero_ways_in_query")
               for v in report["OSM_road_snapshot_audits"].values()):
        raise SystemExit("No historical OSM road data confirmed; no 2017 road absence claim")


if __name__=="__main__":
    main()
