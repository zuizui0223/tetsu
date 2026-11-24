#!/usr/bin/env python3
"""v349b: Actual 2025 iNaturalist plant detections × Overture-2026 railway
centerlines in ALL SIX v348 predeclared continents.

This is an observational *railway proximity* source pilot, not plant relevé
at trackbed. Same frozen d+GPS uncertainty≤80m and background300–1000m
from v348, original species labels and exact positions remain runner-private.
Uses Overture global GeoParquet instead of overloaded public Overpass.
"""
from __future__ import annotations
import argparse,hashlib,importlib.metadata,json,math,sys
from pathlib import Path
from collections import Counter

import pilot_self_extracted_six_continent_rail_flora_v348 as core
from audit_Overture_global_rail_segments_v349 import eligible_segment

CITIES=core.CITIES
MAX_SEGMENTS=300000
SOURCE_RAIL_BUFFER_M=1100.0

def source_rail_query_bbox(lat,lon):
    """Expand source rail coverage beyond fixed plant window by >1km.

    Prevents falsely classifying edge observations as far from a railway
    because the closest track lies just outside the source download bbox.
    """
    a,b,c,d=core.bounds(lat,lon)
    h=SOURCE_RAIL_BUFFER_M/111320.0
    lon_h=h/max(0.05,math.cos(math.radians(lat)))
    return (b-lon_h,a-h,d+lon_h,c+h)


def rail_provider_overture(lat0,lon0):
    from overturemaps import record_batch_reader
    from shapely import wkb
    city=next((c for c in CITIES if c[3]==lat0 and c[4]==lon0),None)
    if city is None:
        raise ValueError("Source city not in frozen original v348 sample frame")
    box=source_rail_query_bbox(lat0,lon0)
    errors=[]
    # Source access retry is *not* a new source selection or biological
    # outcome-driven rule. The Overture STAC metadata can return HTTP429.
    for attempt in (1,2):
        try:
            reader=record_batch_reader("segment",bbox=box)
            if reader is None:raise ValueError("No original Overture segment source")
            table=reader.read_all()
            if table.num_rows>MAX_SEGMENTS:raise ValueError("Source geometry row cap exceeded")
            ways=[];counts=Counter()
            for row in table.to_pylist():
                good,why=eligible_segment(row)
                if not good:
                    counts[why]+=1;continue
                geom=wkb.loads(bytes(row["geometry"]))
                if geom.geom_type=="LineString":geoms=[geom]
                elif geom.geom_type=="MultiLineString":geoms=list(geom.geoms)
                else:
                    counts["other_geometry"]+=1
                    continue
                for ln in geoms:
                    coords=[(float(y),float(x)) for x,y in ln.coords]
                    if len(coords)>=2:ways.append(coords)
            return ways,{"provider":"Overture_global_transportation_2026_rail",
                "source_overall_transportation_feature_count":table.num_rows,
                "source_rail_ways_eligible":len(ways),
                "source_rail_segments_eligible":sum(len(q)-1 for q in ways),
                "source_rail_class_source_def":"standard_gauge/narrow_gauge/unknown, no non-rail/tram/subway, no bridge/tunnel/disused",
                "source_rail_data_download_buffer_around_plant_window_m":SOURCE_RAIL_BUFFER_M,
                "source_no_track_shapes_republished":True,
                "source_attempt_n":attempt,
                "source_error_first_attempt":errors}
        except Exception as e:
            errors.append({"attempt":attempt,"type":type(e).__name__,"short_message":str(e)[:180]})
            # No back-to-back metadata hits when STAC requests are rate limited.
            if attempt==1:
                import time;time.sleep(20)
    return [],{"provider":"Overture_failed_not_zero_biology",
               "source_rail_ways_eligible":0,"source_rail_segments_eligible":0,
               "source_errors":errors}

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",required=True)
    args=parser.parse_args()
    core.source_rail=rail_provider_overture
    report=core.run()
    report["version"]="v349b-six-continents-original-Overture-rail-vs-iNat-plant-proximity"
    report["alternative_network_source"]="official Overture global rail segments (not API Overpass)"
    report["source_geometry_window_buffer_m"]=SOURCE_RAIL_BUFFER_M
    report["original_Overture_client_version"]=importlib.metadata.version("overturemaps")
    report["source_month"]="2026-11-18 candidate Overture release, actual STAC latest must be source pinned for final claims"
    report["country_city_centres_are_fixed_sampling_windows_not_probability_global_representative"]=True
    report["no_evidenced_rail_microclimate_or_plant_temperature_effect"]=True
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
