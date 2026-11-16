#!/usr/bin/env python3
"""Find USA-NPN registered sites near rail BEFORE selecting plant taxa.

Uses ALL registered station coordinates in two prespecified states. This is
a sampling-frame feasibility audit and does not imply stations have 2024
observations, that plants grow beside tracks, or railway-caused phenology.
Station coordinates/IDs omitted from published reports.
"""
from __future__ import annotations
import argparse
import collections
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import requests
from audit_npn_rail_proximity import (SESSION,get_rail_features,bin_for_distance,
                                      UTM_EPSG)

STATIONS_URL="https://services.usanpn.org/npn_portal/stations/getAllStations.json"
STATES=("MD","AZ")


def registered_sites(state):
    if state not in STATES: raise ValueError("state not in predeclared audit")
    r=SESSION.get(STATIONS_URL,params={"state_code":state},timeout=(15,100))
    r.raise_for_status()
    payload=r.json()
    if not isinstance(payload,list):
        raise ValueError(f"NPN station list returned {type(payload).__name__}, not list")
    return payload


def pick_field(row,possible):
    for name in possible:
        if name in row and row[name] not in (None,"","-9999"):
            return row[name]
    return None


def summarize_sites(raw,dist_func):
    keys=sorted({k for row in raw if isinstance(row,dict) for k in row})
    counts=collections.Counter()
    unique_sites=set()
    by_distance=[]
    for row in raw:
        if not isinstance(row,dict):
            counts["invalid_object"]+=1
            continue
        sid=pick_field(row,("station_id","site_id","stationId","id"))
        lat=pick_field(row,("latitude","lat","station_latitude"))
        lon=pick_field(row,("longitude","long","lon","station_longitude"))
        try:
            site=str(sid)
            lat=float(lat)
            lon=float(lon)
            if sid is None or not math.isfinite(lat) or not math.isfinite(lon):
                raise ValueError
            if not (-125<=lon<=-65 and 24<=lat<=50):
                raise ValueError
        except (TypeError,ValueError):
            counts["missing_or_invalid_latlong_or_id"]+=1
            continue
        if site in unique_sites:
            counts["duplicate_station_id"]+=1
            continue
        unique_sites.add(site)
        distance=dist_func(lon,lat)
        by_distance.append(distance)
        counts[bin_for_distance(distance)]+=1
    # Qualifies geographic only; never claim all these had repeat observations.
    return {
        "n_stations_returned":len(raw),
        "n_valid_georeferenced_unique_registered_stations":len(unique_sites),
        "distance_bins":{label:int(counts[label]) for label in
                         ("0-250m","250-1000m","1000m_plus")},
        "invalid_or_duplicated_station_records":{
            k:v for k,v in counts.items() if k not in
            ("0-250m","250-1000m","1000m_plus")},
        "source_fields":keys,
        "min_distance_to_rail_m":round(min(by_distance),1) if by_distance else None,
        "note":"Registered station != active 2024 monitored phenophase site."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/usa_npn_registered_rail_sites.json")
    args=p.parse_args()
    report={"created_utc":datetime.now(timezone.utc).isoformat(),
            "method":"registered NPN site census first, no plant taxon selection",
            "source_station_api":STATIONS_URL,
            "source_fra":"https://services.arcgis.com/xOi1kZaI0eWDREZv/ArcGIS/rest/services/NTAD_North_American_Rail_Network_Lines/FeatureServer/0",
            "states":{},"failures":{}}
    for state in STATES:
        try:
            rows=registered_sites(state)
            nearest,fra=get_rail_features(state)
            result=summarize_sites(rows,nearest)
            result["rail_source"]=fra
            report["states"][state]=result
            print(json.dumps({"state":state,
                              "all_sites":result["n_stations_returned"],
                              "near_far":result["distance_bins"]},
                             ensure_ascii=False),flush=True)
        except Exception as e:
            report["failures"][state]=f"{type(e).__name__}: {str(e)[:350]}"
    report["caveats"]=[
        "No query here confirms which stations have 2024 repeated plant observations.",
        "Distance from registered site coordinate to 2026 FRA rail GIS does not mean vegetation grows on rail property.",
        "Rail vintage 2026 follows phenology 2024. An independently verified historical exposure would strengthen any test.",
        "No access to mowing history or train-gust intensity. This is not evidence for plasticity or local adaptation.",
        "Some public phenology sites may be gardens, visitor centers, urban parks or historic inactive stations.",
    ]
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    if not report["states"]:raise SystemExit("No valid NPN station census for either state")
    print(json.dumps({"success_states":list(report["states"]),"failures":report["failures"]}))


if __name__=="__main__":
    main()
