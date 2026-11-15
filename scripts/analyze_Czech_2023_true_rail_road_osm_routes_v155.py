#!/usr/bin/env python3
"""v15.5 historical Czech OSM RAIL vs ROAD genuine physical-network routes,
for conservative exact 35/39 Czech source original botanical station names.

Source author: Kutlvašr et al (2024), Preslia DOI
10.23855/preslia.2024.247; 2020-22 source plant lists.
2026-08 official SR70-derived STATION CENTROIDS are NOT original plant
plot GPS or historical freight station operation. Exact name only.
Physical OSM country network: Geofabrik Czech Republic 2023-01-01 PBF,
not a source freight traffic flow, true operating schedule or timed
train seed dispersal.
OSM way original ID adjacency is preserved: no false rail crossings.
Full-country shortest route can be truncated at international borders.
Underlying PBF/source XLSX/derived SR70 not republished, only aggregate
source QC and (internal) site-pair features for subsequent analysis.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from itertools import combinations
from pathlib import Path
from statistics import median

from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import (
    original_35_station_sites,real_site_pairs
)
from analyze_Slovenia_true_OSM_rail_road_paths_v150 import (
    read_osm_network,snap_sites,pair_paths,MAIN_RAIL_MAX_SNAP_KM,
    MAIN_ROAD_MAX_SNAP_KM,MIN_INTERSITE_GEO_KM
)

SOURCE_OSM="https://download.geofabrik.de/europe/czech-republic-230101.osm.pbf"
SOURCE_OSM_DATE="2023-01-01"
EXPECTED_SOURCE_PBF_BYTES=816159635
EXPECTED_STATIONS=35
EXPECTED_PAIRS=595

def stream_sha256(path,buf_size=1024*1024*8):
    sha=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(buf_size),b""):
            sha.update(block)
    return sha.hexdigest()

def graph_station_source(sites,pbf,return_pairs=False):
    from statistics import mean
    from math import isfinite,log
    import numpy as np
    size=Path(pbf).stat().st_size
    if size!=EXPECTED_SOURCE_PBF_BYTES:
        raise ValueError("Geofabrik original 2023-01 Czech country PBF bytes do not match published original archive")
    sha=stream_sha256(pbf)
    rail,road=read_osm_network(pbf)
    q={}
    for name,raw,maxdist in (
        ("rail",rail,MAIN_RAIL_MAX_SNAP_KM),
        ("road",road,MAIN_ROAD_MAX_SNAP_KM),
    ):
        graph,coords,component,graphqc=raw.finalize()
        snaps=snap_sites(coords,sites,maxdist)
        routes=pair_paths(graph,snaps)
        all_pairs=list(combinations(range(len(sites)),2))
        q[name]={
          "source_graph":graphqc,
          "n_actual_author_stations_with_matched_2026_point":len(sites),
          "n_source_station_centroids_connected_to_osm_node":sum(
                x["on_supported_source_network"] for x in snaps),
          "n_source_station_centroids_outside_fixed_snap_threshold":sum(
                not x["on_supported_source_network"] for x in snaps),
          "source_median_centroid_to_nearest_network_node_km":median(
                x["distance_km"] for x in snaps),
          "source_maximum_centroid_to_nearest_network_node_km":max(
                x["distance_km"] for x in snaps),
          "n_osm_physically_connected_station_pairs":sum(
                isfinite(routes[i,j]) for i,j in all_pairs),
          "original_predeclared_max_snap_km":maxdist,
          "snaps":snaps,"routes":routes,
        }
        del graph,coords,component,raw
    from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
    pairs=[]
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i],sites[j]
        geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
        x=float(q["rail"]["routes"][i,j])
        y=float(q["road"]["routes"][i,j])
        good=bool(geo>=MIN_INTERSITE_GEO_KM and isfinite(x) and isfinite(y)
                  and x>0 and y>0 and
                  q["rail"]["snaps"][i]["on_supported_source_network"] and
                  q["rail"]["snaps"][j]["on_supported_source_network"] and
                  q["road"]["snaps"][i]["on_supported_source_network"] and
                  q["road"]["snaps"][j]["on_supported_source_network"])
        pairs.append({
          "i":i,"j":j,"source_same_stretch":a["section"]==b["section"],
          "source_station_geodesic_km":geo,
          "network_route_comparison_eligible":good,
          "rail_graph_source_shortest_path_km":x if isfinite(x) else None,
          "road_graph_source_shortest_path_km":y if isfinite(y) else None,
          "rail_over_road_log_distance":math.log(x/y) if good else None,
          "source_both_rail_snap_km":sum(q["rail"]["snaps"][z]["distance_km"] for z in (i,j)),
          "source_both_road_snap_km":sum(q["road"]["snaps"][z]["distance_km"] for z in (i,j)),
        })
    if len(pairs)!=EXPECTED_PAIRS:
        raise ValueError("Czech original field 35 sites must yield 595 physical pair opportunities")
    supported=[p for p in pairs if p["network_route_comparison_eligible"]]
    ratios=[p["rail_graph_source_shortest_path_km"]/p["road_graph_source_shortest_path_km"]
            for p in supported]
    output={
      "version":"v15.5_real_Czech_2023_OSM_graph_35_author_station_centroids",
      "original_botanical_source_DOI":"10.23855/preslia.2024.247",
      "Geofabrik_original_historical_file":SOURCE_OSM,
      "Geofabrik_historical_snapshot_date":SOURCE_OSM_DATE,
      "original_OSM_archive_bytes":size,
      "actual_download_source_SHA256":sha,
      "original_plant_study_years":"2020–2022",
      "source_SR70_station_centroid_reference_year":"2026-08-15",
      "original_author_39_stations":39,
      "conservatively_location_matched_exact_station_names":len(sites),
      "rail_source_graph":{k:v for k,v in q["rail"].items() if k not in ("snaps","routes")},
      "road_source_graph":{k:v for k,v in q["road"].items() if k not in ("snaps","routes")},
      "n_original_35site_pairs":len(pairs),
      "n_original_source_rail_sections":len({s["section"] for s in sites}),
      "source_min_geodesic_distance_km":MIN_INTERSITE_GEO_KM,
      "n_pairs_with_both_real_2023_physical_routes_and_geo_ge_2km":len(supported),
      "n_same_stretch_supported_pairs":sum(p["source_same_stretch"] for p in supported),
      "n_other_stretch_supported_pairs":sum(not p["source_same_stretch"] for p in supported),
      "n_noneligible_source_pairs_NOT_assumed_zero":len(pairs)-len(supported),
      "rail_road_path_length_ratio_median":median(ratios) if ratios else None,
      "rail_road_path_length_ratio_p5_p95":(
          [float(np.percentile(ratios,5)),float(np.percentile(ratios,95))] if ratios else None),
      "proportion_rail_physically_shorter_than_road":sum(
          r["rail_graph_source_shortest_path_km"]<r["road_graph_source_shortest_path_km"]
          for r in supported)/len(supported) if supported else None,
      "research_limits":[
          "Original 2020-22 vegetation site positions unavailable: 2026 SR70 stations are labelled centroids, not individual source field plots.",
          "Physical historic 2023 OSM graph allows only possible connectivity, not confirmed freight trains, train speeds or real seed movement.",
          "National Czech country extract clips paths outside boundaries; physical disconnected pairs do not imply zero seed movement.",
          "The exact 35 matched source stations selected only by name, four others not inferred.",
          "Actual source original plant station areas differ; railway site transect sampling effort not standardized to Slovenia 100m site units.",
          "Road shortest path is not road traffic flux; neither graph is a randomized natural vegetation survey.",
      ],
      "NO_Czech_railway_dispersal_effect_estimated_from_topology_alone":True,
      "copyright_Geofabrik_OSM":"© OpenStreetMap contributors, ODbL; Geofabrik extract",
    }
    if return_pairs:return output,pairs
    return output


def run(s1,s2,disp,life,sr70,pbf,return_internal=False):
    sites,guild,source=original_35_station_sites(s1,s2,disp,life,sr70)
    output,pairs=graph_station_source(sites,pbf,return_pairs=True)
    output["original_author_plant_source_and_station_geographic_provenance"]=source
    if return_internal:return output,pairs,sites,guild
    return output


def main():
    p=argparse.ArgumentParser()
    for n in ("s1","s2","pladias-diaspore","pladias-lifeform","sr70","osm-pbf","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    result=run(a.s1,a.s2,a.pladias_diaspore,a.pladias_lifeform,a.sr70,a.osm_pbf)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in (
       "original_author_plant_source_and_station_geographic_provenance","research_limits")},ensure_ascii=False))


if __name__=="__main__":
    main()
