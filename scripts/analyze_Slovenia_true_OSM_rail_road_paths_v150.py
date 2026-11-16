#!/usr/bin/env python3
"""v15.0 author-original Slovenian railway field sites vs TRUE historic OSM rail/road paths.

Source A: Šabić et al 2026, DOI 10.1007/s44473-026-00155-z;
89 actual BS+BT matched vegetation sites with real 100-m WGS84 transect GPS.
Source B: Geofabrik archive: OSM Slovenia 2024-01-01 PBF, not 2022–23
operational timetables. Original OSM node IDs form only real graph edges.
Crossing rail/road geometries are NEVER spliced into shared nodes.
This is geometric topology ONLY, no reliable freight routing or seed transport.

Only source-original publisher XLSX (SHA validated by existing v10.2/v14.2)
and OSM PBF (full SHA/size captured) are used. No raw PBF or XLSX rehosted.
Source field 89-site GPS and graph node IDs are never published in output.
Cross-border routes may be truncated in the Slovenia-only extraction.
Rail=rail,narrow_gauge. Road=drivable road graph, no private accessibility/
speed validation. Study year and historical graph year differ explicitly.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path
from statistics import median

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx,SOURCE_DOI,SOURCE_SHA256
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    original_gps_coordinates,haversine_km,source_pairs)

GEOGRAPHIC_SOURCE="https://download.geofabrik.de/europe/slovenia-240101.osm.pbf"
SOURCE_DATE="2024-01-01"
RAIL_KINDS=frozenset(("rail","narrow_gauge"))
ROAD_KINDS=frozenset((
    "motorway","motorway_link","trunk","trunk_link","primary","primary_link",
    "secondary","secondary_link","tertiary","tertiary_link","unclassified",
    "residential","living_street","service",
))
MAIN_RAIL_MAX_SNAP_KM=0.50
MAIN_ROAD_MAX_SNAP_KM=2.00
MIN_INTERSITE_GEO_KM=2.0


class SourceNetwork:
    def __init__(self):
        self.coords={}
        self.edge_min={}
        self.eligible_ways=0
        self.skipped_bad_location=0
        self.skipped_bad_edge=0

    def append_way(self,osm_points):
        """Physical OSM way node references with actual lat/lon coordinates.
        Never join ways merely because their geometric segments cross."""
        if len(osm_points)<2:
            return
        self.eligible_ways+=1
        last=None
        for ref,lat,lon in osm_points:
            if not(-90<=lat<=90 and -180<=lon<=180 and math.isfinite(lat) and math.isfinite(lon)):
                self.skipped_bad_location+=1
                last=None
                continue
            self.coords.setdefault(int(ref),(float(lat),float(lon)))
            if last is not None and last[0]!=ref:
                d=haversine_km(last[1],last[2],lat,lon)
                if not(0<d<100):
                    self.skipped_bad_edge+=1
                else:
                    key=(min(last[0],ref),max(last[0],ref))
                    prior=self.edge_min.get(key)
                    if prior is None or d<prior:self.edge_min[key]=d
            last=(int(ref),float(lat),float(lon))

    def finalize(self):
        import numpy as np
        from scipy.sparse import coo_matrix
        if len(self.coords)<2 or len(self.edge_min)<1:
            raise ValueError("Original physical OSM network has no usable way topology")
        ids=list(self.coords)
        ix={osm_id:i for i,osm_id in enumerate(ids)}
        n=len(ids)
        m=len(self.edge_min)
        starts=np.fromiter((ix[u] for u,v in self.edge_min),dtype="int64",count=m)
        ends=np.fromiter((ix[v] for u,v in self.edge_min),dtype="int64",count=m)
        weights=np.fromiter(self.edge_min.values(),dtype="float64",count=m)
        if not np.all(weights>0):raise ValueError("Source physical OSM edge lengths invalid")
        # Distinct OSM-node pairs are guaranteed by edge_min dict.
        # Symmetrize physical line geometry for natural propagule movement
        # diagnostic; this does NOT establish active bi-directional rail traffic.
        graph=coo_matrix((
            np.concatenate((weights,weights)),
            (np.concatenate((starts,ends)),np.concatenate((ends,starts)))
        ),shape=(n,n)).tocsr()
        graph.sort_indices()
        coordinates=np.asarray([self.coords[k] for k in ids],dtype="float64")
        from scipy.sparse.csgraph import connected_components
        count,component=connected_components(graph,directed=False)
        biggest=int(max(Counter(component).values()))
        return graph,coordinates,component,{
            "original_osm_selected_way_count":self.eligible_ways,
            "original_osm_unique_node_count":n,
            "original_osm_unique_undirected_edges":m,
            "original_osm_connected_components":int(count),
            "original_osm_largest_connected_component_node_fraction":biggest/n,
            "source_invalid_node_positions_not_silently_joined":self.skipped_bad_location,
            "source_invalid_length_edges_skipped":self.skipped_bad_edge,
            "route_edges_undirected_for_geometric_opportunity_NOT_actual_oneway_service":True,
        }


def read_osm_network(pbf):
    try:import osmium
    except ImportError as exc:
        raise RuntimeError("Install pyosmium for genuine source PBF processing") from exc
    rail,road=SourceNetwork(),SourceNetwork()
    class OriginalWays(osmium.SimpleHandler):
        def way(self,w):
            mode=w.tags.get("railway")
            highway=w.tags.get("highway")
            target_rail=(mode in RAIL_KINDS and w.tags.get("construction") is None
                         and w.tags.get("disused") is None)
            target_road=highway in ROAD_KINDS
            if not target_rail and not target_road:return
            points=[]
            for n in w.nodes:
                if not n.location.valid():
                    points.append((int(n.ref),math.nan,math.nan))
                else:
                    points.append((int(n.ref),float(n.lat),float(n.lon)))
            if target_rail:rail.append_way(points)
            if target_road:road.append_way(points)
    OriginalWays().apply_file(str(pbf),locations=True,idx="sparse_mem_array")
    return rail,road


def snap_sites(coordinates,sites,maximum_km):
    """Find nearest genuine OSM node; never create artificial rail connections."""
    import numpy as np
    from scipy.spatial import cKDTree
    if coordinates.shape[1]!=2 or len(coordinates)<2:
        raise ValueError("Not valid physical railway/road OSM node coordinates")
    km_lon=111.320*math.cos(math.radians(46.))
    km_lat=111.132
    src=np.asarray([[s["lat"]*km_lat,s["lon"]*km_lon] for s in sites])
    net=np.column_stack((coordinates[:,0]*km_lat,coordinates[:,1]*km_lon))
    k=min(8,len(coordinates))
    _,candidate=cKDTree(net).query(src,k=k)
    if k==1:candidate=candidate.reshape((-1,1))
    result=[]
    for s,opts in zip(sites,candidate):
        distances=[
            (haversine_km(s["lat"],s["lon"],coordinates[int(i),0],coordinates[int(i),1]),
             int(i)) for i in np.atleast_1d(opts)
        ]
        dist,idx=min(distances)
        result.append({"node":idx if dist<=maximum_km else None,
                       "distance_km":dist,
                       "on_supported_source_network":bool(dist<=maximum_km)})
    return result


def pair_paths(graph,snaps):
    """True source-graph shortest route distances, no straight-line fill."""
    import numpy as np
    from scipy.sparse.csgraph import dijkstra
    n=len(snaps)
    out=np.full((n,n),np.inf,dtype="float64")
    eligible=[(i,x["node"]) for i,x in enumerate(snaps) if x["node"] is not None]
    if not eligible:return out
    unique_targets=sorted({v for i,v in eligible})
    for k in range(0,len(unique_targets),6):
        sources=unique_targets[k:k+6]
        lengths=dijkstra(graph,directed=False,indices=sources)
        lengths=np.atleast_2d(lengths)
        for row,source in enumerate(sources):
            for i,start in eligible:
                if start!=source:continue
                for j,dest in eligible:
                    out[i,j]=float(lengths[row,dest])
    return out


def run_network(original,pbf,return_network_matrices=False):
    import numpy as np
    true_sites,guild,meta=source_sites_and_guilds(original)
    sites,gps=original_gps_coordinates(original,true_sites)
    pairs=source_pairs(sites,guild)
    raw=Path(pbf).read_bytes()
    if len(raw)<1_000_000:
        raise ValueError("Truncated or non-genuine OSM extract")
    pbf_sha=hashlib.sha256(raw).hexdigest()
    # Free original PBF bytes before constructing physical node graphs.
    del raw
    rail,road=read_osm_network(pbf)
    graph_data={}
    for name,source,maxsnap in (
        ("rail",rail,MAIN_RAIL_MAX_SNAP_KM),
        ("road",road,MAIN_ROAD_MAX_SNAP_KM),
    ):
        graph,coords,component,graphqc=source.finalize()
        snaps=snap_sites(coords,sites,maxsnap)
        dist=pair_paths(graph,snaps)
        good_sites=[i for i,x in enumerate(snaps) if x["on_supported_source_network"]]
        graph_data[name]={
            "snaps":snaps,"distance":dist,
            "original_OSM_graph_QC":graphqc,
            "n_author_field_sites_snapped_within_threshold":len(good_sites),
            "n_author_field_sites_outside_threshold":len(sites)-len(good_sites),
            "snap_max_km_predeclared":maxsnap,
            "median_actual_snap_km":median(s["distance_km"] for s in snaps),
            "max_actual_snap_km":max(s["distance_km"] for s in snaps),
            "n_source_pairs_physically_connected_with_valid_snap":sum(
                math.isfinite(dist[i,j])
                for i,j in combinations(range(len(sites)),2)),
            "n_source_site_pairs_with_two_valid_snaps":sum(
                snaps[i]["on_supported_source_network"] and snaps[j]["on_supported_source_network"]
                for i,j in combinations(range(len(sites)),2)),
        }
        del source,graph,coords,component
    annotated=[]
    for p in pairs:
        i,j=p["i"],p["j"]
        x=graph_data["rail"]["distance"][i,j]
        y=graph_data["road"]["distance"][i,j]
        rail_ok=(np.isfinite(x) and x>0)
        road_ok=(np.isfinite(y) and y>0)
        valid=bool(rail_ok and road_ok and p["geo_km"]>=MIN_INTERSITE_GEO_KM)
        annotated.append({
            "i":i,"j":j,"same_original_section":p["same"],
            "geo_km":p["geo_km"],
            "geographic_bin":p["geo_bin"],
            "functional_bin":p["functional_bin"],
            "functional_overlap":p["fun"],
            "same_original_year":p["same_year"],
            "same_original_track_type":p["same_track_type"],
            "source_original_jaccard":p["observed_known_species_shared_jaccard"],
            "rail_path_km":float(x) if rail_ok else None,
            "road_path_km":float(y) if road_ok else None,
            "source_rail_both_site_access_km":sum(
                graph_data["rail"]["snaps"][z]["distance_km"] for z in (i,j)),
            "source_road_both_site_access_km":sum(
                graph_data["road"]["snaps"][z]["distance_km"] for z in (i,j)),
            "supported_BOTH_rail_road_routes_and_geo_ge_2km":valid,
            "log_rail_to_road_path_ratio":math.log(x/y) if valid else None,
            "observed_jaccard_only_not_physical_dispersal":True,
        })
    on=[x for x in annotated if x["supported_BOTH_rail_road_routes_and_geo_ge_2km"]]
    ratio=[x["rail_path_km"]/x["road_path_km"] for x in on]
    provenance={
        "version":"v15.0_GENUINE_historical_rail_road_graph_route_feasibility",
        "original_vegetation_DOI":SOURCE_DOI,
        "original_vegetation_source_SHA256":SOURCE_SHA256,
        "original_OSM_pbf_public_archive_URL":GEOGRAPHIC_SOURCE,
        "original_OSM_extract_snapshot_date":SOURCE_DATE,
        "original_OSM_pbf_SHA256":pbf_sha,
        "original_OSM_pbf_bytes":Path(pbf).stat().st_size,
        "source_data_OSM_copyright":"© OpenStreetMap contributors / ODbL, Geofabrik extract",
        "original_author_survey_years":"2022–2023",
        "original_field_sites":89,
        "original_rail_sections":15,
        "original_species_with_both_traits":meta[
            "original_strategy_lifespan_known_taxa_present_in_matched_sites"],
        "actual_source_field_GPS_100m_transitions":gps[
            "verified_original_geographical_WGS84_site_GPS_endpoint_records"],
        "rail_source":{k:v for k,v in graph_data["rail"].items() if k not in ("snaps","distance")},
        "road_source":{k:v for k,v in graph_data["road"].items() if k not in ("snaps","distance")},
        "n_original_field_site_pairs":len(pairs),
        "predeclared_min_inter_site_distance_km":MIN_INTERSITE_GEO_KM,
        "n_geo_and_both_network_paths_supported":len(on),
        "n_same_section_supported":sum(x["same_original_section"] for x in on),
        "n_between_section_supported":sum(not x["same_original_section"] for x in on),
        "n_unavailable_pairs_not_imputed_as_zero":len(pairs)-len(on),
        "fraction_supported_rail_shorter_than_road":sum(x["rail_path_km"]<x["road_path_km"] for x in on)/len(on) if on else None,
        "rail_over_road_path_ratio_median":median(ratio) if ratio else None,
        "rail_over_road_path_ratio_5_95_pct":(
            [float(np.percentile(ratio,5)),float(np.percentile(ratio,95))] if ratio else None),
        "all_site_vegetation_rail_no_nonrail_controls":True,
        "not_active_freight_road_traffic_or_directional_train_services":True,
        "geo_and_road_rail_correlation_may_not_identify_additional_rail_effect":True,
        "OSM_Slovenia_country_border_TRUNCATES_outside_international_routes":True,
        "not_proof_of_seeds_carried_by_train":True,
        "not_mowing_wind_stress_experimental_test":True,
    }
    if return_network_matrices:
        return provenance,annotated,sites,guild
    return provenance


def main():
    parser=argparse.ArgumentParser()
    for option in ("original-xlsx","source-audit102","osm-pbf","output"):
        parser.add_argument("--"+option,required=True)
    a=parser.parse_args()
    source=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    result=run_network(source,a.osm_pbf)
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
