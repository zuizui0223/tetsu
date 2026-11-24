#!/usr/bin/env python3
"""v366a: genuinely independent Czech35 and Slovak19 author railway station
floras mapped to exact Overture 2026 physical GERS connector network.

SOURCE FEASIBILITY ONLY: station-to-source-rail matching and actual original
connected path support; NO plant similarity/network biological effect fit.
Whole Czech botanical stations use 2026 national SR70 representative station
centroids (NOT source plot GPS). Slovak paper prints original 2012-2014
station DMS, not in-track transect coordinates. Source missingness != biology.
"""
from __future__ import annotations
import argparse,hashlib,json,math,time
from collections import Counter
from itertools import combinations
from pathlib import Path

from audit_Czech_2024_true_35_stations_stretch_geodistance_matched_v282 import (
    source_original_czech_station_floras
)
from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import (
    source_list_clauses,parse_clauses,validate_original_source,station_sets
)
from analyze_Slovakia2016_19station_original_physical_distance_v310 import (
    parse_source_dms
)
from analyze_Slovenia_geod_and_rail_distance_matched_v359 import (
    physical_network_distances
)
from recompute_Slovenia_GERS_true_geodesic_network_v364 import (
    corrected_physical_graph,corrected_source_attachments
)

RELEASE="2026-11-18.0"
S1_SHA="df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06"
S2_SHA="00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2"
SR70_SHA="fc43a1969af64f17285c7c6e3a30b00c3994d0ba1cfda6db09d8e4baaa166019"
SK_SHA="99465bdd7ed56985e65b1163f16eca4f702fc821619d84c89bf986d4d3736ae1"
BBOX_MARGIN_DEGREES=.18
MAX_ORIGINAL_SOURCE_TRANSPORT=45_000_000
MAX_RAIL_ROWS=200_000
SOURCE_STUDIES={
 "CZ":{"n_source_botanical_sites":35,"n_original_author_rail_stretches":5,
       "location_type":"2026 THIRD-PARTY-DERIVED national SR70 station centroid, NOT 2020-22 botanical plot GPS",
       "source_actual_rail_near_max_m":300.0,
       "source_strict_accuracy_sensitivity_m":120.0,
       "projected_WGS84_to_metric_UTM":"EPSG:32633",
       "source_survey_years":"2020-2022",
       "source_original_taxa":761},
 "SK":{"n_source_botanical_sites":19,"n_original_author_rail_stretches":1,
       "location_type":"AUTHOR ORIGINAL 2016 printed 19-station DMS coordinates of 2012-14 whole station surveys",
       "source_actual_rail_near_max_m":180.0,
       "source_strict_accuracy_sensitivity_m":80.0,
       "projected_WGS84_to_metric_UTM":"EPSG:32634",
       "source_survey_years":"2012-2014",
       "source_original_taxa":318}
}

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def original_CZ_sites(s1,s2,sr70):
    for p,sha,size in ((s1,S1_SHA,79261),(s2,S2_SHA,13533),(sr70,SR70_SHA,None)):
        if not Path(p).is_file() or digest(p)!=sha:
            raise ValueError("Original Czech publisher or source SR70 station centroid SHA drift")
        if size is not None and Path(p).stat().st_size!=size:
            raise ValueError("Original Czech botanical XLSX source bytes differ")
    src,qa=source_original_czech_station_floras(Path(s1),Path(s2),Path(sr70))
    n_taxa=len(set().union(*(x["taxa"] for x in src)))
    if len(src)!=35 or len({x["source_stretch"] for x in src})!=5 or n_taxa<600:
        raise ValueError("Original Czech source expected true 35 station FLORAS and 5 author railway stretches")
    return [{"lat":x["lat"],"lon":x["lon"],"source_original_group":x["source_stretch"],
             "original_taxa_count":len(x["taxa"])} for x in src],{
               "original_source_39_real_pub_stations_n":39,
               "source_35_geolocated_source_station_centroids_n":35,
               "original_source_4_unlocated_station_GPS_missing_NOT_imputed":4,
               "original_source_stretch_labels":5,
               "original_matched_station_unique_taxa_min":n_taxa,
               "original_publisher_S1_SHA256":S1_SHA,
               "original_publisher_corrected_S2_SHA256":S2_SHA,
               "original_station_registry_derived_SR70_SHA256":SR70_SHA}

def original_SK_sites(pdf):
    import fitz
    original=Path(pdf)
    if not original.is_file() or original.stat().st_size!=249373 or digest(original)!=SK_SHA:
        raise ValueError("Original Slovakia 2016 author paper source not SHA+size matched")
    book=fitz.open(str(original))
    if len(book)!=9:raise ValueError("Original Slovak author paper not nine original pages")
    txt="\n".join(book[i].get_text("text",sort=True) for i in (1,2))
    book.close()
    loc=parse_source_dms(txt)
    clauses,qa=source_list_clauses(original)
    obs,namesqa=parse_clauses(clauses)
    bqa=validate_original_source(obs)
    sets=station_sets(obs)
    if len(loc)!=19 or len(sets)!=19 or bqa["n_original_station_taxon_presence"]!=1763:
        raise ValueError("Real original 19 Slovak station DMS or 1763 occurrence cells missing")
    src=[{"lat":loc[i][0],"lon":loc[i][1],
          "source_original_group":"source_SINGLE_original_line",
          "original_taxa_count":len(sets[i])} for i in range(1,20)]
    return src,{"original_author_PDF_source_sha256":SK_SHA,
        "original_original_station_DMS_coordinates_count":19,
        "original_source_distinct_publisher_original_taxa":318,
        "original_source_station_species_positive_cells":1763,
        "original_author_rail_routes_n":1,
        "source_original_five_printed_richness_benchmarks_validated":True}

def original_overture_rail_source(bbox):
    import pyarrow.compute as pc
    from overturemaps import record_batch_reader
    errors=[]
    for attempt in (1,2,3):
        try:
            source=Counter(); rows=[]
            reader=record_batch_reader("segment",bbox=bbox,release=RELEASE,stac=True)
            if reader is None:raise ValueError("No authentic Overture source scan reader")
            for batch in reader:
                source["original_transport_segment_rows_scanned"]+=batch.num_rows
                if source["original_transport_segment_rows_scanned"]>MAX_ORIGINAL_SOURCE_TRANSPORT:
                    raise ValueError("Original country-bbox source transportation row cap; source unavailable")
                idx=batch.schema.get_field_index("subtype")
                if idx<0:raise ValueError("Original Overture source segment subtype field missing")
                selected=batch.filter(pc.equal(batch.column(idx),"rail"))
                source["source_rail_subtype_rows_in_CZ_SK_bbox"]+=selected.num_rows
                if selected.num_rows:
                    rows.extend(selected.to_pylist())
                if len(rows)>MAX_RAIL_ROWS:raise ValueError("Country bbox original railway source rail rows over bound")
            if not rows:raise ValueError("Original country railway source returned no original rail geometry")
            return rows,{
                "original_Overture_real_release":RELEASE,
                "original_Overture_source_complete":True,
                "original_source_n_all_transport_features":source["original_transport_segment_rows_scanned"],
                "original_source_n_all_rail_subtype":source["source_rail_subtype_rows_in_CZ_SK_bbox"],
                "source_data_latest_snapshot_2026_vs_original_field_date_mismatch":True,
                "source_bbox_WGS84_not_published_as_botanist_GPS":True,
                "source_attempt":attempt}
        except Exception as e:
            errors.append({"attempt":attempt,"reason":type(e).__name__,
                          "short_error":str(e)[:190]})
            if attempt<3:time.sleep(12 if attempt==1 else 35)
    raise RuntimeError("Authentic original Overture spatial source unavailable, NOT no railway: "+str(errors))

def coordinate_source_window(sites):
    lon=[x["lon"] for x in sites]
    lat=[x["lat"] for x in sites]
    bbox=(min(lon)-BBOX_MARGIN_DEGREES,min(lat)-BBOX_MARGIN_DEGREES,
          max(lon)+BBOX_MARGIN_DEGREES,max(lat)+BBOX_MARGIN_DEGREES)
    if not 12<=bbox[0]<=24 or not bbox[0]<bbox[2]<=25 or not 47<=bbox[1]<=51 or not 49<=bbox[3]<=52:
        raise ValueError("Independent Czech/Slovak paper coordinates outside source geography")
    return bbox

def source_rail_gps_coverage(sites,located,network):
    statuses=Counter(x["status"] for x in located)
    offsets=[x["original_source_track_offset_m"] for x in located
             if x.get("original_source_track_offset_m") is not None]
    n=len(sites)
    same=other=0
    reachable=0
    for i,j in combinations(range(n),2):
        if (i,j) not in network:continue
        reachable+=1
        if sites[i]["source_original_group"]==sites[j]["source_original_group"]:same+=1
        else:other+=1
    return {"source_original_true_physical_rail_matched_station_n":statuses[
             "matched_to_real_rail_connector_network"],
            "source_original_sites_ambiguous_disconnected_tracks_n":statuses[
             "ambiguous_multiple_disconnected_source_tracks"],
            "source_original_sites_without_a_qualifying_track_n":statuses[
             "unmatched_to_global_rail_source"],
            "source_botanist_original_rail_station_to_GERS_offset_m_median":
             sorted(offsets)[len(offsets)//2] if offsets else None,
            "source_botanist_original_rail_station_to_GERS_offset_m_max":
             max(offsets) if offsets else None,
            "source_potential_original_site_pairs":n*(n-1)//2,
            "source_real_network_connected_original_station_pairs_within_300km":reachable,
            "source_real_network_connected_pairs_SAME_original_author_stretch":same,
            "source_real_network_connected_pairs_OTHER_original_author_stretch":other,
            "source_no_rail_network_path_for_some_pairs_is_missing_NOT_negative_plant_ecology":True}

def run(study,args,network_source_factory=original_overture_rail_source):
    from pyproj import Transformer
    from analyze_Slovenia_Overture_connector_network_botany_v358 import source_eligible
    from analyze_Slovenia_geod_and_rail_distance_matched_v359 import physical_network_distances
    if study=="CZ":src,srcqa=original_CZ_sites(args.s1,args.s2,args.sr70)
    elif study=="SK":src,srcqa=original_SK_sites(args.pdf)
    else:raise ValueError("Unexpected original independent railway plant study")
    config=SOURCE_STUDIES[study]
    if len(src)!=config["n_source_botanical_sites"] or len({x["source_original_group"] for x in src})!=config["n_original_author_rail_stretches"]:
        raise ValueError("Independent original study source field stations or railway source stretches drift")
    bbox=coordinate_source_window(src)
    rows,sourceQA=network_source_factory(bbox)
    tran=Transformer.from_crs("EPSG:4326",config["projected_WGS84_to_metric_UTM"],always_xy=True)
    old,graph,segments,topqa,metricqa=corrected_physical_graph(rows,tran)
    if graph.number_of_edges()<30 or len(segments)<40:
        raise ValueError("Original source railway graph has insufficient true GERS geometric topology")
    outgraph,attached,siteqa=corrected_source_attachments(
       graph,segments,src,tran,threshold=config["source_actual_rail_near_max_m"])
    shorter,strict,strict_qa=corrected_source_attachments(
       graph,segments,src,tran,threshold=config["source_strict_accuracy_sensitivity_m"])
    physical=physical_network_distances(outgraph,attached)
    phystrict=physical_network_distances(shorter,strict)
    metrics=source_rail_gps_coverage(src,attached,physical)
    metrics_strict=source_rail_gps_coverage(src,strict,phystrict)
    return {"version":"v366-independent-Czech-Slovak-original-flora-physical-railway-source-feasibility",
       "source_study_country":study,"source_publisher_author_years":config["source_survey_years"],
       "source_original_station_coordinate_class_NOT_author_plot_GPS":config["location_type"],
       "source_original_botanical_station_count":len(src),
       "source_original_source_taxon_count_prior_publication":config["source_original_taxa"],
       "source_original_author_rail_stretch_or_line_count":config["n_original_author_rail_stretches"],
       "original_botanical_and_GPS_provenance":srcqa,
       "original_Overture_rail_source_data_QA":sourceQA,
       "original_Overture_rail_GERS_geodetic_component_QA":{
           "source_original_normal_rail_GERS_segments":len(segments),
           "source_original_connector_graph_nodes":graph.number_of_nodes(),
           "source_original_connector_graph_edges":graph.number_of_edges(),
           **metricqa},
       "source_original_primary_point_to_rail_m":config["source_actual_rail_near_max_m"],
       "source_original_strict_point_to_rail_m":config["source_strict_accuracy_sensitivity_m"],
       "source_primary_original_site_to_rail_and_connected_pairs_QA":metrics,
       "source_strict_source_site_to_rail_and_connected_pairs_QA":metrics_strict,
       "botanical_species_compositional_network_effect_NOT_YET_CALCULATED":True,
       "single_Slovak_original_publisher_railway_line_NOT_multiple_line_clusters":study=="SK",
       "CZ_centroid_is_not_botanist_actual_transect_GPS":study=="CZ",
       "source_2026_network_is_temporally_offset_from_field_surveys":True,
       "no_source_GPS_raw_botanical_taxa_original_PDF_or_original_rail_geometry_republished":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--study",choices=("CZ","SK"),required=True)
    p.add_argument("--s1");p.add_argument("--s2");p.add_argument("--sr70")
    p.add_argument("--pdf");p.add_argument("--output",required=True)
    a=p.parse_args()
    obj=run(a.study,a)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(obj,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
