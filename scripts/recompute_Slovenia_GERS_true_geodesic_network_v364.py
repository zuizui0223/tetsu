#!/usr/bin/env python3
"""v364: correct Overture GERS source 'at' as normalized WGS84 geodesic
distance, not local UTM length fraction. Revalidate genuine Slovenian
89-trackside-location source Jaccard vs actual GERS physical path length.

This source-unit audit does NOT change original botanical study source,
OSM railway names, actual connector UUID connectivity or sample identities.
"""
from __future__ import annotations
import argparse,json,math,statistics
from collections import Counter
from pathlib import Path

import analyze_Slovenia_Overture_connector_network_botany_v358 as v358
import analyze_Slovenia_geod_and_rail_distance_matched_v359 as v359
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites
from audit_Slovenia_source_same_line_vs_other_geodistance_matched_v280 import (
    match_candidates,line_equal_aggregate
)
from audit_Slovenia_same_line_geodistance_matched_frequency_null_v281 import (
    original_matched_gap_for_rows,source_groups,experiment
)

SOURCE_RELEASE=v358.RELEASE
SOURCE_OLD_GRAPH_ORIGINAL_SITES=89
SOURCE_OLD_GRAPH_ORIGINAL_SECTIONS=15
SOURCE_OLD_GRAPH_MATCHED_86_AMBIGUOUS_3=(86,3)
SOURCE_V359_OLD_N=27
SOURCE_V359_OLD_SECTIONS=7
SOURCE_V359_OLD_GAP=0.027033311999075067
ANALYSIS_NULL_KINDS=("all_original_sites","within_original_survey_year_and_track_layout")
ANALYSIS_SEED=20261124+359_000

def source_geodetic_line_m(source_lonlat,geod):
    """True WGS84 geodetic polyline length, in metres, matching GERS 'at'."""
    coords=list(source_lonlat.coords)
    if len(coords)<2:raise ValueError("Original Overture rail LineString too short")
    lengths=[]
    for (x0,y0),(x1,y1) in zip(coords,coords[1:]):
        segment=abs(float(geod.inv(x0,y0,x1,y1)[2]))
        if not math.isfinite(segment) or segment<0:
            raise ValueError("Original WGS84 GERS line segment malformed")
        lengths.append(segment)
    total=sum(lengths)
    if total<0.1:raise ValueError("Original global source rail segment geodetic length insufficient")
    return total,lengths

def snapped_WGS84_source_at(lonlat,site_xy,transformer,geod):
    """Closest projected original vertex segment, then WGS84 geodesic cumulative fraction.

    The short local UTM projection determines the nearest original segment
    and its within-vertex position; length along the original rail always
    uses the true source WGS84 ellipsoid, not UTM line-length fraction.
    """
    coords=list(lonlat.coords)
    L,lengths=source_geodetic_line_m(lonlat,geod)
    verts=[transformer.transform(x,y) for x,y in coords]
    cx,cy=site_xy
    candidates=[];before=0.
    for i,((x0,y0),(x1,y1),lm) in enumerate(zip(verts,verts[1:],lengths)):
        vx=x1-x0;vy=y1-y0
        denom=vx*vx+vy*vy
        t=min(1.,max(0.,((cx-x0)*vx+(cy-y0)*vy)/denom)) if denom>0 else 0.
        d=(cx-(x0+t*vx))**2+(cy-(y0+t*vy))**2
        candidates.append((d,i,t,before+lm*t))
        before+=lm
    _,_,_,distance=min(candidates)
    original_at=min(1.,max(0.,distance/L))
    return original_at,math.sqrt(min(x[0] for x in candidates))

def corrected_physical_graph(source_rows,transformer):
    """Recreate baseline exact UUID topology, replace only source edge metres."""
    from shapely import wkb
    from pyproj import Geod
    import networkx as nx
    legacy,segments,qa=v358.build_rail_graph(source_rows,transformer)
    original_by_id={str(r.get("id")):r for r in source_rows if v358.source_eligible(r)}
    geod=Geod(ellps="WGS84")
    graph=nx.Graph()
    ratios=[]
    for s in segments:
        raw=original_by_id.get(s["id"])
        if raw is None:raise ValueError("Source GERS segment lacks immutable raw original")
        rawline=wkb.loads(bytes(raw["geometry"]))
        lm,per_edge=source_geodetic_line_m(rawline,geod)
        s["source_lonlat_line_PRIVATE"]=rawline
        s["source_geodetic_length_m"]=lm
        s["source_per_original_vertex_geodetic_m_PRIVATE"]=per_edge
        ppm=(s["original_projected_length_m"]/lm-1.)*1e6
        ratios.append(abs(ppm))
        v358.add_segment_edges(graph,s["original_connector_chain"],lm)
    if (graph.number_of_nodes()!=legacy.number_of_nodes()
        or graph.number_of_edges()!=legacy.number_of_edges()
        or set(graph.nodes)!=set(legacy.nodes)):
        raise ValueError("Geodetic correction changed original source graph physical UUID connectivity")
    return legacy,graph,segments,qa,{
       "source_segment_count":len(segments),
       "source_projected_vs_WGS84_geodetic_abs_ppm_median":statistics.median(ratios),
       "source_projected_vs_WGS84_geodetic_abs_ppm_max":max(ratios),
       "source_legacy_graph_edges":legacy.number_of_edges(),
       "source_corrected_graph_edges":graph.number_of_edges(),
       "source_exact_same_original_UUID_nodes":True}

def corrected_source_attachments(corrected_graph,source_segments,source_sites,transformer,threshold=v358.MAX_DISTANCE_FROM_TRACK_METRES):
    """Keep legacy species-independent nearest rail assignments, fix site->GERS edge lengths."""
    from shapely.geometry import Point
    from pyproj import Geod
    graph,out=v358.site_attach(corrected_graph,source_segments,source_sites,transformer,threshold=threshold)
    byid={s["id"]:s for s in source_segments}
    geod=Geod(ellps="WGS84")
    fracdelta=[]
    for i,a in enumerate(out):
        if a["status"]!="matched_to_real_rail_connector_network":continue
        seg=byid[a["source_matched_segment_id_private"]]
        xy=transformer.transform(source_sites[i]["lon"],source_sites[i]["lat"])
        f,offset=snapped_WGS84_source_at(seg["source_lonlat_line_PRIVATE"],xy,transformer,geod)
        chain=seg["original_connector_chain"]
        if f < chain[0][0]-1e-10 or f > chain[-1][0]+1e-10:
            raise ValueError("Original botanist site falls outside source GERS connector chain; reject extrapolation")
        k=0
        while k<len(chain)-2 and f>chain[k+1][0]:k+=1
        (start,c0),(end,c1)=chain[k],chain[k+1]
        if not start-1e-10<=f<=end+1e-10:
            raise ValueError("Correct source at not enclosed by GERS connector edge")
        f=min(end,max(start,f))
        # Replace legacy projected-fraction site connections; original node and
        # GERS UUID keys do not change. Other source sites remain intact.
        node=a["source_graph_node"]
        graph.remove_node(node)
        L=seg["source_geodetic_length_m"]
        graph.add_edge(node,c0,weight=max(0.,L*(f-start)))
        graph.add_edge(node,c1,weight=max(0.,L*(end-f)))
        fracdelta.append(abs(f-a["source_fraction"]))
        a["source_old_projected_fraction_PRIVATE"]=a["source_fraction"]
        a["source_fraction"]=f
        a["source_original_length_m"]=L
        a["source_correct_geodetic_linear_referencing"]=True
    return graph,out,{
        "source_original_site_matches_recalibrated":len(fracdelta),
        "original_botanical_GPS_to_track_distance_max_threshold_m":threshold,
        "source_site_source_at_fraction_delta_median":statistics.median(fracdelta) if fracdelta else None,
        "source_site_source_at_fraction_delta_max":max(fracdelta) if fracdelta else None,
        "source_corrected_site_connectors_attached_by_true_WGS84_fraction":True}

def compare_original_paths(old_g,old_a,new_g,new_a):
    old=v359.physical_network_distances(old_g,old_a)
    new=v359.physical_network_distances(new_g,new_a)
    both=set(old)&set(new)
    errors_m=[abs(new[key]-old[key])*1000 for key in both]
    diffs=[new[key]-old[key] for key in both]
    return old,new,{
      "original_directed_genuine_source_GERS_paths_legacy":len(old),
      "corrected_directed_genuine_source_GERS_paths":len(new),
      "original_same_connected_pair_keys_n":len(both),
      "paths_newly_lost_or_gained_at_fixed_300km_cutoff_n":len(set(old)^set(new)),
      "source_WGS84_network_length_change_abs_m_median":statistics.median(errors_m) if errors_m else None,
      "source_WGS84_network_length_change_abs_m_max":max(errors_m) if errors_m else None,
      "source_WGS84_network_length_change_signed_median_m":statistics.median(diffs)*1000 if diffs else None}

def null_for_corrected_matched_panel(sites,matched,source_effect):
    if not matched or len(matched)<12 or len({x["source_author_section"] for x in matched})<4:
        return {"identified":False,"not_source_no_effect":True}
    pairs=[{"anchor_index":x["anchor"],"within_index":x["same"],"other_index":x["different"]}
           for x in matched]
    initial=[set(s["BS"]|s["BT"]) for s in sites]
    obs,_=original_matched_gap_for_rows(initial,sites,pairs)
    if abs(obs-source_effect)>1e-10:
        raise ValueError("Corrected original source plot Jaccard score differs from matched original null evaluation")
    outcomes={}
    for k,kind in enumerate(ANALYSIS_NULL_KINDS):
        outcomes[kind]=experiment(sites,pairs,source_groups(sites,kind),
             seed=ANALYSIS_SEED+99991*k,draws=99,burnin=3000,thin=350)
    return {"identified":True,"source_observed_original_Jaccard":obs,
            "original_source_fixed_species_prevalence_technical_null":outcomes,
            "technical_upper_rank_NOT_field_experimental_p":True}

def run(original,rail_factory=v358.real_overture):
    from pyproj import Transformer
    from networkx.algorithms.components import number_connected_components
    bysec,originalQA=to_original_sites(original)
    sites=[x for key in sorted(bysec) for x in bysec[key]]
    if len(sites)!=SOURCE_OLD_GRAPH_ORIGINAL_SITES or len(bysec)!=SOURCE_OLD_GRAPH_ORIGINAL_SECTIONS:
        raise ValueError("Authentic Slovenian original 89 sites /15 rail sections changed")
    source59,source_match_qa=match_candidates(sites,40.0,2.5,"geo_only")
    if len(source59)!=59:
        raise ValueError("Original botanical v28 geo-only matched 59 anchors source drift")
    meta=line_equal_aggregate(source59,source_match_qa)
    observed59=meta["habitat_comparisons"]["BS_and_BT"]["equal_source_line_same_minus_different_line_similarity"]
    if abs(observed59-v358.SOURCE_OLD_GEO_MATCHED_EFFECT)>1e-10:
        raise ValueError("Original source 59-anchor baseline botanist identities drift")
    box=(min(x["lon"] for x in sites)-.12,min(x["lat"] for x in sites)-.12,
         max(x["lon"] for x in sites)+.12,max(x["lat"] for x in sites)+.12)
    source_rows,sourceQA=rail_factory(box)
    tr=Transformer.from_crs("EPSG:4326",v358.SOURCE_UTM,always_xy=True)
    legacy,correct,segments,graphqa,deltaQA=corrected_physical_graph(source_rows,tr)
    if len(segments)!=6612 or correct.number_of_nodes()!=12346 or correct.number_of_edges()!=14181:
        raise ValueError("Original source 2026-11 regional railway GERS graph changed")
    old_aug,old_site=v358.site_attach(legacy,segments,sites,tr)
    new_aug,new_site,attachQA=corrected_source_attachments(correct,segments,sites,tr)
    for loc in (old_site,new_site):
        nmatch=sum(x["status"]=="matched_to_real_rail_connector_network" for x in loc)
        namb=sum(x["status"]=="ambiguous_multiple_disconnected_source_tracks" for x in loc)
        if (nmatch,namb)!=SOURCE_OLD_GRAPH_MATCHED_86_AMBIGUOUS_3:
            raise ValueError("Original site match changed by WGS84 source distance correction")
    oldpath,newpath,pathqa=compare_original_paths(old_aug,old_site,new_aug,new_site)
    origmatch,_=v359.joint_geo_and_network_triples(sites,oldpath,"geo_and_network")
    if len(origmatch)!=SOURCE_V359_OLD_N or len({x["source_author_section"] for x in origmatch})!=SOURCE_V359_OLD_SECTIONS:
        raise ValueError("Original published v359 27 anchors/7 sections source unchanged check failed")
    oldscore=v359.botanical_comparison(sites,origmatch)["source_habitat_comparisons"][
        "BS_and_BT"]["source_author_section_equal_same_minus_other_Jaccard"]
    if abs(oldscore-SOURCE_V359_OLD_GAP)>1e-10:
        raise ValueError("Original published v359 Jaccard +0.027033 remeasurement failed")
    newmatch,newmatchqa=v359.joint_geo_and_network_triples(sites,newpath,"geo_and_network")
    strictmatch,strictqa=v359.joint_geo_and_network_triples(sites,newpath,"geo_network_year_layout")
    outcome=v359.botanical_comparison(sites,newmatch)
    strictoutcome=v359.botanical_comparison(sites,strictmatch)
    frozen59corrected=v358.summarize_matched_botany(sites,source59,new_site,new_aug,120.)
    oldref={(x["anchor"],x["same"],x["different"]) for x in origmatch}
    newref={(x["anchor"],x["same"],x["different"]) for x in newmatch}
    original_other_comparison=outcome.get("source_habitat_comparisons",{}).get(
         "BS_and_BT",{}).get("source_author_section_equal_same_minus_other_Jaccard")
    freqnull=null_for_corrected_matched_panel(sites,newmatch,original_other_comparison)
    return {"version":"v364-source-true-WGS84-geodesic-rail-graph-corrected-Slovenia-flora",
      "source_original_botanical_DOI":"10.1007/s44473-026-00155-z",
      "source_original_botanical_SHA256":v358.SOURCE_ORIGINAL_SHA,
      "source_Overture_original_release":SOURCE_RELEASE,
      "source_real_region_original_rail_QA":sourceQA,
      "source_original_survey89_sites_section15":True,
      "source_real_6612_segment_12346_connector_14181_edge_topology_preserved":True,
      "source_geodetic_vs_legacy_metric_length_comparison":deltaQA,
      "source_site_attachment_true_geodetic_linear_reference_QA":attachQA,
      "source_original_true_physical_GERS_shortest_path_legacy_vs_geodesic_QA":pathqa,
      "source_first59_unchanged_geographic_plant_triples_remeasured":frozen59corrected,
      "source_v359_original_published_geo_network_Jaccard_gap":oldscore,
      "source_v359_original_published_matched_source_anchors":len(origmatch),
      "source_v364_new_true_WGS84_matched_source_anchors":len(newmatch),
      "source_original_matched_triads_unchanged_n":len(oldref&newref),
      "source_original_matched_triads_displaced_by_correct_source_metres_n":len(oldref-newref),
      "source_v364_geo_and_WGS84_rail_network_match_QA":newmatchqa,
      "source_v364_source_corrected_network_distance_matched_plant_Jaccard":outcome,
      "source_v364_strict_year_layout_source_support":strictqa,
      "source_v364_strict_year_layout_Jaccard_and_support":strictoutcome,
      "source_v364_original_species_frequency_preserving_technical_null":freqnull,
      "source_botanist_site_original_gps_lineIDs_taxon_occurrences_GERS_geometries_not_republished":True,
      "no_causal_train_seed_or_global_plant_network_convergence_claim":True}

def main():
    p=argparse.ArgumentParser()
    for key in ("original-xlsx","source-audit","output"):p.add_argument("--"+key,required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.original_xlsx,a.source_audit)
    r=run(original)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(r,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
