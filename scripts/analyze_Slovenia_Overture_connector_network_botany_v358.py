#!/usr/bin/env python3
"""v358: genuine Slovenian rail flora (89 field sites / 15 original rail sections)
matched onto REAL 2026 Overture connector network rather than only author
source same-line label.

Exploratory source topology *contrast* with fixed original 2026 published
Slovenian source community and v28.0 already-geography-matched triples.
No original source plant/GPS/actual rail geometry published in GitHub.
"""
from __future__ import annotations
import argparse,json,math,statistics,time
from collections import Counter,defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites
from audit_Slovenia_source_same_line_vs_other_geodistance_matched_v280 import (
    match_candidates,line_equal_aggregate
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km

RELEASE="2026-11-18.0"
SOURCE_ORIGINAL_SHA="26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f"
SOURCE_N_SITES=89
SOURCE_N_SECTIONS=15
SOURCE_OLD_GEO_MATCHED_EFFECT=0.01967094157833441
SOURCE_OLD_ANCHORS=59
SOURCE_OLD_SOURCE_LINES=8
MAX_SOURCE_TRANSPORT_ROWS=8_000_000
ALLOWED_CLASSES={"standard_gauge","narrow_gauge","unknown"}
DROP_SOURCE_FLAGS={"is_disused","is_abandoned","is_under_construction"}
MAX_DISTANCE_FROM_TRACK_METRES=120.0
SENSITIVITY_DISTANCE_M=60.0
COMPETING_CANDIDATE_M=20.0
SOURCE_UTM="EPSG:32633"
MAX_GRAPH_ROUTE_KM=300.0

def source_eligible(row):
    if row.get("subtype")!="rail":return False
    if row.get("class") not in ALLOWED_CLASSES:return False
    for f in row.get("rail_flags") or []:
        if isinstance(f,dict) and set(f.get("values") or ()) & DROP_SOURCE_FLAGS:
            return False
    return True

def adjacent_connector_chain(row):
    """Source segment topology edges. Crossings without shared UUID never connect."""
    conns=row.get("connectors") or []
    seen=set(); chain=[]
    for item in conns:
        cid=str(item.get("connector_id") or "")
        val=item.get("at")
        try:at=float(val)
        except (TypeError,ValueError):continue
        if not cid or not math.isfinite(at) or at<0 or at>1:
            continue
        if (cid,round(at,10)) in seen:continue
        seen.add((cid,round(at,10)))
        chain.append((at,cid))
    chain.sort()
    if len(chain)<2:
        return []
    return chain

def add_segment_edges(graph,conns,segment_m):
    if not math.isfinite(segment_m) or segment_m<=0:return 0
    added=0
    for (at0,c0),(at1,c1) in zip(conns,conns[1:]):
        delta=at1-at0
        if c0==c1 or delta<=0:continue
        m=delta*segment_m
        current=graph.get_edge_data(c0,c1)
        if current is None or m<float(current["weight"]):
            graph.add_edge(c0,c1,weight=m)
            added+=1
    return added

def build_rail_graph(rows,transformer):
    import networkx as nx
    from shapely import wkb
    from shapely.ops import transform
    graph=nx.Graph()
    segments=[]
    counts=Counter()
    names=set()
    for row in rows:
        if row.get("subtype")=="rail":counts["all_rail_subtype_rows"]+=1
        if not source_eligible(row):continue
        geom=wkb.loads(bytes(row["geometry"]))
        if geom.geom_type!="LineString" or geom.is_empty:continue
        projected=transform(transformer.transform,geom)
        if projected.length<0.1:continue
        chain=adjacent_connector_chain(row)
        if len(chain)<2:
            counts["source_eligible_rail_rows_missing_usable_connectors"]+=1
            continue
        segid=str(row.get("id") or "")
        if not segid or segid in names:
            raise ValueError("Source Overture original GERS segment ID missing or duplicated")
        names.add(segid)
        added=add_segment_edges(graph,chain,projected.length)
        if added<=0:
            counts["source_original_eligible_rail_with_zero_added_edges"]+=1
        segments.append({"id":segid,"geom":projected,"original_connector_chain":chain,
                         "original_projected_length_m":projected.length})
        counts["source_usable_rail_segments"]+=1
        counts["source_rail_connector_incidences"]+=len(chain)
        counts["source_rail_names_present"]+=bool((row.get("names") or {}).get("primary"))
    return graph,segments,dict(counts)

def site_attach(graph,segments,sites,transformer,threshold=MAX_DISTANCE_FROM_TRACK_METRES):
    import networkx as nx
    from shapely.geometry import Point
    from shapely.strtree import STRtree
    from networkx.algorithms.components import node_connected_component
    shapes=[s["geom"] for s in segments]
    if not shapes:raise ValueError("No authentic Overture rail line geometries for botanical source")
    tree=STRtree(shapes)
    attached=[]
    outgraph=graph.copy()
    base_components={}
    for node in graph:
        if node in base_components:continue
        comp_id=len(base_components) # ordinal only, no reliance on ID meaning
        # use one representative that is fixed across all connectors
        comp=node_connected_component(graph,node)
        stable=min(comp)
        for n in comp:base_components[n]=stable
    for i,site in enumerate(sites):
        pt=Point(*transformer.transform(site["lon"],site["lat"]))
        raw_candidates=tree.query(pt.buffer(threshold))
        ds=sorted((shapes[int(j)].distance(pt),segments[int(j)]["id"],int(j)) for j in raw_candidates)
        if not ds or ds[0][0]>threshold:
            attached.append({"source_site":i,"status":"unmatched_to_global_rail_source",
                            "original_source_track_offset_m":None})
            continue
        dist,seg_id,idx=ds[0]
        segment=segments[idx]
        chain=segment["original_connector_chain"]
        # At sites near junctions and parallel tracks, identify ambiguous
        # *physical connected-component* mappings instead of forcing original line names.
        chosen_component=base_components.get(chain[0][1])
        competitor=set()
        for dist2,_,other_idx in ds[1:]:
            if dist2>dist+COMPETING_CANDIDATE_M:break
            other=segments[other_idx]["original_connector_chain"]
            alt=base_components.get(other[0][1])
            if alt!=chosen_component:competitor.add(alt)
        ambiguous=bool(competitor)
        if ambiguous:
            attached.append({"source_site":i,"status":"ambiguous_multiple_disconnected_source_tracks",
                "original_source_track_offset_m":round(dist,3),
                "source_competing_connected_components_n":len(competitor)+1})
            continue
        fraction=float(segment["geom"].project(pt,normalized=True))
        ix=0
        while ix<len(chain)-1 and fraction>chain[ix+1][0]:
            ix+=1
        ix=min(ix,len(chain)-2)
        (start,a),(end,b)=chain[ix],chain[ix+1]
        node="botanical_site_"+str(i)
        # Always add original site as a new node. The pairwise network
        # shortest route is now geodesic/projection length on physical connectors,
        # not a made-up same-author-section path.
        outgraph.add_edge(node,a,weight=max(0.0,(fraction-start)*segment["original_projected_length_m"]))
        outgraph.add_edge(node,b,weight=max(0.0,(end-fraction)*segment["original_projected_length_m"]))
        attached.append({"source_site":i,"status":"matched_to_real_rail_connector_network",
            "original_source_track_offset_m":round(dist,3),
            "source_matched_segment_id_private":seg_id,
            "source_fraction":fraction,
            "source_original_length_m":segment["original_projected_length_m"],
            "source_graph_node":node,
            "source_component_PRIVATE":chosen_component})
    return outgraph,attached

def summarize_matched_botany(sites,matches,site_attach,graph,threshold):
    import networkx as nx
    matched=Counter(x["status"] for x in site_attach)
    selected=[]
    pair_support=Counter()
    for triad in matches:
        i,j,k=(triad[p] for p in ("anchor_index","within_index","other_index"))
        a,b,c=site_attach[i],site_attach[j],site_attach[k]
        if not all(x["status"]=="matched_to_real_rail_connector_network" and
                   x["original_source_track_offset_m"]<=threshold for x in (a,b,c)):
            pair_support["unmatched_or_ambiguous_source_triads"]+=1
            continue
        # original source graphical rail network can have disconnected
        # regions. Never use absent path as zero network distance.
        if len({a["source_component_PRIVATE"],b["source_component_PRIVATE"],c["source_component_PRIVATE"]})>1:
            pair_support["source_original_disconnected_topology_triads"]+=1
            continue
        names=[a["source_graph_node"],b["source_graph_node"],c["source_graph_node"]]
        try:
            source_dist=nx.single_source_dijkstra_path_length(graph,names[0],cutoff=MAX_GRAPH_ROUTE_KM*1000,
                                                             weight="weight")
            da=source_dist[names[1]]/1000
            db=source_dist[names[2]]/1000
        except KeyError:
            pair_support["source_real_network_paths_longer_than_cutoff"]+=1
            continue
        # For two sites assigned to identical physical original source
        # segment, original direct path can be shorter than via either
        # end connector. Do not overstate rail connectivity length.
        def direct(source_A,source_B):
            if source_A["source_matched_segment_id_private"]==source_B["source_matched_segment_id_private"]:
                return (abs(source_A["source_fraction"]-source_B["source_fraction"])
                        *source_A["source_original_length_m"]/1000)
            return math.inf
        da=min(da,direct(a,b))
        db=min(db,direct(a,c))
        if not(da>=0 and db>=0):raise ValueError("Unphysical source rail path length")
        src_geo_same=triad["same_line_source_geodesic_km"]
        src_geo_other=triad["different_line_source_geodesic_km"]
        if (da+1e-3<src_geo_same or db+1e-3<src_geo_other):
            pair_support["rail_distance_below_geodesic_source_precision_rounding"]+=1
        selected.append({"author_section":sites[i]["section"],
          "source_same_km":da,"source_other_km":db,
          "source_geo_same_km":src_geo_same,
          "source_geo_other_km":src_geo_other,
          "source_Jaccard_same":triad["BS_and_BT"]["same_line"],
          "source_Jaccard_other":triad["BS_and_BT"]["different_line"]})
        pair_support["source_network_triad_complete"]+=1
    groups=defaultdict(list)
    for row in selected:groups[row["author_section"]].append(row)
    identified=(len(selected)>=12 and len(groups)>=4)
    result={"source_original_anchor_candidates":len(matches),
       "source_rail_matched_site_n":matched["matched_to_real_rail_connector_network"],
       "source_rail_unmatched_site_n":matched["unmatched_to_global_rail_source"],
       "source_rail_ambiguous_site_n":matched["ambiguous_multiple_disconnected_source_tracks"],
       "source_network_pairing_support":dict(pair_support),
       "source_true_author_sections_with_qualified_network_triads":len(groups),
       "source_original_matched_triads_with_true_network_paths":len(selected),
       "source_met_minimum_support_four_lines_twelve_anchors":identified,
       "source_original_section_name_NOT_identical_to_overall_connected_component":True}
    if not identified:
        result["no_network_botanical_association_due_to_source_support"]=True
        return result
    lines=[]
    for section,rs in groups.items():
        lines.append({"n_source_anchors":len(rs),
           "Jaccard_advantage":statistics.mean(x["source_Jaccard_same"]-x["source_Jaccard_other"] for x in rs),
           "real_network_same_km":statistics.mean(x["source_same_km"] for x in rs),
           "real_network_other_km":statistics.mean(x["source_other_km"] for x in rs),
           "geodesic_same_km":statistics.mean(x["source_geo_same_km"] for x in rs),
           "geodesic_other_km":statistics.mean(x["source_geo_other_km"] for x in rs)})
    result["source_section_equal_jaccard_same_minus_other"]=statistics.mean(x["Jaccard_advantage"] for x in lines)
    result["source_sections_with_positive_jaccard_advantage"]=sum(x["Jaccard_advantage"]>0 for x in lines)
    result["source_section_equal_mean_network_same_km"]=statistics.mean(x["real_network_same_km"] for x in lines)
    result["source_section_equal_mean_network_other_km"]=statistics.mean(x["real_network_other_km"] for x in lines)
    result["source_section_equal_mean_geodesic_same_km"]=statistics.mean(x["geodesic_same_km"] for x in lines)
    result["source_section_equal_mean_geodesic_other_km"]=statistics.mean(x["geodesic_other_km"] for x in lines)
    result["source_geodesic_matched_same_minus_other_plant_similarity_before_graph_all59"]=SOURCE_OLD_GEO_MATCHED_EFFECT
    result["source_new_graph_complete_site_cohort_not_directly_comparable_to_full59"]=len(selected)!=SOURCE_OLD_ANCHORS
    result["not_a_field_randomized_test_of_train_dispersal"]=True
    return result

def run(original,spatial_rail_source_factory):
    from pyproj import Transformer
    from networkx.algorithms.components import number_connected_components
    by_sec,qa=to_original_sites(original)
    sites=[s for line in sorted(by_sec) for s in by_sec[line]]
    if len(sites)!=SOURCE_N_SITES or len({s["section"] for s in sites})!=SOURCE_N_SECTIONS:
        raise ValueError("Original botanical source sites/section metadata changed")
    triads,matchedQA=match_candidates(sites,40.0,2.5,"geo_only")
    if (len(triads)!=SOURCE_OLD_ANCHORS or
         len({x["source_anchor_line"] for x in triads})!=SOURCE_OLD_SOURCE_LINES):
        raise ValueError("Previously source-authenticated original v28 geography triads differ")
    baseline=line_equal_aggregate(triads,matchedQA)
    y=baseline["habitat_comparisons"]["BS_and_BT"]["equal_source_line_same_minus_different_line_similarity"]
    if abs(y-SOURCE_OLD_GEO_MATCHED_EFFECT)>1e-10:
        raise ValueError("Original published v28 same-vs-different source field Jaccard drift")
    from pyproj import Transformer
    tr=Transformer.from_crs("EPSG:4326",SOURCE_UTM,always_xy=True)
    lat=[x["lat"] for x in sites]; lon=[x["lon"] for x in sites]
    bbox=(min(lon)-.12,min(lat)-.12,max(lon)+.12,max(lat)+.12)
    original_rows,readQA=spatial_rail_source_factory(bbox)
    graph,segments,graphQA=build_rail_graph(original_rows,tr)
    if len(segments)<20 or graph.number_of_nodes()<20:raise ValueError("Original source regional rail graph too small")
    original_network_comp=number_connected_components(graph)
    augmented,located=site_attach(graph,segments,sites,tr)
    thresholds=[MAX_DISTANCE_FROM_TRACK_METRES,SENSITIVITY_DISTANCE_M]
    panels={str(th):summarize_matched_botany(sites,triads,located,augmented,th) for th in thresholds}
    no_site_info=[{k:v for k,v in x.items() if k in ("status","original_source_track_offset_m")} for x in located]
    dists=[x["original_source_track_offset_m"] for x in located if x.get("original_source_track_offset_m") is not None]
    return {"version":"v358-true-Overture-physical-connector-rail-network-Slovenia-flora",
       "original_botanical_source_DOI":"10.1007/s44473-026-00155-z",
       "original_botanical_source_401279_bytes_SHA256":SOURCE_ORIGINAL_SHA,
       "original_botanical_sites":len(sites),"original_author_line_sections":len(by_sec),
       "original_v28_same_line_geography_matched_Jaccard_gap_all59":y,
       "original_rail_subtype_source":readQA,
       "genuine_network_source_feature_count":graphQA,
       "genuine_network_source_nodes":graph.number_of_nodes(),
       "genuine_network_source_edges":graph.number_of_edges(),
       "genuine_regional_physical_connector_components_n":original_network_comp,
       "original_field_sites_matched_to_track_source_status":dict(Counter(x["status"] for x in located)),
       "source_original_GPS_to_Overture_track_median_m":statistics.median(dists) if dists else None,
       "source_original_GPS_to_Overture_track_max_reported_m":max(dists) if dists else None,
       "primary_and_accuracy_sensitivity":panels,
       "rail_graph_source_bbox_ONLY_regional_may_break_long_paths":True,
       "no_original_field_GPS_plant_species_or_rail_network_geometry_republished":True,
       "author_15_line_sections_do_NOT_equate_to_Overture_connected_components":True,
       "no_causal_rail_traffic_seed_dispersal_or_global_similarity_claim":True}

def real_overture(bbox):
    import pyarrow.compute as pc
    from overturemaps import record_batch_reader
    source=Counter(); railrows=[]
    errors=[]
    for attempt in (1,2):
        try:
            reader=record_batch_reader("segment",bbox=bbox,release=RELEASE,stac=True)
            if reader is None:raise ValueError("No original published Overture source returned")
            for batch in reader:
                source["original_transport_records_scanned"]+=batch.num_rows
                if source["original_transport_records_scanned"]>MAX_SOURCE_TRANSPORT_ROWS:
                    raise ValueError("Original rail bbox transport source >8M rows")
                selected=batch.filter(pc.equal(batch.column(batch.schema.get_field_index("subtype")),"rail"))
                source["source_rail_subtype_records"]+=selected.num_rows
                if selected.num_rows:
                    railrows.extend(selected.to_pylist())
                if len(railrows)>100000:raise ValueError("Large original regional Overture rail geometry source exceeded 100k segments")
            return railrows,{"source_Overture_release":RELEASE,
                "source_real_bbox_transport_n":source["original_transport_records_scanned"],
                "source_all_rail_subtype_n":source["source_rail_subtype_records"],
                "source_complete":True,
                "source_attempt":attempt,
                "source_original_bbox_not_published_for_rare_botanical_site_privacy":True}
        except Exception as err:
            errors.append({"attempt":attempt,"error_type":type(err).__name__,
                           "short_message":str(err)[:150]})
            source.clear();railrows=[]
            if attempt==1:time.sleep(12)
    raise RuntimeError("Original complete Overture source unavailable: "+str(errors))

def main():
    p=argparse.ArgumentParser()
    for a in ("original-xlsx","source-audit","output"):
        p.add_argument("--"+a,required=True)
    args=p.parse_args()
    original=extract_original_xlsx(args.original_xlsx,args.source_audit)
    result=run(original,real_overture)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
