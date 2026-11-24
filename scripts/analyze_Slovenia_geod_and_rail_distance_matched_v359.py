#!/usr/bin/env python3
"""v359: original 89 Slovenian botany sites × genuine physical Overture network;
match botanical sites at the SAME geodesic AND rail network distances,
then test section-membership identity Jaccard advantage as a source
association, NOT as causal seed transport.
"""
from __future__ import annotations
import argparse,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import (
    to_original_sites,jaccard_similarity
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from audit_Slovenia_source_same_line_vs_other_geodistance_matched_v280 import (
    match_candidates,line_equal_aggregate
)
from analyze_Slovenia_Overture_connector_network_botany_v358 import (
    build_rail_graph,site_attach,real_overture,RELEASE,SOURCE_ORIGINAL_SHA,
    SOURCE_UTM,MAX_GRAPH_ROUTE_KM,SOURCE_OLD_GEO_MATCHED_EFFECT
)

SOURCE_N=89
SOURCE_SECTIONS=15
SOURCE_V358_MATCHED=86
SOURCE_V358_AMBIGUOUS=3
GEO_MAX_KM=40.
GEO_CALIPER_KM=2.5
NET_CALIPER_KM=5.
MATCH_MODES=("geo_and_network","geo_network_year_layout")
HABITATS=("BS_and_BT","BS","BT")
MAX_SOURCE_SITES_PER_MODE=89
MIN_ANCHORS=12
MIN_SECTIONS=4

def original_community(site,mode):
    if mode=="BS":return site["BS"]
    if mode=="BT":return site["BT"]
    if mode=="BS_and_BT":return site["BS"]|site["BT"]
    raise ValueError("Unknown publisher original rail microhabitat")

def physical_network_distances(graph,attached,cutoff_km=MAX_GRAPH_ROUTE_KM):
    import networkx as nx
    out={}
    for i,a in enumerate(attached):
        if a.get("status")!="matched_to_real_rail_connector_network":continue
        source_d=nx.single_source_dijkstra_path_length(
           graph,a["source_graph_node"],cutoff=cutoff_km*1000,weight="weight")
        for j,b in enumerate(attached):
            if j==i or b.get("status")!="matched_to_real_rail_connector_network":
                continue
            if a["source_component_PRIVATE"]!=b["source_component_PRIVATE"]:continue
            distance=source_d.get(b["source_graph_node"],math.inf)/1000
            if a["source_matched_segment_id_private"]==b["source_matched_segment_id_private"]:
                direct=(abs(a["source_fraction"]-b["source_fraction"])
                   *a["source_original_length_m"]/1000)
                distance=min(distance,direct)
            if math.isfinite(distance) and 0<=distance<=cutoff_km:
                out[i,j]=distance
    return out

def joint_geo_and_network_triples(sites,network_km,mode):
    if mode not in MATCH_MODES:raise ValueError("Unknown source match mode")
    out=[];qa=Counter()
    for i,a in enumerate(sites):
        same=[];other=[]
        for j,b in enumerate(sites):
            if j==i:continue
            if (i,j) not in network_km:
                qa["no_source_Oveture_shortest_path_or_unmatched"]+=1;continue
            geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
            if geo>GEO_MAX_KM:continue
            if mode=="geo_network_year_layout" and (a["year"]!=b["year"] or a["layout"]!=b["layout"]):
                continue
            item=(j,geo,network_km[i,j])
            (same if a["section"]==b["section"] else other).append(item)
        if not same or not other:
            qa["source_anchors_without_both_label_candidate_types"]+=1
            continue
        eligible=[]
        for j,geo_same,net_same in same:
            for k,geo_other,net_other in other:
                dg=abs(geo_same-geo_other)
                dn=abs(net_same-net_other)
                if dg<=GEO_CALIPER_KM and dn<=NET_CALIPER_KM:
                    eligible.append((dg,dn,max(geo_same,geo_other),j,k,
                                     geo_same,geo_other,net_same,net_other))
        if not eligible:
            qa["source_anchors_no_joint_geo_and_network_caliper_match"]+=1
            continue
        q=min(eligible)
        _,_,_,j,k,ga,gb,na,nb=q
        if j==k or a["section"]!=sites[j]["section"] or a["section"]==sites[k]["section"]:
            raise ValueError("Physical route/author source section labels not independent of set")
        out.append({"anchor":i,"same":j,"different":k,
         "source_author_section":a["section"],
         "source_same_geodesic_km":ga,"source_other_geodesic_km":gb,
         "source_same_network_km":na,"source_other_network_km":nb,
         "source_abs_geographic_gap_km":abs(ga-gb),
         "source_abs_network_gap_km":abs(na-nb)})
    qa["source_original_qualifying_anchors"]=len(out)
    qa["original_sections_qualifying"]=len({x["source_author_section"] for x in out})
    return out,dict(qa)

def botanical_comparison(sites,triplets):
    sections=defaultdict(list)
    for x in triplets:sections[x["source_author_section"]].append(x)
    identifiable=len(triplets)>=MIN_ANCHORS and len(sections)>=MIN_SECTIONS
    out={"source_distinct_original_publisher_transects":len(sites),
         "source_qualified_anchors":len(triplets),
         "source_qualified_author_line_sections":len(sections),
         "source_qualifies_predefined_minimum_12_anchors_4_sections":identifiable,
         "geodesic_distance_caliper_km":GEO_CALIPER_KM,
         "source_true_connector_network_distance_caliper_km":NET_CALIPER_KM}
    if not identifiable:
        out["biological_result"]="NOT_IDENTIFIED_due_to_source_support"
        out["do_not_impute_plant_similarity_gap_zero"]=True
        return out
    results={}
    for habitat in HABITATS:
        lines=[]
        for sec,items in sections.items():
            vals=[]
            for x in items:
                a=original_community(sites[x["anchor"]],habitat)
                b=original_community(sites[x["same"]],habitat)
                c=original_community(sites[x["different"]],habitat)
                vals.append((jaccard_similarity(a,b),jaccard_similarity(a,c)))
            means={
               "mean_same_rail_author_section_Jaccard":statistics.mean(y[0] for y in vals),
               "mean_other_rail_author_section_Jaccard":statistics.mean(y[1] for y in vals),
               "gap":statistics.mean(y[0]-y[1] for y in vals)
            }
            lines.append(means)
        diff=[r["gap"] for r in lines]
        loo=[statistics.mean(diff[:i]+diff[i+1:]) for i in range(len(diff))]
        results[habitat]={
          "source_author_section_equal_same_Jaccard":statistics.mean(r["mean_same_rail_author_section_Jaccard"] for r in lines),
          "source_author_section_equal_other_Jaccard":statistics.mean(r["mean_other_rail_author_section_Jaccard"] for r in lines),
          "source_author_section_equal_same_minus_other_Jaccard":statistics.mean(diff),
          "source_author_sections_positive":sum(v>0 for v in diff),
          "source_author_sections_negative":sum(v<0 for v in diff),
          "source_author_section_LOO_gap_range":[min(loo),max(loo)],
          "source_198_source_PSEUDOREPLICATION_guard":True
        }
    out["source_habitat_comparisons"]=results
    for key in ("source_same_geodesic_km","source_other_geodesic_km",
                "source_same_network_km","source_other_network_km",
                "source_abs_geographic_gap_km","source_abs_network_gap_km"):
        bysection=[statistics.mean(v[key] for v in subset) for subset in sections.values()]
        out["line_equal_mean_"+key]=statistics.mean(bysection)
    out["original_source_same_line_label_IS_NOT_physically_shortest_network_distance_or_train_OD"]=True
    out["field_site_pairs_not_independent_sampling_units"]=True
    return out

def run(original,network_factory=real_overture):
    import networkx as nx
    from pyproj import Transformer
    bysec,sourceQA=to_original_sites(original)
    sites=[s for sec in sorted(bysec) for s in bysec[sec]]
    if len(sites)!=SOURCE_N or len(bysec)!=SOURCE_SECTIONS:
        raise ValueError("Original published 89 study sites or 15 author sections changed")
    prior,qa=match_candidates(sites,GEO_MAX_KM,GEO_CALIPER_KM,"geo_only")
    if len(prior)!=59:raise ValueError("Original v28 geography-only matching source baseline differs")
    original_v28=line_equal_aggregate(prior,qa)
    gap=original_v28["habitat_comparisons"]["BS_and_BT"]["equal_source_line_same_minus_different_line_similarity"]
    if abs(gap-SOURCE_OLD_GEO_MATCHED_EFFECT)>1e-10:
        raise ValueError("v359 original author prior study data not reproducible")

    bbox=(min(x["lon"] for x in sites)-.12,min(x["lat"] for x in sites)-.12,
          max(x["lon"] for x in sites)+.12,max(x["lat"] for x in sites)+.12)
    source_rows,railqa=network_factory(bbox)
    utm=Transformer.from_crs("EPSG:4326",SOURCE_UTM,always_xy=True)
    G,segments,graphqa=build_rail_graph(source_rows,utm)
    augmented,located=site_attach(G,segments,sites,utm)
    nmatch=sum(x["status"]=="matched_to_real_rail_connector_network" for x in located)
    namb=sum(x["status"]=="ambiguous_multiple_disconnected_source_tracks" for x in located)
    if (nmatch,namb)!=(SOURCE_V358_MATCHED,SOURCE_V358_AMBIGUOUS):
        raise ValueError("Original v358 physical connector match drift, cannot compare matched source sites")
    network=physical_network_distances(augmented,located)
    panels={}
    for mode in MATCH_MODES:
        triads,triadqa=joint_geo_and_network_triples(sites,network,mode)
        panels[mode]={"source_matching_QA":triadqa,
                      "source_botany":botanical_comparison(sites,triads)}
    return {"version":"v359-original-SI-botany-joint-real-rail-geodistance-matched",
        "original_real_publisher_DOI":"10.1007/s44473-026-00155-z",
        "original_publisher_XLSX_SHA256":SOURCE_ORIGINAL_SHA,
        "source_original_Slovenia_sites":len(sites),
        "source_original_Slovenia_rail_author_sections":len(bysec),
        "source_real_Overture_rail_release":RELEASE,
        "source_Overture_GERS_physical_network_QA":railqa,
        "source_real_graph_nodes":G.number_of_nodes(),
        "source_real_graph_edges":G.number_of_edges(),
        "source_GERS_segment_count":len(segments),
        "source_plants_mapped_to_real_tracks":nmatch,
        "source_plants_not_uniquely_mapped_to_tracks":namb,
        "original_source_true_network_distinct_pairs_with_paths":len(network),
        "original_geo_only_prior_full59_Jaccard_advantage":gap,
        "original_botanical_survival_after_network_matching":panels,
        "original_site_taxonomy_or_GPS_not_redistributed":True,
        "NOT_independent_original_field_randomization_or_caused_by_train_seed_movement":True,
        "no_global_railway_plant_similarity_conclusion":True}

def main():
    p=argparse.ArgumentParser()
    for arg in ("original-xlsx","source-audit","output"):p.add_argument("--"+arg,required=True)
    args=p.parse_args()
    original=extract_original_xlsx(args.original_xlsx,args.source_audit)
    result=run(original)
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
