#!/usr/bin/env python3
"""v360: true OSM `route=railway/tracks` membership on physical Overture rails.

Original Geofabrik Slovenia 2026-10-16 PBF, MD5 and file bytes verified
by workflow. Original botanists 2022-23 sampled 89 sites / 15 line sections.
This only verifies named route coverage; does not calculate Jaccard effects,
re-export raw observations, original GPS, OSM ways/relation IDs or tracks.
"""
from __future__ import annotations
import argparse,hashlib,json,re,statistics
from collections import Counter,defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites
from analyze_Slovenia_Overture_connector_network_botany_v358 import (
    real_overture,build_rail_graph,site_attach,source_eligible,SOURCE_UTM,RELEASE
)

PBF_EXPECTED_BYTES=313730042
OSM_SNAPSHOT="2026-10-16"
PBF_URL="https://download.geofabrik.de/europe/slovenia-261016.osm.pbf"
OSM_WAY_PATTERN=re.compile(r"^w([0-9]+)@[0-9]+$")
RAIL_RELATION_ROUTE_VALUES=frozenset(("railway","tracks"))
AUTHOR_SITES_N=89
AUTHOR_RAIL_SECTIONS_N=15
EXPECTED_OVERTURE_SITE_MAPPING=(86,3)

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def pure_source_eligible_railway_relation(tags):
    return (tags.get("type")=="route"
            and tags.get("route") in RAIL_RELATION_ROUTE_VALUES)

def named_source_railway_relation(tags):
    return (pure_source_eligible_railway_relation(tags)
            and bool(str(tags.get("name") or "").strip()
                     or str(tags.get("ref") or "").strip()))

def source_osm_way_refs(record):
    if not isinstance(record,dict):return set()
    refs=set()
    for s in record.get("sources") or ():
        if not isinstance(s,dict) or str(s.get("dataset","")).lower()!="openstreetmap":
            continue
        m=OSM_WAY_PATTERN.fullmatch(str(s.get("record_id") or ""))
        if m:refs.add(int(m.group(1)))
    return refs

def route_memberships(pbf):
    import osmium
    if Path(pbf).stat().st_size!=PBF_EXPECTED_BYTES:
        raise ValueError("Original Geofabrik 2026-10-16 PBF source byte size drift")
    class Handler(osmium.SimpleHandler):
        def __init__(self):
            super().__init__()
            self.count=Counter()
            self.named=defaultdict(set)
            self.any_source_relation_way=defaultdict(set)
        def relation(self,relation):
            tags=dict(relation.tags)
            self.count["original_OSM_relations_scanned"]+=1
            if tags.get("type")=="route" and tags.get("route")=="train":
                self.count["train_service_route_relations_not_infrastructure"]+=1
            if not pure_source_eligible_railway_relation(tags):return
            self.count["source_infrastructure_route_relations_all"]+=1
            is_named=named_source_railway_relation(tags)
            self.count["source_infrastructure_route_relations_named"]+=int(is_named)
            ways=set()
            for member in relation.members:
                if member.type=="w":
                    ways.add(int(member.ref))
            if not ways:
                self.count["source_route_relations_without_way_members"]+=1
            self.count["source_route_way_member_incidences"]+=len(ways)
            for way in ways:
                self.any_source_relation_way[way].add(relation.id)
                if is_named:self.named[way].add(relation.id)
    h=Handler()
    h.apply_file(str(pbf),locations=False)
    qa={**h.count,
        "distinct_original_OSM_ways_in_any_infrastructure_route":len(h.any_source_relation_way),
        "distinct_original_OSM_ways_in_named_infrastructure_routes":len(h.named),
        "OSM_source_ways_belonging_to_multiple_named_routes":
            sum(len(x)>1 for x in h.named.values())}
    return h.named,h.any_source_relation_way,qa

def source_rail_name_coverage(railrows,named,all_routes):
    n=Counter();by_GERS={}
    for seg in railrows:
        if seg.get("subtype")!="rail":continue
        n["source_all_Overture_rail_segments"]+=1
        if not source_eligible(seg):
            n["source_Overture_rail_segments_rejected_by_normal_class_activity"]+=1
            continue
        n["source_eligible_normal_rail_segments"]+=1
        refs=source_osm_way_refs(seg)
        named_route=set().union(*(named.get(w,set()) for w in refs)) if refs else set()
        anyroute=set().union(*(all_routes.get(w,set()) for w in refs)) if refs else set()
        n["source_eligible_original_GERS_with_OSM_way_references"]+=bool(refs)
        n["source_eligible_original_GERS_with_any_OSM_route_relation"]+=bool(anyroute)
        n["source_eligible_original_GERS_with_named_OSM_route_relation"]+=bool(named_route)
        n["source_eligible_original_GERS_with_multiple_named_OSM_route_relations"]+=len(named_route)>1
        n["source_eligible_GERS_original_OSM_way_provenance_count"]+=len(refs)
        segid=seg.get("id")
        if not isinstance(segid,str) or not segid:
            raise ValueError("Original Overture GERS rail ID absent")
        if segid in by_GERS:
            raise ValueError("Source Overture GERS rail ID duplicated")
        by_GERS[segid]={"named":named_route,"any":anyroute,"source_way_n":len(refs)}
    return by_GERS,dict(n)

def site_name_coverage(attached,by_GERS):
    out=Counter()
    sectionmap=defaultdict(list)
    for site in attached:
        out["original_author_sites_total"]+=1
        status=site["status"]
        if status!="matched_to_real_rail_connector_network":
            out["original_author_sites_not_unambiguously_matched_to_physical_rail"]+=1
            continue
        seg=site["source_matched_segment_id_private"]
        if seg not in by_GERS:
            raise ValueError("Matched physical rail GERS segment not source cross-joined")
        row=by_GERS[seg]
        n=len(row["named"])
        out["original_author_sites_matched_to_genuine_overature_physical_rails"]+=1
        out["original_sites_with_zero_named_OSM_route_relations"]+=n==0
        out["original_sites_with_exactly_one_named_OSM_route_relation"]+=n==1
        out["original_sites_with_multiple_named_OSM_route_relations"]+=n>1
        out["original_sites_with_any_OSM_infrastructure_relations"]+=bool(row["any"])
    return dict(out)

def run(pbf,author_original):
    from pyproj import Transformer
    import networkx as nx
    sections,botQA=to_original_sites(author_original)
    sites=[x for line in sorted(sections) for x in sections[line]]
    if len(sites)!=AUTHOR_SITES_N or len(sections)!=AUTHOR_RAIL_SECTIONS_N:
        raise ValueError("Original publisher 89 source sites / 15 source sections missing")
    named,anyroute,osmQA=route_memberships(pbf)
    bbox=(min(x["lon"] for x in sites)-.12,min(x["lat"] for x in sites)-.12,
          max(x["lon"] for x in sites)+.12,max(x["lat"] for x in sites)+.12)
    rail,railQA=real_overture(bbox)
    tr=Transformer.from_crs("EPSG:4326",SOURCE_UTM,always_xy=True)
    graph,segments,graphQA=build_rail_graph(rail,tr)
    augmented,site=site_attach(graph,segments,sites,tr)
    mapped=sum(x["status"]=="matched_to_real_rail_connector_network" for x in site)
    ambiguous=sum(x["status"]=="ambiguous_multiple_disconnected_source_tracks" for x in site)
    if (mapped,ambiguous)!=EXPECTED_OVERTURE_SITE_MAPPING:
        raise ValueError("Original v358 86/89 Overture source botanical site mapping drift")
    byGERS,sourceCoverage=source_rail_name_coverage(rail,named,anyroute)
    sourceSites=site_name_coverage(site,byGERS)
    if (sourceSites["original_author_sites_matched_to_genuine_overature_physical_rails"]!=86
        or sourceSites["original_author_sites_not_unambiguously_matched_to_physical_rail"]!=3):
        raise ValueError("Original botanical source line membership counts disagree with physical source")
    return {"version":"v360-true-OSM-named-infrastructure-railway-relation-coverage",
      "source_pbf_original_Geofabrik_snapshot":OSM_SNAPSHOT,
      "source_pbf_original_official_website":PBF_URL,
      "source_original_PBF_bytes":PBF_EXPECTED_BYTES,
      "source_original_PBF_sha256":sha256(pbf),
      "source_Overture_release":RELEASE,
      "source_Overture_original_bbox_complete":railQA,
      "source_original_publisher_botanical_plots":89,
      "source_original_author_sections":15,
      "source_original_Overture_physical_graph_segments":graphQA["source_usable_rail_segments"],
      "original_OSM_railway_route_relation_source_QA":osmQA,
      "original_Overture_eligible_rail_source_covered_by_OSM_route_relations":sourceCoverage,
      "source_field_site_real_named_route_coverage":sourceSites,
      "physical_connector_component_is_not_real_named_line":True,
      "route_railway_membership_does_not_measure_train_service_or_seed_transport":True,
      "OSM_Geofabrik_2026_10_vs_Overture_2026_11_source_mismatch":True,
      "source_named_routes_many_to_many_not_one_line_per_segment":True,
      "no_botanical_Jaccard_or_global_railway_line_similarity_calculated":True,
      "source_original_ways_names_relationIDs_GPS_species_usernames_and_geometries_not_republished":True}

def main():
    p=argparse.ArgumentParser()
    for name in ("source-osm-pbf","original-xlsx","source-audit","output"):
        p.add_argument("--"+name,required=True)
    args=p.parse_args()
    original=extract_original_xlsx(args.original_xlsx,args.source_audit)
    result=run(args.source_osm_pbf,original)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
