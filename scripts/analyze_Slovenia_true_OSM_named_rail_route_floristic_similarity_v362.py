#!/usr/bin/env python3
"""v362: TRUE infrastructure railway route relations, not author "section" labels.

The 2026 Oct 1 OSM route=railway/tracks member WAY relations are joined using
actual Overture 2026 Nov source way IDs (not fabricated nearest line names)
to 89 genuine 2022/23 Slovenian BS+BT botanical sites.
Compare SAME named OSM route versus DISJOINT route source communities at
BOTH original physical rail network and geographic distance matched windows.
No original route names/IDs, precise GPS, plant taxon lists exported.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import (
    to_original_sites,jaccard_similarity
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from analyze_Slovenia_Overture_connector_network_botany_v358 import (
    build_rail_graph,site_attach,real_overture,SOURCE_UTM,RELEASE
)
from analyze_Slovenia_geod_and_rail_distance_matched_v359 import (
    physical_network_distances,original_community
)
from audit_original_OSM_named_rail_route_Slovenia_v360 import (
    route_memberships,source_rail_name_coverage,PBF_EXPECTED_BYTES,sha256
)
from audit_Slovenia_same_line_geodistance_matched_frequency_null_v281 import (
    source_groups,exact_source_margins,swap_within_source_groups
)

PBF_ORIGINAL_SHA="cce34d7bdfd310238dd302f50b9ad61c3ba9e66bc376856bf2cf56e002dbd9f7"
OSM_ORIGINAL_NAMED_RELATIONS=73
EXPECTED_SITE_COVERAGE=(43,42,1,3)
YEAR="2026-10-01"
SITE_N=89
AUTHOR_SECTIONS=15
RAIL_SEGMENTS=6612
GEO_MAX=40.0
GEO_CAL=2.5
NETWORK_CAL=5.0
PANELS=("unique_named_route","any_named_route_overlap")
MODES=("all_year_layout","exact_year_and_layout")
HABITATS=("BS_and_BT","BS","BT")
N_MIN=12
LINE_MIN=4
SEED=3621124
CHAINS=2
DRAWS=99
BURN=3000
THIN=350

def route_mode(a,b,method):
    """True iff site names original membership directly defines SAME route."""
    if method=="unique_named_route":
        if len(a)!=1 or len(b)!=1:return None
        return next(iter(a))==next(iter(b))
    if method=="any_named_route_overlap":
        if not a or not b:return None
        return bool(a&b)
    raise ValueError("Unknown, non-frozen actual OSM route comparison")

def geometry_only_matches(sites,network_dist,site_routes,method,strict):
    if len(sites)!=SITE_N or len(site_routes)!=SITE_N:
        raise ValueError("Original botanical and source site route cohort drift")
    grouped=[]
    qa=Counter()
    for i,a in enumerate(sites):
        ra=site_routes[i]
        if not ra or (method=="unique_named_route" and len(ra)!=1):
            qa["source_anchor_missing_or_multiple_named_routes"]+=1
            continue
        same=[];other=[]
        for j,b in enumerate(sites):
            if j==i:continue
            classification=route_mode(ra,site_routes[j],method)
            if classification is None:continue
            dist=network_dist.get((i,j))
            if dist is None:continue
            geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
            if geo>GEO_MAX:continue
            if strict and (a["year"]!=b["year"] or a["layout"]!=b["layout"]):
                continue
            (same if classification else other).append((j,geo,dist))
        if not same or not other:
            qa["original_anchor_lacks_both_named_route_candidates"]+=1
            continue
        choices=[]
        for j,gj,nj in same:
            for k,gk,nk in other:
                if abs(gj-gk)>GEO_CAL or abs(nj-nk)>NETWORK_CAL:continue
                choices.append((abs(gj-gk),abs(nj-nk),max(gj,gk),j,k,gj,gk,nj,nk))
        if not choices:
            qa["source_named_route_matches_outside_frozen_geodesic_or_network_caliper"]+=1
            continue
        _,_,_,j,k,gj,gk,nj,nk=min(choices)
        if j==k:raise AssertionError("Source route same/other positions cannot be same botanical site")
        signature=tuple(sorted(ra))
        grouped.append({"anchor_index":i,"within_index":j,"other_index":k,
           "original_anchor_route_signature_PRIVATE":signature,
           "original_named_route_same_geodesic_km":gj,
           "original_named_route_other_geodesic_km":gk,
           "original_named_route_same_network_km":nj,
           "original_named_route_other_network_km":nk,
           "abs_geo_gap_km":abs(gj-gk),"abs_network_gap_km":abs(nj-nk)})
        qa["source_qualified_true_named_route_anchors"]+=1
    qa["source_named_route_signature_groups_with_comparison_n"]=len({
        z["original_anchor_route_signature_PRIVATE"] for z in grouped})
    return grouped,dict(qa)

def score_incidence(rows,triplets,source_sites,habitat,route_as_groups=True):
    by_line=defaultdict(list)
    for t in triplets:
        i,j,k=(t[p] for p in ("anchor_index","within_index","other_index"))
        a,b,c=rows[i],rows[j],rows[k]
        if not a or not b or not c:raise ValueError("Original site species set missing")
        if route_as_groups:
            key=t["original_anchor_route_signature_PRIVATE"]
        else:
            key=source_sites[i]["section"]
        by_line[key].append((jaccard_similarity(a,b),jaccard_similarity(a,c)))
    if not by_line:return None
    perline={name:{
          "source_anchor_n":len(tup),
          "same":statistics.mean(z[0] for z in tup),
          "other":statistics.mean(z[1] for z in tup)}
          for name,tup in by_line.items()}
    diff={k:v["same"]-v["other"] for k,v in perline.items()}
    return {"source_n_independent_route_signatures":len(perline),
       "source_line_equal_same_Jaccard":statistics.mean(v["same"] for v in perline.values()),
       "source_line_equal_other_Jaccard":statistics.mean(v["other"] for v in perline.values()),
       "source_line_equal_same_minus_other_Jaccard":statistics.mean(diff.values()),
       "source_line_signature_positive_count":sum(v>0 for v in diff.values()),
       "source_line_signature_negative_count":sum(v<0 for v in diff.values())}

def source_null(triplets,sites,observed,group_kind,seed):
    full=[set(original_community(s,"BS_and_BT")) for s in sites]
    source_sections=source_groups(sites,group_kind)
    before=exact_source_margins(full,source_sections)
    nulls=[];accept=[]
    for chain in range(CHAINS):
        rows=[set(s) for s in full]
        rng=random.Random(seed+9001*chain)
        nacc=swap_within_source_groups(rows,source_sections,rng,BURN)
        for trial in range(DRAWS):
            nacc+=swap_within_source_groups(rows,source_sections,rng,THIN)
            val=score_incidence(rows,triplets,sites,"BS_and_BT")
            nulls.append(val["source_line_equal_same_minus_other_Jaccard"])
        if exact_source_margins(rows,source_sections)!=before:
            raise AssertionError("Original source site species richness or exact source stratum species occurrence changed")
        accept.append({"source_chain":chain+1,"original_checkerboard_swaps_accepted":nacc})
    v=sorted(nulls);m=statistics.mean(v)
    def q(p):
        f=(len(v)-1)*p;j=int(f)
        return v[j]*(1-(f-j))+v[min(j+1,len(v)-1)]*(f-j)
    return {"original_null_species_frequency_constraints":group_kind,
      "observed_original_same_route_Jaccard_advantage":observed,
      "original_null_gap_mean":m,
      "observed_minus_source_null_mean":observed-m,
      "original_null_technical_q025":q(.025),"original_null_technical_q975":q(.975),
      "technical_upper_rank_NOT_randomized_field_p":
        (1+sum(x>=observed for x in nulls))/(1+len(nulls)),
      "source_correlated_draws_NOT_independent_botanical_studies":len(nulls),
      "original_source_plant_incidence_pool_site_n":len(full),
      "source_survey_year_and_track_layout_species_occurrence_preserved_within_strata":
          group_kind=="within_original_survey_year_and_track_layout",
      "source_checkerboard_acceptance":accept}

def apply_panels(sites,network,route_sets):
    outcomes={}
    for panel in PANELS:
        for mode in MODES:
            strict=(mode=="exact_year_and_layout")
            pairs,qa=geometry_only_matches(sites,network,route_sets,panel,strict)
            kval=f"{panel}__{mode}"
            qualifies=(len(pairs)>=N_MIN and
                qa["source_named_route_signature_groups_with_comparison_n"]>=LINE_MIN)
            out={"source_geometry_only_matched_route_pair_QA":qa,
                 "source_qualifies_12_anchors_4_original_route_signatures":qualifies,
                 "source_named_route_site_triplets_n":len(pairs),
                 "source_line_signature_count":qa["source_named_route_signature_groups_with_comparison_n"],
                 "analysis_is_specific_to_true_OSM_named_infrastructure_routes":True}
            if not qualifies:
                out["biological_outcome"]="NOT_IDENTIFIABLE_LOW_GENUINE_ROUTE_SUPPORT"
                out["no_original_botanical_effect_imputed_when_unsupported"]=True
            else:
                scores={}
                for habitat in HABITATS:
                    rows=[set(original_community(s,habitat)) for s in sites]
                    scores[habitat]=score_incidence(rows,pairs,sites,habitat)
                out["original_named_route_floristic_similarity"]=scores
                out["source_matched_geographical_km_same_vs_other"]=[
                    statistics.mean(x["original_named_route_same_geodesic_km"] for x in pairs),
                    statistics.mean(x["original_named_route_other_geodesic_km"] for x in pairs)]
                out["source_matched_real_rail_network_km_same_vs_other"]=[
                    statistics.mean(x["original_named_route_same_network_km"] for x in pairs),
                    statistics.mean(x["original_named_route_other_network_km"] for x in pairs)]
                obs=scores["BS_and_BT"]["source_line_equal_same_minus_other_Jaccard"]
                out["original_species_occurrence_fixed_nulls"]={
                  kind:source_null(pairs,sites,obs,kind,seed=SEED+717*i)
                  for i,kind in enumerate(("all_original_sites",
                     "within_original_survey_year_and_track_layout"))}
            outcomes[kval]=out
    return outcomes

def run(original,pbf,rail_factory=real_overture):
    from pyproj import Transformer
    bysection,originalQA=to_original_sites(original)
    sites=[s for section in sorted(bysection) for s in bysection[section]]
    if len(sites)!=SITE_N or len(bysection)!=AUTHOR_SECTIONS:
        raise ValueError("Original 89 sites/15 publisher sections source drift")
    if Path(pbf).stat().st_size!=PBF_EXPECTED_BYTES or sha256(pbf)!=PBF_ORIGINAL_SHA:
        raise ValueError("Source official OSM monthly original 2026-10-01 PBF size/SHA drift")
    named,allnamed,osmqa=route_memberships(pbf)
    if osmqa["source_infrastructure_route_relations_named"]!=OSM_ORIGINAL_NAMED_RELATIONS:
        raise ValueError("Original OSM route railway/tracks named relationship source drift")
    original_bbox=(min(x["lon"] for x in sites)-.12,min(x["lat"] for x in sites)-.12,
                   max(x["lon"] for x in sites)+.12,max(x["lat"] for x in sites)+.12)
    original_rail,railqa=rail_factory(original_bbox)
    tr=Transformer.from_crs("EPSG:4326",SOURCE_UTM,always_xy=True)
    G,segments,QA=build_rail_graph(original_rail,tr)
    if len(segments)!=RAIL_SEGMENTS:
        raise ValueError("Original 6612 physical Overture source segments changed")
    augmented,located=site_attach(G,segments,sites,tr)
    by_GERS,coverage=source_rail_name_coverage(original_rail,named,allnamed)
    routes=[]
    for row in located:
        if row["status"]=="matched_to_real_rail_connector_network":
            routes.append(by_GERS[row["source_matched_segment_id_private"]]["named"])
        else:routes.append(set())
    sourceCount=Counter()
    for row,route in zip(located,routes):
        if row["status"]!="matched_to_real_rail_connector_network":
            sourceCount["ambiguous_or_unmatched"]+=1
        elif not route:sourceCount["zero_named_route"]+=1
        elif len(route)==1:sourceCount["exactly_one_named_route"]+=1
        else:sourceCount["multiple_named_routes"]+=1
    if tuple(sourceCount[x] for x in ("exactly_one_named_route","multiple_named_routes",
            "zero_named_route","ambiguous_or_unmatched"))!=EXPECTED_SITE_COVERAGE:
        raise ValueError("Original v360 43/42/1/3 genuine named route coverage source drift")
    network=physical_network_distances(augmented,located)
    panels=apply_panels(sites,network,routes)
    return {"version":"v362-original-Slovenia-true-OSM-named-infrastructure-route-floristic-beta",
        "source_original_biological_DOI":"10.1007/s44473-026-00155-z",
        "source_original_osm_pbf_sha256":PBF_ORIGINAL_SHA,
        "source_original_osm_monthly_snapshot":YEAR,
        "source_real_Overture_release":RELEASE,
        "source_Overture_rail_network_6612_segments":len(segments),
        "source_railway_route_relations_named":osmqa["source_infrastructure_route_relations_named"],
        "source_railway_relation_has_many_to_many_way_membership":True,
        "source_original_site_named_route_status":dict(sourceCount),
        "source_real_route_and_physical_network_matched_panels":panels,
        "source_study_is_ONE_country_not_global_route_botany":True,
        "taxon_identity_did_NOT_select_route_or_geodistance_matched_reference_sites":True,
        "original_actual_railway_train_seed_movement_NOT_observed":True,
        "original_source_botanical_taxa_site_GPS_route_ID_or_rail_way_not_reexported":True}

def main():
    p=argparse.ArgumentParser()
    for key in ("osm-pbf","original-xlsx","source-audit","output"):p.add_argument("--"+key,required=True)
    args=p.parse_args()
    original=extract_original_xlsx(args.original_xlsx,args.source_audit)
    result=run(original,args.osm_pbf)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
