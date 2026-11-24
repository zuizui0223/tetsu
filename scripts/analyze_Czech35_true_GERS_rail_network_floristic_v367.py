#!/usr/bin/env python3
"""v367 genuine independent Czech original station floras (35 source centroids,
5 author railway stretches) × physical GERS WGS84 railway network lengths.

Match same-author-stretch and other-author-stretch stations using BOTH
geographic and real network path distances, before looking at plant species.
Preserve source station richness and original taxon occurrence frequencies
under fixed checkerboard null. NOT a train seed transport causal test.
"""
from __future__ import annotations
import argparse,json,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path

from audit_Czech_2024_true_35_stations_stretch_geodistance_matched_v282 import (
    source_original_czech_station_floras,candidates as old_geo_candidates,
    aggregate as old_aggregate
)
from audit_Czech_35station_distance_matched_species_prevalence_null_v283 import (
    source_groups,exact_counts
)
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import swaps
from audit_original_CZ_SK_station_GERS_topology_v366 import (
    original_overture_rail_source,coordinate_source_window,
    SOURCE_STUDIES,source_rail_gps_coverage
)
from recompute_Slovenia_GERS_true_geodesic_network_v364 import (
    corrected_physical_graph,corrected_source_attachments
)
from analyze_Slovenia_geod_and_rail_distance_matched_v359 import physical_network_distances
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import jaccard_similarity

SOURCE_AUTHORS_DOI="10.23855/preslia.2024.247"
EXPECTED_SOURCE_CZ_35_5=(35,5)
EXPECTED_PREVIOUS_GEO_ONLY_ANCHORS=33
EXPECTED_GEO_ONLY_OBSERVED_JACCARD=-0.00173118063895798
EXPECTED_GERS_CONNECTED_STATION_PAIRS=568
MAX_GEODESIC_KM=160.
MAX_ABS_GEO_GAP_KM=10.
MAX_ABS_GERS_NETWORK_GAP_KM=20.
MIN_MATCHED_STATION_ANCHORS=12
MIN_AUTHOR_STRETCHES=4
NULL_SCOPE=("original_35_stations_total_species_frequency",
            "original_5_stretches_species_frequency")
CHAINS=2
DRAWS=99
BURN=2500
THIN=350
SEED=3671124

def source_geo_km(a,b):
    from pyproj import Geod
    return abs(Geod(ellps="WGS84").inv(a["lon"],a["lat"],b["lon"],b["lat"])[2])/1000.

def matched_geographic_and_physical_rail_triplets(sites,network):
    """Do NOT inspect source botanical species identities while selecting triads."""
    out=[];qa=Counter()
    for i,a in enumerate(sites):
        same=[];other=[]
        for j,b in enumerate(sites):
            if i==j:continue
            rail=network.get((i,j))
            if rail is None:
                qa["source_true_GERS_path_missing"]+=1
                continue
            geo=source_geo_km(a,b)
            if geo>MAX_GEODESIC_KM:continue
            item=(j,geo,rail)
            (same if a["source_stretch"]==b["source_stretch"] else other).append(item)
        if not same or not other:
            qa["original_source_anchors_lack_same_and_other_stretch_options"]+=1
            continue
        qualifying=[]
        for j,gs,ns in same:
            for k,go,no in other:
                dgeo=abs(gs-go)
                dnet=abs(ns-no)
                if dgeo<=MAX_ABS_GEO_GAP_KM and dnet<=MAX_ABS_GERS_NETWORK_GAP_KM:
                    qualifying.append((dgeo,dnet,max(gs,go),j,k,gs,go,ns,no))
        if not qualifying:
            qa["source_original_anchors_no_simultaneous_geo_and_network_caliper_match"]+=1
            continue
        (_,_,_,j,k,gs,go,ns,no)=min(qualifying)
        if j==k or a["source_stretch"]!=sites[j]["source_stretch"] or a["source_stretch"]==sites[k]["source_stretch"]:
            raise AssertionError("True Czech source stretch membership cannot be selected from botanical outcome")
        out.append({"anchor":i,"same":j,"other":k,"original_source_stretch_PRIVATE":a["source_stretch"],
                    "geo_same":gs,"geo_other":go,"network_same":ns,"network_other":no,
                    "abs_geo_gap":abs(gs-go),"abs_network_gap":abs(ns-no)})
    qa["n_real_original_Czech_station_anchors"]=len(out)
    qa["original_author_stretches_supported"]=len({t["original_source_stretch_PRIVATE"] for t in out})
    return out,dict(qa)

def line_equal_gap(rows,matched,sites):
    bystretch=defaultdict(list)
    for m in matched:
        i,j,k=(m[z] for z in ("anchor","same","other"))
        a,b,c=rows[i],rows[j],rows[k]
        if not a or not b or not c:raise ValueError("Source original field botanist plant taxa not present")
        bystretch[sites[i]["source_stretch"]].append(
            jaccard_similarity(a,b)-jaccard_similarity(a,c))
    if len(bystretch)<MIN_AUTHOR_STRETCHES:
        raise ValueError("Original independent source CZ railway stretch support below four")
    byline={group:statistics.mean(v) for group,v in bystretch.items()}
    return statistics.mean(byline.values()),byline

def exact_frequency_null(sites,matched,mode,seed):
    originals=[set(s["taxa"]) for s in sites]
    groups=source_groups(sites,mode)
    target=exact_counts(originals,groups)
    observed,observed_by_stretch=line_equal_gap(originals,matched,sites)
    draws=[];orig_sign=Counter();line_draws=defaultdict(list);chains=[]
    for chain in range(CHAINS):
        rows=[set(x) for x in originals]
        rng=random.Random(seed+7019*chain)
        accepts=0
        def mutate(n):
            accepted=0
            for group in groups:
                if len(group)<2:continue
                subset=[rows[i] for i in group]
                accepted+=swaps(subset,rng,n)
            return accepted
        accepts+=mutate(BURN)
        for _ in range(DRAWS):
            accepts+=mutate(THIN)
            value,stretch=line_equal_gap(rows,matched,sites)
            draws.append(value)
            for k,x in stretch.items():line_draws[k].append(x)
        if exact_counts(rows,groups)!=target:
            raise ValueError("Original Czech 35 field station richness/source taxon frequency margins not preserved")
        chains.append({"chain":chain+1,"source_checkerboard_swaps_accepted":accepts,
                       "proposal_count_per_original_null_group":BURN+DRAWS*THIN})
    null=sorted(draws)
    def quantile(p):
        x=p*(len(null)-1)
        i=int(x)
        return null[i] if i==len(null)-1 else null[i]*(1-x+i)+null[i+1]*(x-i)
    return {"original_source_station_sites_for_frequency_pool_n":len(sites),
      "original_matched_jaccard_same_stretch_minus_other":observed,
      "source_frequency_null_mode":mode,
      "source_null_mean":statistics.mean(draws),
      "observed_minus_fixed_species_frequency_null":observed-statistics.mean(draws),
      "technical_rank_upper_NOT_field_p":(1+sum(x>=observed-1e-12 for x in draws))/(len(draws)+1),
      "technical_q025_q975_NOT_field_CI":[quantile(.025),quantile(.975)],
      "original_qualifying_stretches_observed_exceeding_own_source_null_mean_n":
         sum(value>statistics.mean(line_draws[sec]) for sec,value in observed_by_stretch.items()),
      "original_all35_station_site_richness_and_species_occurrence_margins_exact":True,
      "original_5_stretch_frequency_null_conditions_on_line_signal":mode==NULL_SCOPE[1],
      "technical_correlated_MCMC_draws_not_independent_studies":len(draws),
      "technical_mcmc_chains":chains}

def run(source_s1,source_s2,sr70,rail_factory=original_overture_rail_source):
    from pyproj import Transformer
    source,siteqa=source_original_czech_station_floras(
        Path(source_s1),Path(source_s2),Path(sr70))
    if (len(source),len({x["source_stretch"] for x in source}))!=EXPECTED_SOURCE_CZ_35_5:
        raise ValueError("Expected Czech original source study 35/5 changed")
    old,oldqa=old_geo_candidates(source,160.,10.,None)
    oldagg=old_aggregate(old,oldqa)
    if len(old)!=EXPECTED_PREVIOUS_GEO_ONLY_ANCHORS or abs(
        oldagg["real_stretch_equal_same_minus_other_jaccard"]-
        EXPECTED_GEO_ONLY_OBSERVED_JACCARD)>1e-10:
        raise ValueError("Old Czech v282 source 33 geography-only triads/Jaccard baseline drift")
    bbox=coordinate_source_window(source)
    real,geoqa=rail_factory(bbox)
    T=Transformer.from_crs("EPSG:4326","EPSG:32633",always_xy=True)
    _,G,segments,graphqa,metric=corrected_physical_graph(real,T)
    if len(segments)!=19783 or G.number_of_nodes()!=46570 or G.number_of_edges()!=52014:
        raise ValueError("Original source Czech v366 19783 GERS rail segments source drift")
    outgraph,attached,siteqa2=corrected_source_attachments(
        G,segments,source,T,threshold=300.)
    net=physical_network_distances(outgraph,attached)
    sitecounts=source_rail_gps_coverage(
        [{**x,"source_original_group":x["source_stretch"]} for x in source],
        attached,net)
    if (sitecounts["source_original_true_physical_rail_matched_station_n"]!=35
        or sitecounts["source_real_network_connected_original_station_pairs_within_300km"]!=EXPECTED_GERS_CONNECTED_STATION_PAIRS):
        raise ValueError("True previous original Czech v366 physical connected station support drift")
    matched,qa=matched_geographic_and_physical_rail_triplets(source,net)
    identified=(len(matched)>=MIN_MATCHED_STATION_ANCHORS and qa["original_author_stretches_supported"]>=MIN_AUTHOR_STRETCHES)
    out={"version":"v367-CZ-authentic-Preslia-35station-botany-GERS-geo-network-matched-source",
        "original_publisher_DOI":SOURCE_AUTHORS_DOI,
        "original_publisher_source_S1_S2_SHA256":{
            "S1":siteqa["original_author_preslia_plant_source_QA"]["source_S1_SHA256"],
            "S2":siteqa["original_author_preslia_plant_source_QA"]["source_S2_corrected_erratum_SHA256"]},
        "original_2026_SR70_thirdparty_derived_station_GPS_NOT_field_plot_GPS":True,
        "original_source_35_botanical_station_centroids_n":len(source),
        "original_source_author_rail_stretches_n":5,
        "original_source_Overture_release":geoqa["original_Overture_real_release"],
        "original_GERS_rail_segments":len(segments),
        "original_true_source_rail_connected_station_pairs":sitecounts["source_real_network_connected_original_station_pairs_within_300km"],
        "original_geo_only_previous_v282_baseline":{"source_anchors":len(old),
           "original_author_stretch_Jaccard_gap":oldagg["real_stretch_equal_same_minus_other_jaccard"]},
        "real_geo_and_GERS_network_site_matching_QA":qa,
        "original_source_min12_anchors_4_stretches_supported":identified,
        "new_original_source_network_matched_beta_and_null":None,
        "no_source_station_GPS_raw_taxa_matrix_or_source_GERS_geometry_republished":True,
        "no_train_seed_transport_traffic_or_global_replication_inference":True}
    if identified:
        rows=[set(s["taxa"]) for s in source]
        gap,perstretch=line_equal_gap(rows,matched,source)
        loo=[statistics.mean([v for k,v in perstretch.items() if k!=omit]) for omit in perstretch]
        out["new_original_source_network_matched_beta_and_null"]={
            "original_matched_station_anchors_n":len(matched),
            "original_author_stretches_matched_n":len(perstretch),
            "source_equal_stretch_same_minus_other_Jaccard":gap,
            "source_author_stretches_observed_positive_n":sum(x>0 for x in perstretch.values()),
            "source_equal_stretch_LOO_range":[min(loo),max(loo)],
            "source_matched_geographic_km_same_other":[
                statistics.mean(x["geo_same"] for x in matched),
                statistics.mean(x["geo_other"] for x in matched)],
            "source_matched_physical_GERS_network_km_same_other":[
                statistics.mean(x["network_same"] for x in matched),
                statistics.mean(x["network_other"] for x in matched)],
            "source_site_grouped_exact_original_species_prevalence_nulls":{
                mode:exact_frequency_null(source,matched,mode,SEED+771*i)
                for i,mode in enumerate(NULL_SCOPE)}}
    else:
        out["source_biological_effect_status"]="NOT_IDENTIFIED_support_below_frozen_12_anchors_4_author_stretches"
    return out

def main():
    p=argparse.ArgumentParser()
    for key in ("s1","s2","sr70","output"):p.add_argument("--"+key,required=True)
    x=p.parse_args()
    data=run(x.s1,x.s2,x.sr70)
    dest=Path(x.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(data,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
