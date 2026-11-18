#!/usr/bin/env python3
"""v28.2 independent Czech original railway 35 located station FLORAS:
geographically matched same-author-railway-STRETCH versus other stretch.

Source Kutlvašr et al. (2024) Preslia S1+corrected S2, 39 station-scale
original survey floras 2020-22, 761 recent observed taxon names / 6941
positive station×taxon marks. Official Czech station SR70 2026 third-party
DERIVED centroids full-name matched to 35/39 actual original stations,
exact SHA-gated in original source readers. The four unlocated stations
are UNKNOWN, never zero. Centroids are NOT actual field transect GPS.

For each original station i find j from same original author STRETCH and
k from a DIFFERENT original STRETCH whose source WGS84 great-circle
distances from i are within a fixed absolute caliper and max distance.
Match geography independently of taxa (or habitat composition).
Compute taxonomic Jaccard (full 761 original recently observed taxa):
  J(i,j) - J(i,k), one pair/actual anchor station.
The true study has only FIVE stretches; source anchor/site pairs repeat
observations, so equal station within stretch, equal original stretch
estimates; leave-one-STRETCH-out. Never count 595 pair records as
independent biological replicates.

PREDECLARED coarse physical station-centroid ranges 40/80/160/300km
and distance matching difference thresholds 2.5/5/10/20km.
Primary is 160km+10km, before source results.
Sensitivity: original S2 source station area differences, optional
strict |log area_j - log area_k| <= .5 or 1.0; this is author
STATION SCALE, and cannot be assumed equal to the real 100-m BS/BT
Slovenian railway transect or to non-rail natural plots.
When fewer than 4 independently sampled source stretches and 12
matched actual anchor stations, return NOT IDENTIFIABLE.

This is independent source country relative to Slovenian original 89
railway microhabitat locations; different field survey grain/year,
different coordinate precision and biodiversity pool, so don't pool
raw Jaccard effects. Never infer actual train movement, seed
arrival, freight OD or mowing/wind adaptation from section label.
No raw original plant records or station centroid individual details
are serialized.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_Czech_2024_39station_functional_identity_v143 import (
    read_original_czech_publisher,jaccard
)
from audit_Czech_2024_S2_to_SR70_original_station_points_v153 import (
    original_s2_names,sr70_records,match_czech_stations
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km

DISTANCES=(40.,80.,160.,300.)
CALIPERS=(2.5,5.,10.,20.)
AREA_THRESHOLDS=(None,.5,1.)
PRIMARY=("all","160km","10km")
MIN_STRETCHES=4
MIN_STATIONS=12


def source_original_czech_station_floras(s1:Path,s2:Path,sr70:Path)->tuple[list[dict],dict]:
    floras,metadata,floraqa=read_original_czech_publisher(s1,s2)
    original=original_s2_names(s2)
    sr70_rows,geoqa=sr70_records(sr70)
    matched=match_czech_stations(original,sr70_rows)
    bycode={x["original_locality_code"]:x for x in matched}
    if len(bycode)!=39 or len(set(bycode))!=39:
        raise ValueError("Original S2 station code keys not unique")
    source=[]
    no_coordinate=[]
    for code,meta in sorted(metadata.items()):
        r=bycode[code]
        coord=r["matched_czech_2026_source_station_centroid_WGS84"]
        if coord is None:
            no_coordinate.append(code)
            continue
        if meta["stretch"]!=r["source_section"]:
            raise ValueError("Original Preslia author S2 stretch differs from 2026 matching key")
        area=meta["original_source_station_total_area_ha"]
        if area is not None and (area<=0 or not math.isfinite(area)):
            raise ValueError("Source station botanical area invalid")
        source.append({
            "source_site":code,
            "source_stretch":meta["stretch"],
            "lat":coord["lat"],"lon":coord["lon"],
            "area_ha":area,
            "taxa":floras[code]
        })
    if len(source)!=35 or len(no_coordinate)!=4 or len({x["source_stretch"] for x in source})!=5:
        raise ValueError("True Czech original study 39 stations / 35 independently centroids / five stretches gate failed")
    n_same=sum(source[i]["source_stretch"]==source[j]["source_stretch"]
               for i in range(35) for j in range(i+1,35))
    if n_same!=113:raise ValueError("Original 35 geolocated source within-stretch pair count differs from v15.4 (113)")
    return source,{
        "original_author_preslia_plant_source_QA":floraqa,
        "independent_S2_registry_derived_official_station_centroid_QA":geoqa,
        "source_real_original_stations":39,
        "source_geolocated_by_FULL_exact_station_name":35,
        "n_real_source_unlocated_not_imputed_stations":4,
        "n_original_railway_stretches":5,
        "source_same_stretch_unordered_station_pairs":113,
        "source_other_stretch_unordered_pairs":595-113,
        "original_botanical_survey_years":"2020-2022",
        "station_representative_point_year":"2026",
        "source_761_contemporary_not_original_764_paper_full_habitat_taxon_names":True,
    }


def candidates(
    sites:list[dict],max_distance:float,caliper:float,area_delta:float|None
)->tuple[list[dict],dict]:
    if max_distance<=0 or caliper<=0 or area_delta is not None and area_delta<=0:
        raise ValueError("Source geography/area matching controls invalid")
    matched=[]
    n=len(sites)
    same_total=0;other_total=0
    all_dist=[[0.]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            a,b=sites[i],sites[j]
            all_dist[i][j]=all_dist[j][i]=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
    for i,a in enumerate(sites):
        same=[];other=[]
        for j,b in enumerate(sites):
            if i==j or all_dist[i][j]>max_distance:continue
            d=all_dist[i][j]
            if a["source_stretch"]==b["source_stretch"]:same.append((d,j))
            else:other.append((d,j))
        same_total+=len(same);other_total+=len(other)
        if not same or not other:continue
        allowed=[]
        for da,j in same:
            for db,k in other:
                if abs(da-db)>caliper:continue
                if area_delta is not None:
                    x,y=sites[j]["area_ha"],sites[k]["area_ha"]
                    if x is None or y is None:continue
                    if abs(math.log(x)-math.log(y))>area_delta:continue
                allowed.append((abs(da-db),max(da,db),da,db,j,k))
        if not allowed:continue
        _,_,da,db,j,k=min(allowed)
        if i in (j,k) or j==k or sites[j]["source_stretch"]!=a["source_stretch"] or sites[k]["source_stretch"]==a["source_stretch"]:
            raise AssertionError("Source geographic grouping not independent of original plant identities")
        matched.append({
            "source_stretch":a["source_stretch"],
            "source_station_i":i,"source_station_same":j,"source_station_other":k,
            "geographic_source_same_km":da,
            "geographic_source_other_km":db,
            "geographic_matching_gap_km":abs(da-db),
            "original_source_jaccard_same":jaccard(a["taxa"],sites[j]["taxa"]),
            "original_source_jaccard_other":jaccard(a["taxa"],sites[k]["taxa"]),
            "source_all_taxa_same_minus_other":
                jaccard(a["taxa"],sites[j]["taxa"])-jaccard(a["taxa"],sites[k]["taxa"]),
            "original_ref_station_area_same_minus_other":
                (sites[j]["area_ha"]-sites[k]["area_ha"] if
                 sites[j]["area_ha"] is not None and sites[k]["area_ha"] is not None else None),
            "original_ref_station_species_count_same_minus_other":
                len(sites[j]["taxa"])-len(sites[k]["taxa"])
        })
    return matched,{
        "source_n_actual_geolocated_author_stations":len(sites),
        "source_n_original_stretch_labels":len({x["source_stretch"] for x in sites}),
        "source_matched_geographic_maximum_km":max_distance,
        "source_matched_difference_caliper_km":caliper,
        "source_reference_station_area_log_caliper":area_delta,
        "source_matched_anchors":len(matched),
        "directed_same_stretch_candidate_pairs_within_distance":same_total,
        "directed_other_stretch_candidate_pairs_within_distance":other_total,
        "species_identity_never_used_to_select_source_reference_sites":True
    }


def aggregate(pairs:list[dict],qc:dict)->dict:
    groups=defaultdict(list)
    for p in pairs:groups[p["source_stretch"]].append(p)
    n=len(pairs)
    g=len(groups)
    out={
        "identified":n>=MIN_STATIONS and g>=MIN_STRETCHES,
        "source_matched_original_station_anchors":n,
        "source_original_stretches_with_real_geographic_match":g,
        "source_actual_stretch_anchor_count_minmax":
            [min(map(len,groups.values())),max(map(len,groups.values()))] if groups else None,
        "median_original_station_geodesic_pair_gap_km":
            median(p["geographic_matching_gap_km"] for p in pairs) if pairs else None,
        "source_reference_station_richness_diff_equal_original_stretch":
            mean(mean(x["original_ref_station_species_count_same_minus_other"] for x in p)
                 for p in groups.values()) if groups else None,
        "original_source_source_QA":qc
    }
    if not out["identified"]:
        out["stop_reason"]="Fewer than four original author stretches OR twelve true stations with same/other geographic matching"
        return out
    section={stretch:{
        "source_original_stations":len(pp),
        "source_same_sim":mean(x["original_source_jaccard_same"] for x in pp),
        "source_other_sim":mean(x["original_source_jaccard_other"] for x in pp),
        "source_delta":mean(x["source_all_taxa_same_minus_other"] for x in pp)
    } for stretch,pp in groups.items()}
    scores=[x["source_delta"] for x in section.values()]
    l=[mean(scores[:i]+scores[i+1:]) for i in range(len(scores))]
    out.update({
        "real_stretch_equal_source_same_jaccard":mean(x["source_same_sim"] for x in section.values()),
        "real_stretch_equal_source_different_jaccard":mean(x["source_other_sim"] for x in section.values()),
        "real_stretch_equal_same_minus_other_jaccard":mean(scores),
        "original_stretches_positive_source_matched_gap":sum(x>0 for x in scores),
        "original_stretches_negative_source_matched_gap":sum(x<0 for x in scores),
        "original_stretch_leave_one_out_minmax":[min(l),max(l)],
        "original_source_stretches_and_stationpairs_not_independent_countries":True
    })
    return out


def run(s1:Path,s2:Path,sr70:Path)->dict:
    sites,qc=source_original_czech_station_floras(s1,s2,sr70)
    modes={"all":None,"station_area_log05":.5,"station_area_log10":1.}
    result={}
    for mode,area in modes.items():
        result[mode]={}
        for distance in DISTANCES:
            row={}
            for caliper in CALIPERS:
                pairs,qa=candidates(sites,distance,caliper,area)
                row[f"{caliper:g}km"]=aggregate(pairs,qa)
            result[mode][f"{distance:g}km"]=row
    return {
        "version":"v28.2",
        "genuine_source":"Kutlvašr et al (2024), DOI 10.23855/preslia.2024.247",
        "genuine_2026_station_centroid_IS_NOT_original_2020_2022_field_plot_GPS":True,
        "source_QA":qc,
        "predeclared_source_primary":"160km max geographic radius, absolute source distance match <=10km, all 761 original recent taxa",
        "Czech_35_real_station_matched_same_stretch_vs_other_stretch":result,
        "original_source_source_geography_distinct_country_from_Slovenia":True,
        "not_actual_train_service_or_freight_OD_OR_individual_plant_seed_dispersal":True,
        "no_original_taxon_by_station_station_names_or_GPS_serialized":True
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("s1","s2","sr70","output"):p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    result=run(a.s1,a.s2,a.sr70)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for mode in result["Czech_35_real_station_matched_same_stretch_vs_other_stretch"]:
        for radius,calipers in result["Czech_35_real_station_matched_same_stretch_vs_other_stretch"][mode].items():
            for caliper,panel in calipers.items():
                print(json.dumps({
                    "version":"v28.2","mode":mode,"radius":radius,"caliper":caliper,
                    "identified":panel["identified"],
                    "source_original_stretches":panel["source_original_stretches_with_real_geographic_match"],
                    "source_original_station_anchors":panel["source_matched_original_station_anchors"],
                    "source_equal_stretch_same_minus_other":panel.get("real_stretch_equal_same_minus_other_jaccard"),
                    "original_positive_stretches":panel.get("original_stretches_positive_source_matched_gap"),
                    "original_LOO":panel.get("original_stretch_leave_one_out_minmax")
                },ensure_ascii=False))


if __name__=="__main__":main()
