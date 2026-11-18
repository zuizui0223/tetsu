#!/usr/bin/env python3
"""v28.0 Slovenia original railway: same national line versus a DIFFERENT
national line at closely matched real author GPS great-circle distances.

This tests a distinct level from v27.7-27.9 LOCAL BS_i~BT_i coupling.
All comparisons are between two DISTINCT genuine original physical sites:
  anchor i -> same-source-national-line site j
  anchor i -> different-source-national-line site k
j and k selected solely by physical source WGS84 geodesic distances with
predeclared absolute <=2.5 km matching caliper (sensitivity <=5 km)
and max geographic radius 10/20/40/80 km; source plant similarity is NEVER
used to choose candidate sites. No repeated 173 side-pair pseudoreplication.

Source data: Šabić et al (2026), exact SHA-gated 401279-byte original
railway XLSX, 89 original sites/15 source national line identifiers.
Plants are pooled within each genuine source location: original L/R matched
BS, BT plot taxa; output habitat-specific BS-only, BT-only and BS+BT union
similarities to guard arbitrary one-zone selection.

PRIMARY: maximum geographical separation 40km and absolute same/different
geodesic distance difference <=2.5km (geo_only). Fixed original source
line, section label, 2022/23 year and SINGLE/DOUBLE layout sensitivity,
plus tighter 1km/looser 5km geographical caliper checks. Optional
richness-matched site pair sensitivity (candidate BS+BT union richness
difference <=5 taxa, support permitting). Richness matching is based on
observed source plants, not an independent terrain covariate; secondary.

For each anchor location i, use EXACTLY ONE source same-line reference j
and one other-line reference k, chosen before taxon similarity and reuse
same matching for BS-only/BT-only/union. Real source 15 line SECTION
clusters, not 3916 independent pair experiments. Equal anchors per
source section, source section equal overall. Include all supported
source line sections, leave-one-source-section-out descriptive
sensitivity; if <4 source lines or <12 anchor sites, return
NOT IDENTIFIABLE rather than impute an effect or retune thresholds.

This is a source-geography/line LABEL controlled *association* only.
Same line does not prove operational train connection (switches/service
frequency/rail OD absent) and geographical proximity/functional guild
composition and local origin pools can still be confounded. No inference
of wind adaptation, viability, transport, evolution, or railway's
global biogeographic override. Raw source taxa, GPS, station kilometrage,
exact source site/section IDs and pair records not exported.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import (
    to_original_sites,jaccard_similarity
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    haversine_km
)

N_SITES=89
N_LINES=15
N_SAME_UNORDERED=340
DISTANCES=(10.,20.,40.,80.)
CALIPERS=(1.,2.5,5.)
MODES=("geo_only","same_year","same_year_layout","source_richness")
HABITATS=("BS","BT","BS_and_BT")
MIN_ANCHORS=12
MIN_LINES=4


def source_community(site:dict, habitat:str)->set:
    if habitat=="BS":return site["BS"]
    if habitat=="BT":return site["BT"]
    if habitat=="BS_and_BT":return site["BS"]|site["BT"]
    raise ValueError("Unknown original railway source habitat")


def match_candidates(
    sites:list[dict],max_dist:float,caliper:float,mode:str
)->tuple[list[dict],dict]:
    if max_dist<=0 or caliper<=0 or mode not in MODES:
        raise ValueError("Predeclared spatial controls invalid")
    if len(sites)!=N_SITES or len({s["section"] for s in sites})!=N_LINES:
        raise ValueError("Original 89 actual sites / 15 line labels changed")
    # Precompute all physical distances independent of species identities.
    n=len(sites)
    distance=[[0.]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            d=haversine_km(sites[i]["lat"],sites[i]["lon"],
                           sites[j]["lat"],sites[j]["lon"])
            distance[i][j]=distance[j][i]=d
    trials=[]
    source_within_candidates=0
    source_other_candidates=0
    for i,src in enumerate(sites):
        same=[];different=[]
        for j,ref in enumerate(sites):
            if i==j:continue
            d=distance[i][j]
            if d>max_dist:continue
            if mode in ("same_year","same_year_layout") and src["year"]!=ref["year"]:
                continue
            if mode=="same_year_layout" and src["layout"]!=ref["layout"]:
                continue
            t=(d,j)
            if src["section"]==ref["section"]:
                same.append(t)
            else:
                different.append(t)
        source_within_candidates+=len(same)
        source_other_candidates+=len(different)
        if not same or not different:continue
        allowed=[]
        for da,j in same:
            for db,k in different:
                if abs(da-db)>caliper:continue
                if mode=="source_richness":
                    # Secondary: candidate communities have comparable
                    # source ALL-BS+BT richness (not a fitness measurement).
                    ra=len(source_community(sites[j],"BS_and_BT"))
                    rb=len(source_community(sites[k],"BS_and_BT"))
                    if abs(ra-rb)>5:continue
                allowed.append((abs(da-db),max(da,db),(da+db)/2,da,db,j,k))
        if not allowed:continue
        _,_,_,d_same,d_diff,j,k=min(allowed)
        if j==k or i in (j,k) or sites[j]["section"]!=src["section"] or sites[k]["section"]==src["section"]:
            raise ValueError("Plant species similarity may not define same-line vs other-line reference site")
        source_record={
            "anchor_index":i,
            "within_index":j,"other_index":k,
            "source_anchor_line":src["section"],
            "source_anchor_year":src["year"],
            "source_anchor_layout":src["layout"],
            "same_line_source_geodesic_km":d_same,
            "different_line_source_geodesic_km":d_diff,
            "source_pair_physical_distance_gap_km":abs(d_same-d_diff),
            "same_line_source_candidate_count":len(same),
            "other_line_source_candidate_count":len(different)
        }
        for habitat in HABITATS:
            base=source_community(src,habitat)
            near=source_community(sites[j],habitat)
            other=source_community(sites[k],habitat)
            simsame=jaccard_similarity(base,near)
            simother=jaccard_similarity(base,other)
            source_record[habitat]={
                "same_line":simsame,
                "different_line":simother,
                "matched_same_line_minus_other_line":simsame-simother,
                "ref_source_species_richness_same_line":len(near),
                "ref_source_species_richness_different_line":len(other)
            }
        trials.append(source_record)
    qa={
        "original_sites":n,
        "source_lines":len({s["section"] for s in sites}),
        "source_pairs_same_line_undirected":
            sum(sites[i]["section"]==sites[j]["section"]
                for i in range(n) for j in range(i+1,n)),
        "source_distance_limit_km":max_dist,
        "source_geodesic_matching_absolute_caliper_km":caliper,
        "source_eligibility_constraint":mode,
        "source_anchor_sites_with_both_matched_reference_options":len(trials),
        "source_directed_same_line_opportunities_before_matching":source_within_candidates,
        "source_directed_different_line_opportunities_before_matching":source_other_candidates,
        "matched_source_distances_only_not_species_similarity":mode!="source_richness",
        "matching_uses_plant_richness_secondary_only":mode=="source_richness",
    }
    if qa["source_pairs_same_line_undirected"]!=N_SAME_UNORDERED:
        raise ValueError("Original 340 within-section pair source frequency changed")
    return trials,qa


def line_equal_aggregate(trials:list[dict],qa:dict)->dict:
    sections=defaultdict(list)
    for trial in trials:
        sections[trial["source_anchor_line"]].append(trial)
    out={
        "identified":len(sections)>=MIN_LINES and len(trials)>=MIN_ANCHORS,
        "source_anchor_sites":len(trials),
        "source_lines_with_comparable_same_and_other_reference":len(sections),
        "source_anchor_reference_line_cluster_unit":True,
        "source_median_actual_same_line_geodesic_distance_km":
            median(x["same_line_source_geodesic_km"] for x in trials) if trials else None,
        "source_median_actual_other_line_geodesic_distance_km":
            median(x["different_line_source_geodesic_km"] for x in trials) if trials else None,
        "source_median_actual_absolute_distance_gap_km":
            median(x["source_pair_physical_distance_gap_km"] for x in trials) if trials else None,
        "source_max_actual_absolute_distance_gap_km":
            max((x["source_pair_physical_distance_gap_km"] for x in trials),default=None),
        "original_source_matching_QA":qa
    }
    if not trials:return out
    section_summary={}
    for sec,rs in sections.items():
        section_summary[sec]={
            "source_anchor_count":len(rs),
            "mean_source_distance_gap_km":mean(x["source_pair_physical_distance_gap_km"] for x in rs),
            "habitats":{
                habitat:{
                    "mean_same_line_similarity":mean(x[habitat]["same_line"] for x in rs),
                    "mean_different_line_similarity":mean(x[habitat]["different_line"] for x in rs),
                    "source_line_gap":mean(x[habitat]["matched_same_line_minus_other_line"] for x in rs),
                    "mean_reference_richness_difference_same_minus_other":
                        mean(x[habitat]["ref_source_species_richness_same_line"]-
                             x[habitat]["ref_source_species_richness_different_line"] for x in rs)
                } for habitat in HABITATS}
        }
    if not out["identified"]:
        out["stop_rule"]="BELOW_FOUR_SOURCE_LINES_OR_TWELVE_ANCHORS"
        out["no_positive_or_negative_railway_network_claim_due_to_support"]=True
        return out
    overall={}
    for habitat in HABITATS:
        a=[x["habitats"][habitat] for x in section_summary.values()]
        pairs=[x["source_line_gap"] for x in a]
        loovalues=[mean(pairs[:i]+pairs[i+1:]) for i in range(len(pairs))]
        overall[habitat]={
            "equal_original_line_mean_same_line_similarity":
                mean(x["mean_same_line_similarity"] for x in a),
            "equal_original_line_mean_different_line_similarity":
                mean(x["mean_different_line_similarity"] for x in a),
            "equal_source_line_same_minus_different_line_similarity":
                mean(pairs),
            "source_lines_with_positive_matched_gap":
                sum(x>0 for x in pairs),
            "source_lines_with_negative_matched_gap":
                sum(x<0 for x in pairs),
            "original_line_leave_one_out_range":
                [min(loovalues),max(loovalues)],
            "original_line_leave_one_out_sign_consistent":
                all(x>0 for x in loovalues) if mean(pairs)>0 else
                all(x<0 for x in loovalues) if mean(pairs)<0 else False,
            "equal_original_line_reference_richness_difference":
                mean(x["mean_reference_richness_difference_same_minus_other"] for x in a),
            "NOT_field_randomized_treatment_or_independent_anchor_pairs":True,
        }
    # Aggregate only: no source section names, site IDs or locations.
    out["habitat_comparisons"]=overall
    out["original_line_anchor_distribution_counts_only"]={
        "min":min(map(len,sections.values())),
        "median":median(map(len,sections.values())),
        "max":max(map(len,sections.values()))
    }
    return out


def run(xlsx:Path,audit:Path)->dict:
    source=extract_original_xlsx(xlsx,audit)
    bysection,sourceqa=to_original_sites(source)
    sites=[v for section in sorted(bysection) for v in bysection[section]]
    if len(sites)!=N_SITES:raise ValueError("True author original 89 source sites changed")
    results={}
    for mode in MODES:
        results[mode]={}
        for maximum in DISTANCES:
            results[mode][str(int(maximum))+"km"]={}
            for caliper in CALIPERS:
                matched,qa=match_candidates(sites,maximum,caliper,mode)
                outcome=line_equal_aggregate(matched,qa)
                results[mode][str(int(maximum))+"km"][str(caliper)+"km"]=outcome
    return {
        "version":"v28.0",
        "source_author_original_slovenia_railway_DOI":"10.1007/s44473-026-00155-z",
        "source_original_89_site_15_line_cohort":sourceqa["source_cohort"],
        "predeclared_primary":"geo_only maximum geodesic radius 40km and matched distance difference <=2.5km",
        "source_matched_same_source_section_vs_different_section_geodesic_design":results,
        "corrected_prior_v14_2_null_message":"Geographic+functional matched same-line species similarity loses robust positive excess when source year x track species frequencies are also fixed; this audit never treats same line as train service",
        "same_line_label_not_authorized_actual_physical_train_service_graph":True,
        "not_exact_traffic_freight_OD_or_seed_carriage":True,
        "site_and_line_nesting_not_ordinary_independent_dyads":True,
        "source_genuine_railway_BS_BT_survey_area_equivalence_not_verified":True,
        "source_ephemeral_matched_plant_site_species_and_GPS_NOT_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",required=True,type=Path)
    p.add_argument("--source-audit",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    out=run(a.original_xlsx,a.source_audit)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for mode in MODES:
        for radius in DISTANCES:
            for caliper in CALIPERS:
                r=out["source_matched_same_source_section_vs_different_section_geodesic_design"][mode][str(int(radius))+"km"][str(caliper)+"km"]
                print(json.dumps({
                    "version":"v28.0","mode":mode,"radius_km":radius,"matched_caliper_km":caliper,
                    "identified":r["identified"],
                    "source_lines":r["source_lines_with_comparable_same_and_other_reference"],
                    "source_anchor_locations":r["source_anchor_sites"],
                    "median_geo_gap":r["source_median_actual_absolute_distance_gap_km"],
                    "matched_union_species":r.get("habitat_comparisons",{}).get("BS_and_BT"),
                    "BS_only":r.get("habitat_comparisons",{}).get("BS"),
                    "BT_only":r.get("habitat_comparisons",{}).get("BT")
                },ensure_ascii=False))
if __name__=="__main__":main()
