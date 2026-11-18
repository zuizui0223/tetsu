#!/usr/bin/env python3
"""v27.9 source-authentic Slovenia railway plant local coupling:
original rail route chainage vs STRAIGHT-LINE WGS84 distance as candidate
reference-exchange neighborhood, evaluated on identical eligible SITE SETS.

Original publisher Šabić et al. 2026 author XLSX, exactly 89 true sites,
15 real Slovenian national railway sections (14 multi-site), 173 original
same-side BS shoulder–BT toe pairs previously unioned per field location.
No inferred train OD or service network routes.

For each source national line and predetermined radius 5,10,20km:
  E_route={(i,j): |km_i-km_j| <= radius};
  E_geodesic={(i,j): haversine(GPS_i,GPS_j) <= radius}
  Targets = sites with at least one neighbor in BOTH original graphs.
  Observed BS_i–BT_i Jaccard beta averaged over exact same Targets.
  Independently swap BT species within chosen E set, while preserving
  *all* original sites' BT richness and section-specific species
  occupancy frequencies; BS community remains fixed.
  Compare candidate source association:
    excess_route = E_null_route[beta(Targets)] - observed beta(Targets)
    excess_geo = E_null_geo[beta(Targets)] - observed beta(Targets).
  Same-site source target cohort holds across both metrics, eliminating
  sample composition bias due to extra candidate edges.

ALWAYS also report number of source sites with no candidate edge in the
road/geo exchange graph. All sites remain in the fixed source taxon margins;
we score only sites that participate in both neighborhood graphs. This
avoids v27.8 conservative diluted stats from unexchangeable sites.

A positive extra local coupling after 5km route AND 5km geographic
reference null is evidence for local site matching *beyond sampled nearby
community membership*, not for train-mediated seed transport.
A difference between geodesic and route-limited null is NOT network
causation: distances and source site coverage are highly correlated,
and source route km is not a graph of real freight OD. MCMC 79/99
correlated draws are technical, never experimental p-values.

Original locality routekm/GPS, EPPO taxon×site, line-site names, original
source record vectors are never exported.
"""
from __future__ import annotations

import argparse,json,math,random
from collections import defaultdict,Counter
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from audit_Slovenia_rail_nearby_route_distance_fixed_frequency_null_v278 import (
    propose,eligible_edges
)
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    beta,source_freq
)

RADII=(5.,10.,20.)
MODES=("all","same_survey_year")
SYSTEMS=("route","geodesic")
DRAWS=79
BURNIN=2500
SPACING=300
SEED=279118
MIN_SECTIONS=4


def geodesic_edges(items:list[dict],radius:float,mode:str)->list[tuple[int,int]]:
    if radius<=0 or mode not in MODES:raise ValueError("Unapproved original source spatial matching condition")
    edges=[]
    for i in range(len(items)):
        for j in range(i+1,len(items)):
            a,b=items[i],items[j]
            if a["section"]!=b["section"]:
                raise AssertionError("Cannot combine reference rail floras across different national lines")
            d=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
            if d>radius:continue
            if mode=="same_survey_year" and a["year"]!=b["year"]:
                continue
            edges.append((i,j))
    return edges


def endpoints(edges:list[tuple[int,int]])->set[int]:
    return {i for a,b in edges for i in (a,b)}


def null_same_targets(items:list[dict],edges:list[tuple[int,int]],
                      indices:list[int],draws:int,burnin:int,spacing:int,seed:int)->dict:
    if draws<20 or len(items)<2 or len(indices)<2:
        raise ValueError("No actual original source test sites")
    frozen=[set(x["BT"]) for x in items]
    old_freq=source_freq(frozen)
    old_rich=[len(x) for x in frozen]
    observed=mean(beta(items[i]["BS"],frozen[i])["beta_jac"] for i in indices)
    source=[set(x) for x in frozen]
    rng=random.Random(seed)
    accepted=propose(source,edges,rng,burnin)
    nulls=[]
    for d in range(draws):
        accepted+=propose(source,edges,rng,spacing)
        nulls.append(mean(beta(items[i]["BS"],source[i])["beta_jac"] for i in indices))
    if old_rich!=[len(x) for x in source] or source_freq(source)!=old_freq:
        raise ValueError("Original publisher source richness or per-species national-line frequency changed")
    return {
        "source_test_site_count":len(indices),
        "source_swap_edges":len(edges),
        "source_swap_acceptances":accepted,
        "source_total_swap_proposals":burnin+spacing*draws,
        "identified":accepted>0,
        "original_source_observed_jaccard_beta":observed,
        "source_MCMC_null_mean_jaccard_beta":mean(nulls),
        "source_null_minus_observed_beta_local_coupling":mean(nulls)-observed,
        "original_source_site_and_taxon_frequencies_fixed_exact":True,
        "no_source_off_neighborhood_reference_swaps":True,
    }


def section_pair(items:list[dict],radius:float,mode:str,
                 draws:int,burnin:int,spacing:int,seed:int)->dict:
    e_route,_=eligible_edges(items,radius,mode)
    e_geo=geodesic_edges(items,radius,mode)
    targets=sorted(endpoints(e_route)&endpoints(e_geo))
    meta={
        "source_original_section_sites":len(items),
        "source_route_window_candidate_edges":len(e_route),
        "source_geodesic_window_candidate_edges":len(e_geo),
        "source_sites_with_candidate_route_neighborhood":len(endpoints(e_route)),
        "source_sites_with_candidate_geodesic_neighborhood":len(endpoints(e_geo)),
        "same_source_sites_with_BOTH_real_distance_neighborhoods":len(targets),
        "original_sites_without_both_type_neighborhoods":len(items)-len(targets),
        "source_window_radius_km":radius,
        "original_year_requirement":mode,
        "native_site_richness_and_source_taxon_frequency_remains_fixed":True,
    }
    if len(targets)<2:
        return {"identified":False,"support":meta}
    sub={}
    for k,(distance_kind,edges) in enumerate((("route",e_route),("geodesic",e_geo))):
        sub[distance_kind]=null_same_targets(items,edges,targets,draws,burnin,
                                              spacing,seed+701*k)
    if not all(x["identified"] for x in sub.values()):
        return {"identified":False,"support":meta,"one_reference_null_cannot_exchange":True}
    if abs(sub["route"]["original_source_observed_jaccard_beta"]-
           sub["geodesic"]["original_source_observed_jaccard_beta"])>1e-12:
        raise ValueError("Both spatial reference configurations must score the EXACT same genuine original plant site set")
    return {
        "identified":True,
        "support":meta,
        "original_source_same_site_jaccard_beta":sub["route"]["original_source_observed_jaccard_beta"],
        "source_route_limited_reference":sub["route"],
        "source_geodesic_limited_reference":sub["geodesic"],
        "source_extra_local_coupling_route_minus_geo":{
            "local_excess_route":sub["route"]["source_null_minus_observed_beta_local_coupling"],
            "local_excess_geodesic":sub["geodesic"]["source_null_minus_observed_beta_local_coupling"],
            "route_excess_minus_geo_excess":
                sub["route"]["source_null_minus_observed_beta_local_coupling"]-
                sub["geodesic"]["source_null_minus_observed_beta_local_coupling"]
        }
    }


def aggregate_panel(bysection:dict,radius:float,mode:str,draws:int,burnin:int,
                    spacing:int,seed:int)->dict:
    source={}
    all_multisite=0
    for k,(sec,items) in enumerate(sorted(bysection.items())):
        if len(items)<2:continue
        all_multisite+=1
        source[sec]=section_pair(items,radius,mode,draws,burnin,spacing,
                                 seed+131*k)
    usable={sec:v for sec,v in source.items() if v["identified"]}
    if not usable:
        return {"identified":False,"original_true_sections_with_usable_source_null":0}
    g=list(usable.values())
    obs=[x["original_source_same_site_jaccard_beta"] for x in g]
    route=[x["source_route_limited_reference"]["source_null_minus_observed_beta_local_coupling"] for x in g]
    geo=[x["source_geodesic_limited_reference"]["source_null_minus_observed_beta_local_coupling"] for x in g]
    differential=[a-b for a,b in zip(route,geo)]
    n=len(g)
    def loovalues(z):
        return [mean(z[:i]+z[i+1:]) for i in range(len(z))] if len(z)>1 else []
    route_loo=loovalues(route);geo_loo=loovalues(geo);d_loo=loovalues(differential)
    return {
        "identified":n>=MIN_SECTIONS,
        "source_original_multisite_sections_attempted":all_multisite,
        "source_original_railway_sections_with_both_nulls_informative":n,
        "source_original_actual_matched_neighborhood_sites":sum(x["support"]["same_source_sites_with_BOTH_real_distance_neighborhoods"] for x in g),
        "original_all_sites_in_source_sections_with_valid_comparison":sum(x["support"]["source_original_section_sites"] for x in g),
        "original_rail_chainage_exchangeable_edges":sum(x["support"]["source_route_window_candidate_edges"] for x in g),
        "original_geodesic_exchangeable_edges":sum(x["support"]["source_geodesic_window_candidate_edges"] for x in g),
        "original_section_equal_observed_Jaccard_beta":mean(obs),
        "original_section_equal_route_limited_extra_local_coupling":mean(route),
        "original_section_equal_geodesic_limited_extra_local_coupling":mean(geo),
        "original_section_equal_route_minus_geodesic_extra_coupling":mean(differential),
        "original_sections_positive_route_excess":sum(x>0 for x in route),
        "original_sections_positive_geodesic_excess":sum(x>0 for x in geo),
        "original_sections_route_excess_greater_than_geo":sum(x>0 for x in differential),
        "route_excess_leave_one_original_section_out_range":
            [min(route_loo),max(route_loo)] if route_loo else None,
        "geodesic_excess_leave_one_original_section_out_range":
            [min(geo_loo),max(geo_loo)] if geo_loo else None,
        "route_minus_geo_leave_original_section_out_range":
            [min(d_loo),max(d_loo)] if d_loo else None,
        "same_source_site_set_used_to_evaluate_route_and_geodesic_effects":True,
        "not_actual_train_network_OD_or_velocity_or_vehicle_transport":True,
    }


def run(xlsx:Path,audit:Path,draws:int=DRAWS,burnin:int=BURNIN,spacing:int=SPACING)->dict:
    source=extract_original_xlsx(xlsx,audit)
    section,qa=to_original_sites(source)
    out={}
    for mi,mode in enumerate(MODES):
        out[mode]={}
        for ri,r in enumerate(RADII):
            out[mode][f"{r:g}km"]=aggregate_panel(section,r,mode,draws,burnin,
                                                    spacing,SEED+7919*mi+1709*ri)
    if not out["all"]["5km"]["identified"]:
        raise ValueError("Author original source matched 5km source route and geodesic reference panels unsupported")
    return {
        "version":"v27.9",
        "author_original_source_DOI":"10.1007/s44473-026-00155-z",
        "source_original_89_site_15_real_national_rail_line_QA":qa,
        "source_exposure_as_real_publisher_route_km_vs_midpoint_WGS84":True,
        "source_common_original_site_target_set_for_both_distance_metrics":True,
        "source_real_rail_window5_10_20_km_same_site_beta":out,
        "source_small_radius_BS_BT_both_true_rail_habitats":True,
        "MCMC_not_randomized_field_experiment_or_train_transport":True,
        "source_plant_species_names_site_ids_routekm_and_GPS_NOT_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",type=Path,required=True)
    p.add_argument("--source-audit",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    p.add_argument("--burnin",type=int,default=BURNIN)
    p.add_argument("--spacing",type=int,default=SPACING)
    a=p.parse_args()
    output=run(a.original_xlsx,a.source_audit,a.draws,a.burnin,a.spacing)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    for mode,results in output["source_real_rail_window5_10_20_km_same_site_beta"].items():
        for window,r in results.items():
            print(json.dumps({
                "version":"v27.9","mode":mode,"radius":window,
                "identified":r["identified"],
                "original_source_line_sections":r.get("source_original_railway_sections_with_both_nulls_informative"),
                "actual_source_sites_with_both_neighborhoods":r.get("source_original_actual_matched_neighborhood_sites"),
                "original_source_route_window_edges":r.get("original_rail_chainage_exchangeable_edges"),
                "original_source_geodesic_window_edges":r.get("original_geodesic_exchangeable_edges"),
                "observed_beta":r.get("original_section_equal_observed_Jaccard_beta"),
                "route_nullminusobserved":r.get("original_section_equal_route_limited_extra_local_coupling"),
                "geodesic_nullminusobserved":r.get("original_section_equal_geodesic_limited_extra_local_coupling"),
                "route_minus_geodesic_excess":r.get("original_section_equal_route_minus_geodesic_extra_coupling"),
                "source_positive_route_sections":r.get("original_sections_positive_route_excess"),
                "source_positive_geo_sections":r.get("original_sections_positive_geodesic_excess"),
            },ensure_ascii=False))


if __name__=="__main__":
    main()
