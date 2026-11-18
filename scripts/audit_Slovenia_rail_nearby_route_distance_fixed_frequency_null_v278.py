#!/usr/bin/env python3
"""v27.8 genuine Slovenian RAIL local flora under *distance-restricted*
source checkerboard occupancy- and site-richness-fixed null.

Author original Šabić et al. 2026 railway source 89 matched BS/BT
locations, 15 national rail sections; 14 multi-site sections/88 sites.
Original true route kilometer chainage and survey year from authentic
publisher XLSX checksum verification in workflow.

Freeze every true BS source species assemblage. Only reshuffle BT
source species between sites in the SAME original national rail
section and within 5/10/20/50km source route-kilometer separation,
or ALL within-section pairs (v27.3 prior method). Accept a checkerboard
swap only if each BT site has a distinct exclusive taxon; maintain
EXACT source BT richness per site and each taxon's BT incidence
frequency within original section. Survey-year-matched restriction
is a separate fixed sensitivity, not a selected best result.
No sample with >radius is exchanged, no out-of-line source hopping.

The original same-site within-location mean Jaccard dissimilarity is
compared to the fixed-margins spatially-local reference null.
A positive null beta - actual beta is extra local species similarity
beyond this conditional source plant-frequency/richness expectation.

IMPORTANT: with small local windows, an original section may have
zero exchangeable edges or no accepted exchanges; this is UNINFORMATIVE,
not zero ecological effect. Report exact source radius support and
accepted/proposed swaps; compare radii also on the INTERSECTION of
the same informative original sections so no sample-composition
reversal is hidden. Technical 79/99 correlated MCMC draws are not
independent field experiments, and matrix swap mixing is only a
diagnostic. Nearest route chainage is not train OD, freight flows
or actual traversed travel distance.

Do not redistribute raw original EPPO taxa, locations, route
chainages, GPS or site-by-plot species arrays.
"""
from __future__ import annotations
import argparse,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    beta,source_freq,finish_null,METRICS
)

RADII=(5.,10.,20.,50.,math.inf)
MODES=("all","same_survey_year")
DRAWS=99
BURNIN=2500
SPACING=300
SEED=278118
SOURCE_14_SECTION_JAC=.6367636451153651


def eligible_edges(items:list[dict],radius:float,mode:str)->tuple[list[tuple[int,int]],dict]:
    if radius<=0 or mode not in MODES:
        raise ValueError("Source route-distance checkerboard restrictions invalid")
    edges=[]
    all_spatial_pairs=0
    excluded_by_year=0
    source_route_vs_geo_issues=0
    for i in range(len(items)):
        for j in range(i+1,len(items)):
            a,b=items[i],items[j]
            if a["section"]!=b["section"]:
                raise AssertionError("Source reference swap cannot cross different national railway sections")
            gap=abs(a["routekm"]-b["routekm"])
            if gap>radius:continue
            all_spatial_pairs+=1
            geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
            if gap<geo-.5:source_route_vs_geo_issues+=1
            if mode=="same_survey_year" and a["year"]!=b["year"]:
                excluded_by_year+=1
                continue
            edges.append((i,j))
    return edges,{
        "source_true_within_line_pairs_within_route_distance":all_spatial_pairs,
        "source_year_mismatch_edges_excluded":excluded_by_year,
        "source_eligible_nearby_reference_swap_edges":len(edges),
        "source_rail_chainage_less_than_geodesic_by_halfkm_pairs_flagged":source_route_vs_geo_issues,
        "no_nonrail_or_out_of_line_reference_swapping":True
    }


def propose(rows:list[set],edges:list[tuple[int,int]],
            rng:random.Random,n_proposals:int)->int:
    if not edges:return 0
    n=0
    for _ in range(n_proposals):
        i,j=rng.choice(edges)
        left=rows[i]-rows[j]
        right=rows[j]-rows[i]
        if not left or not right:continue
        a=rng.choice(sorted(left));b=rng.choice(sorted(right))
        rows[i].remove(a);rows[i].add(b)
        rows[j].remove(b);rows[j].add(a)
        n+=1
    return n


def original_rail_section_null(items:list[dict],radius:float,mode:str,
                               draws:int=DRAWS,burnin:int=BURNIN,
                               spacing:int=SPACING,seed:int=SEED)->dict:
    if len(items)<2 or draws<20:
        raise ValueError("No genuine field-site original section for local source frequency null")
    original_bs=[set(x["BS"]) for x in items]
    original_bt=[set(x["BT"]) for x in items]
    if any(not x for x in original_bs+original_bt):
        raise ValueError("Original source rail BS/BT plot cannot be missing or zero")
    widths=[len(x) for x in original_bt]
    colfreq=source_freq(original_bt)
    edges,spatialqa=eligible_edges(items,radius,mode)
    observed={m:mean(beta(x,y)[m] for x,y in zip(original_bs,original_bt))
              for m in METRICS}
    # A section with no distance-compatible or exchangeable pairs cannot
    # identify the source-reference spatially-local null.
    if not edges:
        return {"informative":False,"source_original_sites":len(items),
                "source_edges":spatialqa,
                "reason":"NO_DISTANCE_ELIGIBLE_NEARBY_REAL_SITE_PAIR",
                "observed":observed}
    bt=[set(x) for x in original_bt]
    rng=random.Random(seed)
    swaps=propose(bt,edges,rng,burnin)
    technical=[]
    for _ in range(draws):
        swaps+=propose(bt,edges,rng,spacing)
        technical.append({m:mean(beta(a,b)[m] for a,b in zip(original_bs,bt))
                          for m in METRICS})
    if [len(x) for x in bt]!=widths or source_freq(bt)!=colfreq:
        raise ValueError("Spatially local MIREN-null method did NOT preserve Slovenian actual site richness/section taxon occupancy")
    informative=swaps>0
    if not informative:
        return {"informative":False,"source_original_sites":len(items),
                "source_edges":spatialqa,
                "reason":"SOURCE_INCIDENCE_MATRIX_NOT_EXCHANGEABLE_AT_ROUTE_RADIUS",
                "observed":observed}
    return {
        "informative":True,
        "source_original_sites":len(items),
        "source_edges":spatialqa,
        "source_original_BT_richness_and_section_occupancy_exactly_conserved":True,
        "accepted_swaps":swaps,
        "source_swap_proposals":burnin+spacing*draws,
        "low_source_MCMC_exchange_acceptance_flag":swaps<50,
        "observed":observed,
        "null_draws":technical,
        "null_summary":finish_null(observed,technical)
    }


def mean_metric(records:list[dict],m:str)->float:
    return mean(x[m] for x in records)


def output_panel(bysection:dict,*,mode:str,radius:float,draws:int=DRAWS,
                 burnin:int=BURNIN,spacing:int=SPACING,seed:int=SEED)->tuple[dict,dict]:
    all_sections={}
    n_attempted=0
    for k,(section,items) in enumerate(sorted(bysection.items())):
        if len(items)<2:continue
        n_attempted+=1
        val=original_rail_section_null(items,radius,mode,draws,
                                       burnin,spacing,seed=seed+1409*k)
        all_sections[section]=val
    informative={k:v for k,v in all_sections.items() if v["informative"]}
    if not informative:
        return {"identified":False,"n_informative_original_sections":0},all_sections
    observed={m:mean(v["observed"][m] for v in informative.values())
              for m in METRICS}
    randoms=[{m:mean(v["null_draws"][j][m] for v in informative.values())
              for m in METRICS} for j in range(draws)]
    effect=finish_null(observed,randoms)
    low=sum(v["low_source_MCMC_exchange_acceptance_flag"] for v in informative.values())
    qc={
        "identified":len(informative)>=4 and sum(v["source_original_sites"] for v in informative.values())>=12,
        "radius_km":"ALL" if math.isinf(radius) else radius,
        "source_constraint_mode":mode,
        "source_14_original_multisite_sections_attempted":n_attempted,
        "n_informative_original_sections":len(informative),
        "source_original_sites_in_informative_sections":sum(v["source_original_sites"] for v in informative.values()),
        "n_uninformative_original_sections":n_attempted-len(informative),
        "original_source_sections_with_weak_swaps_under_50":low,
        "source_actual_road_chainage_exchangeable_pairs":sum(
            v["source_edges"]["source_eligible_nearby_reference_swap_edges"] for v in informative.values()),
        "n_informative_original_sections_actual_beta_below_null":sum(
            v["null_summary"]["beta_jac"]["observed_minus_frequency_null"]<0
            for v in informative.values()),
        "section_equal_original_observed_and_local_frequency_fixed_null":effect,
        "technical_null_draws_are_MCMC_correlated_not_field_p_values":True
    }
    return qc,all_sections


def intersection_stats(panels:dict,mode:str,draws:int=DRAWS)->dict:
    # Same exact original source sections across *all 5 radii*.
    name=lambda x: "ALL" if math.isinf(x) else f"{x:g}km"
    members=[{k for k,v in panels[mode][name(r)]["source_internal"].items()
               if v["informative"]} for r in RADII]
    shared=set.intersection(*members)
    if not shared:
        return {"identified":False,"original_sections_in_all_radii":0}
    result={}
    for radius in RADII:
        key=name(radius)
        sub=[panels[mode][key]["source_internal"][section] for section in sorted(shared)]
        observed={m:mean(x["observed"][m] for x in sub) for m in METRICS}
        nul=[{m:mean(x["null_draws"][j][m] for x in sub) for m in METRICS}
             for j in range(draws)]
        val=finish_null(observed,nul)
        result[key]={
            "actual_source_original_section_n":len(sub),
            "real_source_locations_n":sum(x["source_original_sites"] for x in sub),
            "original_site_source_beta_jac":val["beta_jac"]["observed"],
            "source_spatially_restricted_null_beta_jac_mean":val["beta_jac"]["frequency_null_mean"],
            "source_local_coupling_excess_nullminusobserved":
                -val["beta_jac"]["observed_minus_frequency_null"],
            "original_section_count_local_coupling_positive":sum(
                x["null_summary"]["beta_jac"]["observed_minus_frequency_null"]<0 for x in sub
            ),
        }
    return {
        "identified":len(shared)>=4,
        "original_sections_in_all_radii":len(shared),
        "same_fixed_original_section_panel_all_radii":True,
        "results":result
    }


def run(xlsx:Path,sourceaudit:Path,draws:int=DRAWS,
        burnin:int=BURNIN,spacing:int=SPACING)->dict:
    original=extract_original_xlsx(xlsx,sourceaudit)
    sections,sourceqa=to_original_sites(original)
    allout={}
    for mode_i,mode in enumerate(MODES):
        allout[mode]={}
        for r_i,radius in enumerate(RADII):
            key="ALL" if math.isinf(radius) else f"{radius:g}km"
            out,internal=output_panel(sections,mode=mode,radius=radius,
                draws=draws,burnin=burnin,spacing=spacing,
                seed=SEED+99*mode_i+17003*r_i)
            allout[mode][key]={"source_aggregate":out,"source_internal":internal}
    all_unrestricted=allout["all"]["ALL"]["source_aggregate"]
    if not all_unrestricted["identified"] or all_unrestricted["n_informative_original_sections"]!=14:
        raise ValueError("Prior v27.3 authentic source 14-original-section unrestricted fixed species-frequency null failed")
    if abs(all_unrestricted["section_equal_original_observed_and_local_frequency_fixed_null"][
        "beta_jac"]["observed"]-.6367636451153651)>1e-10:
        raise ValueError("True source 14-line equal mean observed Jaccard drift from v27.3")
    common={mode:intersection_stats(allout,mode,draws) for mode in MODES}
    # Remove per-original-railway-section internal ephemeral source plant
    # maps/coupled MCMC draws; only section-level aggregate metrics remain.
    compressed={
        mode:{k:v["source_aggregate"] for k,v in allout[mode].items()}
        for mode in MODES
    }
    return {
        "version":"v27.8",
        "source_genuine_Slovenia_author_DOI":"10.1007/s44473-026-00155-z",
        "source_original_89site_15national_line_geographic_QA":sourceqa,
        "source_radii_km":[5,10,20,50,"ALL"],
        "distance_restricted_original_section_BT_species_incidence_fixed_margins":compressed,
        "strict_same_original_section_panel_all_radii":common,
        "source_real_site_richness_and_BT_species_occupancy_fixed":True,
        "no_swaps_beyond_original_chainage_cutoff":True,
        "source_Slovenia_author_true_rail_not_ROAD":True,
        "technical_MCMC_not_independent_field_permutations_or_pvalues":True,
        "no_rail_frequency_freight_OD_mowing_wind_selection_demonstrated":True,
        "no_publisher_site_ID_GPS_EPPO_taxon_matrix_or_original_plant_observations_exported":True
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
    report=run(a.original_xlsx,a.source_audit,a.draws,a.burnin,a.spacing)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for mode,panels in report["distance_restricted_original_section_BT_species_incidence_fixed_margins"].items():
        for label,p in panels.items():
            print(json.dumps({
                "version":"v27.8","mode":mode,"radius":label,
                "identified":p["identified"],
                "source_sections":p.get("n_informative_original_sections"),
                "source_actual_locations":p.get("source_original_sites_in_informative_sections"),
                "source_eligible_original_reference_pairs":p.get("source_actual_road_chainage_exchangeable_pairs"),
                "weakly_mixed_original_sections":p.get("original_source_sections_with_weak_swaps_under_50"),
                "observed_jaccard":p.get("section_equal_original_observed_and_local_frequency_fixed_null",{}).get("beta_jac",{}).get("observed"),
                "null_jaccard":p.get("section_equal_original_observed_and_local_frequency_fixed_null",{}).get("beta_jac",{}).get("frequency_null_mean"),
                "sections_source_beta_below_null":p.get("n_informative_original_sections_actual_beta_below_null")
            },ensure_ascii=False))
    print("V278_COMMON_SUPPORT="+json.dumps(report["strict_same_original_section_panel_all_radii"],ensure_ascii=False))


if __name__=="__main__":
    main()
