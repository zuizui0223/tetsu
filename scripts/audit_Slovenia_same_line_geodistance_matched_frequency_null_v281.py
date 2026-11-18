#!/usr/bin/env python3
"""v28.1 — is the slight source same-national-line plant similarity excess
after PREDECLARED geodesic matching more than site richness and taxon
prevalence predict? Genuine author-original Slovenian RAIL study.

v28.0 independently selected anchor site i, other same-line j, and other
different-line k ONLY by source WGS84 distance (max 40km and absolute
distance difference <=2.5km). Here DO NOT change source matched
site triples after seeing taxonomic outcomes. Each site is represented
by the union of authentic matched BS+BT source taxa from its physical
original location. Freeze matched source triads, NOT their plant records.

Checkerboard 2x2 incidence swaps across 89 true site×species original
taxon sets preserve EXACTLY each original site's species richness and
each species' total occupancy across original source sample.
  A) all genuine source sites (geography/year/track NOT preserved)
  B) source survey YEAR × SINGLE/DOUBLE original track type strata,
     preserving every taxon's original occupancy within each stratum.
Both intentionally DO NOT preserve original national railway line
per-species occurrence: that is the potential explanatory association
under test. The matched source geography remains fixed in both nulls,
but site-specific climate/soil/floristic region and sampling effort
are NOT controlled by random species reshuffles.

Original source anchored same-vs-other-line Jaccard similarity gap is
summarized anchor -> national-line equal, with line LOO. Three
PRE-DECLARED source matching designs:
  PRIMARY: all-year/layout 40km, absolute source geodistance caliper2.5km.
  strict: same-source-year & same SINGLE/DOUBLE track source layout.
  richness: candidate reference communities source union richness delta<=5
            with source distance caliper2.5km (not a causal richness control).
No support retuning after source results. Full original 89×15 authentic
source fixed and 173 original L/R same-side pairs pooled within site.

Report technical 99 correlated MCMC samples and two independently seeded
chains for each source matching design / frequency stratification, no
randomized-field-experiment p-values. A positive source matched
same-line observed gap that is not bigger than within-year×layout
source occupancy null would NOT establish rail network-specific
dispersal; same line label is not route topology or actual train service.

No original plant×site records, taxon names, exact locality coordinates,
publisher railway chainage or matched triples exported.
"""
from __future__ import annotations

import argparse,json,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites
from audit_Slovenia_source_same_line_vs_other_geodistance_matched_v280 import (
    match_candidates,line_equal_aggregate,source_community
)
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    swaps,source_freq
)

CONFIGS={
    "primary":{"mode":"geo_only","max_dist":40.,"caliper":2.5},
    "same_year_layout":{"mode":"same_year_layout","max_dist":40.,"caliper":2.5},
    "source_richness_sensitivity":{"mode":"source_richness","max_dist":40.,"caliper":2.5},
}
NULLS=("all_original_sites","within_original_survey_year_and_track_layout")
CHAIN_COUNT=2
DRAWS=99
BURN=3000
THIN=350
SEED=281118
EXPECTED_OBS={
    "primary":.01967094157833441,
    "same_year_layout":.01116389394622206,
    "source_richness_sensitivity":.028690616483159237
}


def source_groups(sites:list[dict],type_:str)->list[list[int]]:
    if type_=="all_original_sites":
        return [list(range(len(sites)))]
    if type_=="within_original_survey_year_and_track_layout":
        grouped=defaultdict(list)
        for i,site in enumerate(sites):
            grouped[(site["year"],site["layout"])].append(i)
        return [sorted(ids) for key,ids in sorted(grouped.items())]
    raise ValueError("Do not create arbitrary post hoc source flora strata")


def exact_source_margins(rows:list[set],groups:list[list[int]])->list[dict]:
    return [{
        "source_sites":len(group),
        "site_richness":[len(rows[i]) for i in group],
        "species_frequency":source_freq([rows[i] for i in group])
    } for group in groups]


def swap_within_source_groups(
    rows:list[set],groups:list[list[int]],rng:random.Random,n_proposals:int
)->int:
    if n_proposals<0:raise ValueError("Unknown source null draw count")
    active=[g for g in groups if len(g)>=2]
    if not active:return 0
    # Draw genuine site-pair uniforms among each stratum only.
    # No shuffle between different source years or SINGLE/DOUBLE layouts.
    total=0
    for group in active:
        subgroup=[rows[i] for i in group]  # views of original mutable sets
        accepted=swaps(subgroup,rng,n_proposals)
        total+=accepted
        if any(subgroup[j] is not rows[site_i] for j,site_i in enumerate(group)):
            raise AssertionError("MCMC original species sets no longer source-linked")
    return total


def original_matched_gap_for_rows(
    rows:list[set],source_sites:list[dict],matched:list[dict]
)->tuple[float,dict]:
    bysection=defaultdict(list)
    for triad in matched:
        i,j,k=(triad[q] for q in ("anchor_index","within_index","other_index"))
        # SOURCE index selection already frozen by geography, never by flora.
        def jac(a,b):
            return len(a&b)/len(a|b)
        score=jac(rows[i],rows[j])-jac(rows[i],rows[k])
        bysection[source_sites[i]["section"]].append(score)
    values={section:mean(scores) for section,scores in bysection.items()}
    if not values:
        raise ValueError("Genuine source 89-site matched original support absent")
    return mean(values.values()),values


def quantile(vals:list[float],p:float)->float:
    q=sorted(vals);position=(len(q)-1)*p;i=int(position)
    return q[i] if i==len(q)-1 else q[i]*(1-(position-i))+q[i+1]*(position-i)


def experiment(
    source_sites:list[dict],matched:list[dict],groups:list[list[int]],
    seed:int,draws:int=DRAWS,burnin:int=BURN,thin:int=THIN
)->dict:
    if draws<50 or burnin<0 or thin<=0:raise ValueError("Invalid fixed MCMC specifications")
    initial=[set(source_community(s,"BS_and_BT")) for s in source_sites]
    margins=exact_source_margins(initial,groups)
    observed,lines=original_matched_gap_for_rows(initial,source_sites,matched)
    if len(lines)<4 or len(matched)<12:
        raise ValueError("Source matched class comparison cannot support randomization")
    scores=[];chainresults=[];linechanges=defaultdict(list)
    for chain in range(CHAIN_COUNT):
        rows=[set(x) for x in initial]
        rng=random.Random(seed+9193*chain)
        accepted=swap_within_source_groups(rows,groups,rng,burnin)
        local=[]
        for _ in range(draws):
            accepted+=swap_within_source_groups(rows,groups,rng,thin)
            score,individual=original_matched_gap_for_rows(rows,source_sites,matched)
            local.append(score)
            for sec,val in individual.items():
                linechanges[sec].append(val)
        if exact_source_margins(rows,groups)!=margins:
            raise ValueError("Authentic taxon occupancy/source site richness changed in exact-frequency null")
        scores+=local
        chainresults.append({
            "chain":chain+1,
            "source_initial_genuine_sites":len(rows),
            "proposals_per_source_group":burnin+thin*draws,
            "accepted_taxon_checkerboard_swaps":accepted,
            "chain_mean_geodistance_matched_line_gap":mean(local)
        })
    result={
        "source_original_site_count":len(initial),
        "original_lines_with_distance_matched_anchors":len(lines),
        "source_matched_anchor_site_count":len(matched),
        "source_matched_observed_line_gap":observed,
        "source_line_equal_mcmc_null_gap_mean":mean(scores),
        "source_observed_minus_fixed_frequency_null_line_gap":observed-mean(scores),
        "source_MCMC_technical_q025":quantile(scores,.025),
        "source_MCMC_technical_q975":quantile(scores,.975),
        "n_source_correlated_MCMC_draws":len(scores),
        "technical_upper_rank_NOT_real_field_p":
            (1+sum(x>=observed-1e-12 for x in scores))/(len(scores)+1),
        "source_original_lines_where_matched_gap_exceeds_null_mean":
            sum(val>mean(linechanges[sec]) for sec,val in lines.items()),
        "source_original_lines_where_matched_gap_less_than_null_mean":
            sum(val<mean(linechanges[sec]) for sec,val in lines.items()),
        "two_independent_seed_chain_diagnostics":chainresults,
        "site_richness_and_taxon_frequency_preserved_per_defined_group":True,
        "original_source_site_plant_community_richness_not_rarefied":True,
        "mcmc_null_not_field_experiment_or_rail_network_effect":True
    }
    return result


def run(xlsx:Path,audit:Path,draws:int=DRAWS,burnin:int=BURN,thin:int=THIN)->dict:
    source=extract_original_xlsx(xlsx,audit)
    original,qa=to_original_sites(source)
    sites=[s for line in sorted(original) for s in original[line]]
    if len(sites)!=89:raise ValueError("Source original Slovenian 89 railway locations not present")
    outcome={}
    for ci,(panel,settings) in enumerate(CONFIGS.items()):
        candidates,matchqa=match_candidates(sites,settings["max_dist"],settings["caliper"],settings["mode"])
        raw=line_equal_aggregate(candidates,matchqa)
        if not raw["identified"]:
            raise ValueError("Predeclared original source main/strict matched reference population unsupported")
        actual=raw["habitat_comparisons"]["BS_and_BT"][
            "equal_source_line_same_minus_different_line_similarity"]
        if abs(actual-EXPECTED_OBS[panel])>1e-10:
            raise ValueError("Frozen true source v28.0 matched site setup not reproduced")
        outcome[panel]={
            "original_match_support":{k:raw[k] for k in
                ("identified","source_anchor_sites","source_lines_with_comparable_same_and_other_reference",
                 "source_median_actual_absolute_distance_gap_km")},
            "original_source_true_matched_gap":actual,
            "fixed_community_frequency_nulls":{}
        }
        for ni,kind in enumerate(NULLS):
            grouped=source_groups(sites,kind)
            value=experiment(sites,candidates,grouped,
                seed=SEED+100003*ci+1021*ni,
                draws=draws,burnin=burnin,thin=thin)
            outcome[panel]["fixed_community_frequency_nulls"][kind]=value
    return {
        "version":"v28.1",
        "original_source_DOI":"10.1007/s44473-026-00155-z",
        "original_author_railway_source_89site_15line_QA":qa,
        "source_original_baseline_matched_line_geodistance_40km_caliper2_5km":outcome,
        "source_actual_survey_year_and_track_layout_occupancy_when_stratified":True,
        "same_line_membership_not_actual_train_od_or_freight_link":True,
        "all_results_conditional_source_country_sample_and_technical_correlated_draws":True,
        "never_causal_train_dispersal_wind_heritable_adaptation":True,
        "no_source_taxon_site_pairs_gps_km_or_plant_observations_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",required=True,type=Path)
    p.add_argument("--source-audit",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS)
    p.add_argument("--burnin",type=int,default=BURN)
    p.add_argument("--thin",type=int,default=THIN)
    args=p.parse_args()
    result=run(args.original_xlsx,args.source_audit,args.draws,args.burnin,args.thin)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    for panel,v in result["source_original_baseline_matched_line_geodistance_40km_caliper2_5km"].items():
        for null,values in v["fixed_community_frequency_nulls"].items():
            print(json.dumps({
                "version":"v28.1","panel":panel,"null_stratification":null,
                "site_anchor_cohort":v["original_match_support"],
                "observed_Jaccard_line_gap":values["source_matched_observed_line_gap"],
                "frequency_null_mean":values["source_line_equal_mcmc_null_gap_mean"],
                "observed_minus_frequency_null":values["source_observed_minus_fixed_frequency_null_line_gap"],
                "source_lines_observed_greater_than_null":
                    values["source_original_lines_where_matched_gap_exceeds_null_mean"],
                "technical_upper_rank":values["technical_upper_rank_NOT_real_field_p"],
                "chain_results":values["two_independent_seed_chain_diagnostics"]
            },ensure_ascii=False))


if __name__=="__main__":
    main()
