#!/usr/bin/env python3
"""v28.3 — INDEPENDENT Czech 35 genuine railway STATIONS (2020–22):
is matched same-author-stretch floristic identity beyond source-species
frequency expected, and is this robust to fixing each study stretch pool?

Source Kutlvašr et al (2024) PRESLIA original S1+corrected S2: 39
stations, 761 actual contemporary taxa / 6941 original marks.
Author independent exact full-name matching to 2026 Czech SR70 third-party
station representative points yields 35 true geolocated stations across
5 author field sampling stretches; source S1/S2/SR70 SHA locked upstream.

Pre-fixed original matching from v28.2:
  MAIN geographical max 160km, absolute same-stretch-vs-other-stretch
       distance difference <=10km; same set of 33 source station anchors
       in 5 source stretches, full published S1 recent 761-taxon pool.
  SOURCE AREA sensitivity: reference station area log-difference <=0.5
       (30 anchors 5 stretches), possible strong sampling area confound.

Each v28.2 source anchor i + true j on same original author sampling
stretch and true k from different stretch are FROZEN before any
species-null permutations. Original Jaccard gap J(i,j)−J(i,k), all original
station's full recorded current taxa.

2x2 checkerboard swaps keep EVERY actual station's exact original
species richness and each recorded species' original occupancy:
  - globally across 35 stations, species frequency preserved but not
    study-stretch taxon frequency (test cross-stretch structure).
  - within each of the 5 genuine original stretches, same strict
    site richness and source stretch-specific species occupancy.
    This conditions on the very original stretch species-pool
    composition under investigation, so may be over-adjustment;
    report as *sensitivity*, not valid control proving no dispersal.

Independent 2 chains×99 draws, each technical MCMC autocorrelated;
not field-experiment p-values. DO NOT reuse station-centroid source
positions as if original botanical GPS, road vs rail matched control,
direct freight OD, train seed capture, evolved wind/grazing resistance.

Output original source support, aggregate test statistics and
within-source-stretch direction counts. Never serialize taxon names,
individual station positions, station IDs, or raw original Excel.
"""
from __future__ import annotations

import argparse,json,random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from audit_Czech_2024_true_35_stations_stretch_geodistance_matched_v282 import (
    source_original_czech_station_floras,candidates,aggregate
)
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    swaps,source_freq
)

CONFIGS={
    "primary":{"dist":160.,"caliper":10.,"area":None},
    "same_1970s_source_station_area":{"dist":160.,"caliper":10.,"area":.5}
}
SCOPE=("original_35_stations_total_species_frequency",
       "original_5_stretches_species_frequency")
DRAWS=99
CHAINS=2
BURNIN=2500
THIN=350
SEED=283118
EXPECTED_PRIMARY_OBS=-0.00173118063895798
EXPECTED_AREA_OBS=.014063812315168824


def source_groups(sites:list[dict],mode:str)->list[list[int]]:
    if mode==SCOPE[0]:
        return [list(range(len(sites)))]
    if mode==SCOPE[1]:
        out=defaultdict(list)
        for i,site in enumerate(sites):
            out[site["source_stretch"]].append(i)
        if len(out)!=5:
            raise ValueError("Original five Czech botanic stretches absent")
        return [z for _,z in sorted(out.items())]
    raise ValueError("Undocumented null source plant frequency stratum")


def exact_counts(rows:list[set],groups:list[list[int]])->list[dict]:
    return [{
        "source_site_richness":[len(rows[i]) for i in group],
        "source_species_occupancies":source_freq([rows[i] for i in group])
    } for group in groups]


def fixed_source_matched_gap(rows:list[set],sites:list[dict],matched:list[dict])->tuple[float,dict]:
    grouped=defaultdict(list)
    for p in matched:
        i,j,k=(p[x] for x in ("source_station_i","source_station_same","source_station_other"))
        a,b,c=rows[i],rows[j],rows[k]
        delta=len(a&b)/len(a|b)-len(a&c)/len(a|c)
        grouped[sites[i]["source_stretch"]].append(delta)
    if len(grouped)!=5:
        raise ValueError("Fewer than five real source stretches supported in frozen full-size main/area matching")
    means={k:mean(v) for k,v in grouped.items()}
    return mean(means.values()),means


def one_source_null(sites:list[dict],matched:list[dict],mode:str,seed:int,
                    draws:int=DRAWS,burnin:int=BURNIN,thin:int=THIN)->dict:
    if draws<50:raise ValueError("At least 50 source technical draws per chain required")
    original=[set(x["taxa"]) for x in sites]
    group=source_groups(sites,mode)
    target=exact_counts(original,group)
    observed,real= fixed_source_matched_gap(original,sites,matched)
    chains=[];values=[];line_null=defaultdict(list)
    for ch in range(CHAINS):
        rows=[set(x) for x in original]
        rng=random.Random(seed+ch*7019)
        acc=0
        def proposals(n:int)->int:
            # Each stratum chooses its own site pairs, SOURCE species
            # identities can move only inside source frequency block.
            total=0
            for indices in group:
                if len(indices)<2:continue
                sets=[rows[i] for i in indices]
                total+=swaps(sets,rng,n)
            return total
        acc+=proposals(burnin)
        local=[]
        for d in range(draws):
            acc+=proposals(thin)
            score,means=fixed_source_matched_gap(rows,sites,matched)
            local.append(score)
            for sec,value in means.items():line_null[sec].append(value)
        if exact_counts(rows,group)!=target:
            raise ValueError("Original source Czech station richness or species occupancy MCMC margin drift")
        values+=local
        chains.append({
            "chain":ch+1,"proposals_per_valid_original_stratum":burnin+draws*thin,
            "accepted_source_true_checkerboard_swaps":acc,
            "chain_mean_original_stretch_Jaccard_gap":mean(local)
        })
    ordered=sorted(values)
    def quantile(p):
        i=(len(ordered)-1)*p
        a=int(i)
        return ordered[a] if a==len(ordered)-1 else ordered[a]*(1-i+a)+ordered[a+1]*(i-a)
    return {
        "source_real_original_station_cohort":len(sites),
        "source_observed_original_stretch_equal_gap":observed,
        "source_frequency_fixed_null_mean":mean(values),
        "original_observed_minus_frequency_null":observed-mean(values),
        "technical_MCMC_q025":quantile(.025),
        "technical_MCMC_q975":quantile(.975),
        "technical_upper_tail_rank_NOT_field_p":
            (1+sum(z>=observed-1e-12 for z in values))/(len(values)+1),
        "source_stretches_actual_above_within_scope_null_mean":
            sum(real[x]>mean(line_null[x]) for x in real),
        "n_source_correlated_MCMC_draws":len(values),
        "two_independent_chain_acceptance_and_mean":chains,
        "original_stretch_treatment_not_actual_rail_transport":True,
        "original_source_richness_and_species_frequency_exactly_fixed":True
    }


def run(s1:Path,s2:Path,sr70:Path)->dict:
    original,qa=source_original_czech_station_floras(s1,s2,sr70)
    res={}
    for ki,(panel,meta) in enumerate(CONFIGS.items()):
        matched,matqa=candidates(original,meta["dist"],meta["caliper"],meta["area"])
        report=aggregate(matched,matqa)
        if not report["identified"] or report["source_original_stretches_with_real_geographic_match"]!=5:
            raise ValueError("Original 5-Czech-stretch source primary geographic support drift")
        expected=EXPECTED_PRIMARY_OBS if ki==0 else EXPECTED_AREA_OBS
        if abs(report["real_stretch_equal_same_minus_other_jaccard"]-expected)>1e-10:
            raise ValueError("Original genuine v28.2 source distance matched difference not reproduced")
        null={}
        for j,scope in enumerate(SCOPE):
            null[scope]=one_source_null(original,matched,scope,SEED+12701*ki+319*j)
        res[panel]={
            "source_station_and_original_stretch_design":{
                "n_original_anchor_sites":len(matched),
                "n_original_field_stretches":report["source_original_stretches_with_real_geographic_match"],
                "distance_gap_median_km":report["median_original_station_geodesic_pair_gap_km"],
                "source_Jaccard_same_line":report["real_stretch_equal_source_same_jaccard"],
                "source_Jaccard_other_line":report["real_stretch_equal_source_different_jaccard"],
                "source_Jaccard_same_minus_other":report["real_stretch_equal_same_minus_other_jaccard"],
                "source_original_stretch_LOO_range":report["original_stretch_leave_one_out_minmax"]},
            "taxon_incidence_source_nulls":null
        }
    return {
        "version":"v28.3",
        "source_original_Czech_2024_Preslia_DOI":"10.23855/preslia.2024.247",
        "same_author_original_Czech_survey_and_station_centroids_QA":qa,
        "source_real_czech_study_matched_source_frequency_null":res,
        "original_5_stretch_conditional_null_can_OVERCONTROL_hypothesized_network_signal":True,
        "independent_Czech_country_study_from_Slovenia_but_distinct_survey_grain":True,
        "source_true_station_centroids_2026_NOT_original_flora_field_coordinates":True,
        "no_train_service_freight_OD_seed_capture_or_genetic_adaptation":True,
        "source_taxon_names_station_coordinates_or_individual_station_plant_observations_NOT_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("s1","s2","sr70","output"):p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args()
    report=run(args.s1,args.s2,args.sr70)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for key,value in report["source_real_czech_study_matched_source_frequency_null"].items():
        for scope,r in value["taxon_incidence_source_nulls"].items():
            print(json.dumps({
                "version":"v28.3","panel":key,"original_species_frequency_null_scope":scope,
                "original_true_station_anchors":value["source_station_and_original_stretch_design"]["n_original_anchor_sites"],
                "source_original_5_stretch_observed_gap":r["source_observed_original_stretch_equal_gap"],
                "source_frequency_null_mean":r["source_frequency_fixed_null_mean"],
                "observed_minus_frequency_null":r["original_observed_minus_frequency_null"],
                "source_stretches_above_null":r["source_stretches_actual_above_within_scope_null_mean"],
                "technical_rank":r["technical_upper_tail_rank_NOT_field_p"]
            },ensure_ascii=False))
if __name__=="__main__":
    main()
