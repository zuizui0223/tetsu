#!/usr/bin/env python3
"""v27.1 MIREN ROAD within-survey-year and independent-chain robustness.

v27.0 observed corridor/reference within-site beta was 0.71497, vs source
species-frequency plus per-reference-site richness fixed within-ROAD
checkerboard null beta=0.83224 (16 Region codes, 48 roads, 669 original
100m2 sites). Here test whether that pattern survives:
  A) independent MCMC proposal seed/chain (same exact 16/48/669 panel);
  B) more restrictive source ROAD x FIRST-SURVEY-YEAR strata, preserving
     reference species occupancies per ROAD *and* original first survey year,
     with >=4 sites per road-year, >=2 roads per Region.

All specimen source is MD5-authentic 2025 MIREN Zenodo v3, only earliest
complete same site-year original 100m2 Plot1 road / Plot3 interior,
original first-year 2012+ only. Follow v27.0 same algebra/Baselga pairwise
Jaccard partition. Does NOT equate study-specific Kashmir railway stations
with true MIREN ROADS or identify seed dispersal, mowing or adaptation.

Fixed-margin swap Markov chain draws are technically autocorrelated; report
descriptive technical sensitivity, NOT independent random-experiment p-values.
No plot/site-level raw plant data or GPS re-exported.
"""
from __future__ import annotations
import argparse,json
from collections import defaultdict
from pathlib import Path
from statistics import mean

from audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 import source_read
from analyze_original_RAIL_vs_global_ROAD_taxonomic_microhabitat_beta_v269 import (
    road_first_surveys,select_miren_sites
)
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    sampler,finish_null,modern_roads,METRICS,SEED as OLD_SEED
)

DRAWS=79
BURNIN=3000
SPACING=450
SEED2=271118
MAIN_SOURCE_EXPECTED={"regions":16,"roads":48,"sites":669}
OLD_OBS_BETA_JAC=.7149696640083497
OLD_NULL_BETA_JAC=.8322392749739689
OLD_NULL_SHARE=6.714518137309704
MIN_SITES_YEAR=4
MIN_ROADS_REGION=2


def average_m(samples:list[dict])->dict:
    if not samples:raise ValueError("No independently surveyed source strata to average")
    return {metric:mean(z[metric] for z in samples) for metric in METRICS}


def summarize_groups(regional:dict,seed:int,draws:int=DRAWS,
                     burnin:int=BURNIN,spacing:int=SPACING)->dict:
    regional_means={}
    no_accept=0
    group_count=0
    for ri,(reg,roads) in enumerate(sorted(regional.items())):
        road_outputs=[]
        for j,(road,groups) in enumerate(sorted(roads.items())):
            # One road may contain multiple first-survey-year strata
            # in the sensitivity; each group is source true original sites.
            st=[]
            for k,records in enumerate(groups):
                a=[r["corridor"] for r in records]
                b=[r["reference"] for r in records]
                st.append(sampler(a,b,draws=draws,burnin=burnin,spacing=spacing,
                                  seed=seed+1561*ri+53*j+11*k))
                group_count+=1
            no_accept+=sum(s["accepted_swaps_total"]==0 for s in st)
            road_outputs.append({
                "n_sites":sum(s["original_site_count"] for s in st),
                "obs":average_m([s["observed"] for s in st]),
                "null":[average_m([s["null_draws"][q] for s in st])
                        for q in range(draws)],
            })
        regional_means[reg]={
            "n_roads":len(road_outputs),
            "n_original_sites":sum(z["n_sites"] for z in road_outputs),
            "obs":average_m([z["obs"] for z in road_outputs]),
            "null":[average_m([z["null"][q] for z in road_outputs]) for q in range(draws)]
        }
    if not regional_means:
        return {"identified":False,"n_region_codes":0}
    actual=average_m([r["obs"] for r in regional_means.values()])
    null=[average_m([r["null"][q] for r in regional_means.values()]) for q in range(draws)]
    differences=finish_null(actual,null)
    return {
        "identified":True,
        "n_region_codes":len(regional_means),
        "n_roads":sum(r["n_roads"] for r in regional_means.values()),
        "n_sites":sum(r["n_original_sites"] for r in regional_means.values()),
        "n_source_road_year_strata":group_count,
        "n_road_year_strata_without_accepted_swaps":no_accept,
        "source_fixed_site_and_species_frequency_results":differences,
        "regions_below_fixed_frequency_null_beta_jac":sum(
            r["obs"]["beta_jac"]<mean(z["beta_jac"] for z in r["null"])
            for r in regional_means.values()),
        "regions_with_more_shared_taxa_than_fixed_frequency_null":sum(
            r["obs"]["shared_taxa"]>mean(z["shared_taxa"] for z in r["null"])
            for r in regional_means.values()),
        "MCMC_chain_draws_correlated_NOT_field_CI":True
    }


def blocks_by_source_year(first:list[dict])->tuple[dict,dict]:
    """Only true same ORIGINAL ROAD x first survey year strata with >=4 sites."""
    recent=select_miren_sites(first,"first_2012plus")
    grouped=defaultdict(list)
    for row in recent:
        grouped[row["region"],row["road"],row["year"]].append(row)
    blocks=defaultdict(lambda:defaultdict(list))
    for (reg,road,year),rows in grouped.items():
        if len(rows)>=MIN_SITES_YEAR:
            blocks[reg][road].append(sorted(rows,key=lambda x:str(x["site"])))
    qualified={reg:roads for reg,roads in blocks.items() if len(roads)>=MIN_ROADS_REGION}
    qa={
        "n_source_original_first_2012plus_site_rows":len(recent),
        "n_road_year_strata_with_four_original_sites":sum(len(years) for roads in blocks.values() for years in roads.values()),
        "n_road_year_strata_supported_in_two_road_regions":sum(len(g) for roads in qualified.values() for g in roads.values()),
        "n_source_regions_with_min_two_original_roads":len(qualified),
        "n_original_source_roads":sum(len(x) for x in qualified.values()),
        "n_original_spatial_sites":sum(len(y) for roads in qualified.values() for g in roads.values() for y in g),
        "same_author_survey_year_within_each_fixed_frequency_reference_shuffle":True
    }
    return qualified,qa


def primary_groups(first:list[dict])->tuple[dict,dict]:
    panel,support=modern_roads(first)
    # Group each road's original modern first-site-year records together;
    # preserve original 16region 48road 669site source population exactly.
    return {reg:{road:[sites] for road,sites in roads.items()}
            for reg,roads in panel.items()},support


def run(miren:Path,draws:int=DRAWS,burnin:int=BURNIN,spacing:int=SPACING)->dict:
    source,qc=source_read(miren)
    first,firstqa=road_first_surveys(source)
    modern,source_support=primary_groups(first)
    if source_support!=MAIN_SOURCE_EXPECTED:
        raise ValueError("Frozen v27.0 source cohort changed")
    full=summarize_groups(modern,SEED2,draws,burnin,spacing)
    if abs(full["source_fixed_site_and_species_frequency_results"]["beta_jac"]["observed"]-OLD_OBS_BETA_JAC)>1e-12:
        raise ValueError("v27.0 original source real Jaccard observed did not reproduce")
    if abs(full["source_fixed_site_and_species_frequency_results"]["beta_jac"]["frequency_null_mean"]-OLD_NULL_BETA_JAC)>.025:
        raise ValueError("Alternative fixed-margin MCMC chain differs >.025 from source v27.0; mixing requires examination")
    yr,yrqc=blocks_by_source_year(first)
    years=summarize_groups(yr,SEED2+7711,draws,burnin,spacing)
    return {
        "version":"v27.1",
        "original_source_md5_QA":qc,
        "original_source_first_year_support":firstqa,
        "unchanged_original_2012plus_16_region_48_road_669_site_panel":{
            "source_support":source_support,
            "independent_rng_seed_fixed_margin_chain":full,
            "v27_0_original_observed_beta_jac_reproduced_exact_1e_minus_12":True,
            "v27_0_original_chain_null_beta_jac":OLD_NULL_BETA_JAC,
            "v27_0_independent_chain_max_predeclared_tolerance":.025
        },
        "strict_road_by_original_survey_year_frequency_fixed_reference_null":{
            "year_blocked_source_support_QA":yrqc,
            "site_and_original_source_species_frequency_fixed_by_survey_year":years
        },
        "roadsite_source_real_species_lists_not_republished":True,
        "no_railway_site_in_MIREN":True,
        "no_original_year_global_composition_confound_removed_in_full_chain":True,
        "no_causal_train_mowing_gust_genetic_adaptation":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--miren-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    p.add_argument("--burnin",type=int,default=BURNIN)
    p.add_argument("--spacing",type=int,default=SPACING)
    a=p.parse_args()
    output=run(a.miren_csv,a.draws,a.burnin,a.spacing)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for label,sub in (
        ("independent_full",output["unchanged_original_2012plus_16_region_48_road_669_site_panel"]["independent_rng_seed_fixed_margin_chain"]),
        ("survey_year_blocked",output["strict_road_by_original_survey_year_frequency_fixed_reference_null"]["site_and_original_source_species_frequency_fixed_by_survey_year"])
    ):
        print(json.dumps({"version":"v27.1","mode":label,
            "identified":sub["identified"],
            "regions":sub.get("n_region_codes"),
            "roads":sub.get("n_roads"),
            "sites":sub.get("n_sites"),
            "n_road_year_blocks":sub.get("n_source_road_year_strata"),
            "zero_acceptance_strata":sub.get("n_road_year_strata_without_accepted_swaps"),
            "beta_jac":sub.get("source_fixed_site_and_species_frequency_results",{}).get("beta_jac"),
            "shared":sub.get("source_fixed_site_and_species_frequency_results",{}).get("shared_taxa"),
            "n_region_with_below_null_beta":sub.get("regions_below_fixed_frequency_null_beta_jac")
        },ensure_ascii=False))


if __name__=="__main__":main()
