#!/usr/bin/env python3
"""v27.0 source authentic fixed-margin null for local transport-corridor flora.

Two genuinely independent publisher FIELD STUDIES are kept separate:
- MIREN mountain ROADS 2012+ first complete plot year, 16 source Region codes /
  48 real roads / 669 original 100-m2 site pairs (source v26.9 modern cohort).
- Kashmir RAIL source 2014 and 2017, 15 station + 12 between-station original
  spatial sites with 4 genuine 100-m2 zones per repeated year.

Per original stratum, the REAL corridor-side species sets stay FROZEN and the
reference-zone incidence matrix is independently shuffled by classic 2x2
checkerboard swaps, preserving EXACTLY (i) EACH original reference site's
species richness, (ii) EACH reference species' occupancy frequency within
original road (MIREN) or station/between stratum×year (Kashmir), and (iii)
number/identity of sites/taxa within the actual original local source pool.
This does not fix local plant traits, source environmental fields or original
individual/site detection. It is a source-conditioned frequency null, not
the causal absence of railway habitat/transport/mowing.

Technically correlated MCMC swap samples are NOT independent field studies;
reported null quantiles are descriptive Monte Carlo diagnostics NOT
experimental CIs or confirmatory p-values. No plot rarefaction: unlike
sitewise richness equalization, original beta_nestedness remains estimable.

Source access/hashes handled by existing SHA-gated original XLSX/CSV readers.
This analysis does not upload any raw plant-by-site rows, taxon names or GPS.
"""
from __future__ import annotations

import argparse,json,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_original_RAIL_vs_global_ROAD_taxonomic_microhabitat_beta_v269 import (
    beta,road_first_surveys,select_miren_sites,
    MIN_SITES,MIN_ROADS,KASHMIR_ZONES,KASHMIR_REFERENCES
)
from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix,matching_sites
from audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 import source_read

DRAWS=99
BURNIN=3000
SPACING=350
SEED=270118
METRICS=("beta_jac","beta_jtu_replacement","beta_jne_nestedness",
         "shared_taxa","corridor_only_taxa","corridor_unshared_fraction")
EXPECTED_ROAD={"regions":16,"roads":48,"sites":669}


def source_freq(rows:list[set])->Counter:
    return Counter(x for s in rows for x in s)


def swaps(rows:list[set],rng:random.Random,proposals:int)->int:
    """Fisher-Yates style Markov checkerboard trades with exact row/col margins."""
    if len(rows)<2:raise ValueError("At least two genuine source spatial sites needed")
    accepted=0
    for _ in range(proposals):
        i,j=rng.sample(range(len(rows)),2)
        left=rows[i]-rows[j]
        right=rows[j]-rows[i]
        if not left or not right:continue
        a=rng.choice(sorted(left));b=rng.choice(sorted(right))
        rows[i].remove(a);rows[i].add(b)
        rows[j].remove(b);rows[j].add(a)
        accepted+=1
    return accepted


def site_mean(corridor:list[set],reference:list[set])->dict:
    if len(corridor)!=len(reference) or len(corridor)<2:
        raise ValueError("Original corridor/reference site pairing unmatched")
    scores=[beta(a,b) for a,b in zip(corridor,reference)]
    return {m:mean(x[m] for x in scores) for m in METRICS}


def sampler(corridor:list[set],reference:list[set],
            *,draws:int=DRAWS,burnin:int=BURNIN,spacing:int=SPACING,
            seed:int=SEED)->dict:
    if draws<10 or burnin<0 or spacing<1:
        raise ValueError("Original incidence frequency null needs valid MCMC settings")
    if not corridor or not reference or any(not a or not b for a,b in zip(corridor,reference)):
        raise ValueError("Original source requires actually recorded positive plots")
    rows=[set(s) for s in reference]
    before_sizes=[len(x) for x in rows]
    before_freq=source_freq(rows)
    observed=site_mean(corridor,rows)
    rng=random.Random(seed)
    accepted=swaps(rows,rng,burnin)
    samples=[]
    for _ in range(draws):
        accepted+=swaps(rows,rng,spacing)
        samples.append(site_mean(corridor,rows))
    if [len(x) for x in rows]!=before_sizes or source_freq(rows)!=before_freq:
        raise ValueError("Original site reference richness or source species frequency margin changed")
    if any(abs(x["beta_jac"]-x["beta_jtu_replacement"]-x["beta_jne_nestedness"])>1e-10
           for x in (observed,*samples)):
        raise AssertionError("Source original Baselga beta partition invalid")
    return {
        "original_site_count":len(rows),
        "original_reference_species_pool":len(before_freq),
        "source_species_frequency_margins_fixed_exact":True,
        "each_original_site_reference_richness_fixed_exact":True,
        "frozen_original_corridor_species_occurrences":True,
        "accepted_swaps_total":accepted,
        "proposed_swap_count":burnin+spacing*draws,
        "observed":observed,
        "null_draws":samples,
    }


def quantile(xs:list[float],q:float)->float:
    v=sorted(xs);z=(len(v)-1)*q;lo=int(z)
    return v[lo] if lo==len(v)-1 else v[lo]*(1-(z-lo))+v[lo+1]*(z-lo)


def finish_null(obs:dict,draws:list[dict])->dict:
    if not draws:raise ValueError("Missing source frequency-preserving null draws")
    result={}
    for metric in METRICS:
        values=[x[metric] for x in draws]
        original=obs[metric]
        result[metric]={
            "observed":original,
            "frequency_null_mean":mean(values),
            "observed_minus_frequency_null":original-mean(values),
            "null_q025_technical":quantile(values,.025),
            "null_q975_technical":quantile(values,.975),
            "upper_tail_rank_correlation_sensitive_NOT_field_p":
                (1+sum(x>=original-1e-12 for x in values))/(len(values)+1),
            "lower_tail_rank_correlation_sensitive_NOT_field_p":
                (1+sum(x<=original+1e-12 for x in values))/(len(values)+1)
        }
    return result


def modern_roads(first:list[dict])->tuple[dict,dict]:
    source=select_miren_sites(first,"first_2012plus")
    byroad=defaultdict(list)
    for x in source:
        byroad[(x["region"],x["road"])].append(x)
    byreg=defaultdict(dict)
    for (reg,road),items in sorted(byroad.items()):
        if len(items)>=MIN_SITES:
            byreg[reg][road]=sorted(items,key=lambda x:str(x["site"]))
    qualified={r:groups for r,groups in byreg.items() if len(groups)>=MIN_ROADS}
    count={
        "regions":len(qualified),
        "roads":sum(len(g) for g in qualified.values()),
        "sites":sum(len(v) for g in qualified.values() for v in g.values())
    }
    if count!=EXPECTED_ROAD:
        raise ValueError(f"Frozen TRUE 2012+ MIREN modern roads 16/48/669 source cohort drift {count}")
    return qualified,count


def road_analysis(source_file:Path,draws:int=DRAWS,
                  burnin:int=BURNIN,spacing:int=SPACING)->dict:
    source,source_qc=source_read(source_file)
    first,first_qc=road_first_surveys(source)
    regional,count=modern_roads(first)
    output={}
    all_draws=[]
    for reg,roads in sorted(regional.items()):
        summary=[]
        for i,(road,rows) in enumerate(sorted(roads.items())):
            corridor=[r["corridor"] for r in rows]
            reference=[r["reference"] for r in rows]
            # Fixed deterministic seed by source region/road ordinal, not
            # environment-derived outcomes.
            st=sampler(corridor,reference,draws=draws,burnin=burnin,
                       spacing=spacing,seed=SEED+1669*len(output)+i*37)
            summary.append(st)
        observed={m:mean(x["observed"][m] for x in summary) for m in METRICS}
        null=[{m:mean(x["null_draws"][d][m] for x in summary)
               for m in METRICS} for d in range(draws)]
        output[reg]={
            "n_original_roads":len(summary),
            "n_original_road_sites":sum(x["original_site_count"] for x in summary),
            "n_source_roads_without_accepted_compositional_swaps":
                sum(x["accepted_swaps_total"]==0 for x in summary),
            "mean_accepted_swap_fraction":mean(x["accepted_swaps_total"]/x["proposed_swap_count"] for x in summary),
            "frequency_matched_per_metric":finish_null(observed,null)
        }
        all_draws.append((observed,null))
    real={m:mean(x[0][m] for x in all_draws) for m in METRICS}
    null=[{m:mean(x[1][d][m] for x in all_draws) for m in METRICS}
          for d in range(draws)]
    overall=finish_null(real,null)
    signs={m:{
        "source_regions_observed_below_frequency_null":sum(
            x["frequency_matched_per_metric"][m]["observed_minus_frequency_null"]<0 for x in output.values()),
        "source_regions_observed_above_frequency_null":sum(
            x["frequency_matched_per_metric"][m]["observed_minus_frequency_null"]>0 for x in output.values())
    } for m in METRICS}
    return {
        "source_original_MIREN_metadata":source_qc,
        "source_original_first_complete_site_year_QA":first_qc,
        "modern_original_source_support":count,
        "observed_minus_source_frequency_fixed_margin_results":overall,
        "source_region_direction_support":signs,
        "source_regions_aggregate_only":output,
        "randomization_unit":"reference plant incidence checkerboard swaps WITHIN same original road; corridor frozen",
        "no_same_site_reference_richness_equalized":True,
        "MIREN_ROADS_not_RAIL":True,
        "all_original_2012plus_sites_used_only_once":True
    }


def kashmir_analysis(source_file:Path,draws:int=DRAWS,
                     burnin:int=BURNIN,spacing:int=SPACING)->dict:
    plots,source_qc=parse_xlsx_matrix(source_file)
    result={}
    for group in ("station","between"):
        sites=matching_sites(plots,group,KASHMIR_ZONES)
        if len(sites)!=(15 if group=="station" else 12):
            raise ValueError("Original genuine Kashmir 27 source complete spatial sites drift")
        result[group]={}
        for year in (2014,2017):
            result[group][str(year)]={}
            for index,ref in enumerate(KASHMIR_REFERENCES):
                corridor=[plots[s,year,"Track"] for s in sites]
                other=[plots[s,year,ref] for s in sites]
                if any(x is None for x in corridor+other):
                    raise ValueError("Actual source railway plot absent")
                null=sampler(corridor,other,draws=draws,
                    burnin=burnin,spacing=spacing,
                    seed=SEED+(0 if group=="station" else 3901)+year*101+index*37)
                result[group][str(year)][ref]={
                    "source_original_spatial_sites":len(sites),
                    "species_frequency_fixed_margin":finish_null(
                        null["observed"],null["null_draws"]),
                    "accepted_swaps":null["accepted_swaps_total"],
                    "proposed_swaps":null["proposed_swap_count"],
                    "ref_species_count":null["original_reference_species_pool"],
                    "year_and_source_group_are_not_train_frequency_or_region":True
                }
    return {
        "source_original_authentic_kashmir_QA":source_qc,
        "source_observed_station_and_between":{"station":15,"between":12},
        "actual_original_stratum_not_a_separate_european_country_study":True,
        "null_requires_same_station_between_year_not_actual_geo_matched_distance":True,
        "results":result
    }


def run(miren:Path,kashmir:Path,draws:int=DRAWS,burnin:int=BURNIN,
        spacing:int=SPACING)->dict:
    road=road_analysis(miren,draws,burnin,spacing)
    rail=kashmir_analysis(kashmir,draws,burnin,spacing)
    return {
        "version":"v27.0",
        "source_original_studies":2,
        "MIREN_original_ROAD_within_road_species_frequency_fixed_reference_null":road,
        "Kashmir_original_RAIL_within_source_group_species_frequency_fixed_reference_null":rail,
        "MCMC_draws_per_original_stratum":draws,
        "MCMC_source_proposals_before_sampling":burnin,
        "MCMC_source_proposals_between_draws":spacing,
        "randomization_preserves_each_real_reference_site_richness_and_ref_species_occupancy":True,
        "MCMC_draws_autocorrelated_not_exact_independent_field_permutations":True,
        "original_species_pools_and_road_regional_climate_not_controlled_between_studies":True,
        "no_rail_road_causal_effect_or_train_movement_OR_plant_heritability_inferred":True,
        "source_individual_site_taxon_observations_not_republished":True
    }


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--miren-csv",type=Path,required=True)
    p.add_argument("--kashmir-xlsx",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    p.add_argument("--burnin",type=int,default=BURNIN)
    p.add_argument("--spacing",type=int,default=SPACING)
    args=p.parse_args()
    result=run(args.miren_csv,args.kashmir_xlsx,args.draws,args.burnin,args.spacing)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    road=result["MIREN_original_ROAD_within_road_species_frequency_fixed_reference_null"]
    print(json.dumps({"version":"v27.0","source_study":"MIREN_ROAD",
      "source_support":road["modern_original_source_support"],
      "source_road_frequency_null":road["observed_minus_source_frequency_fixed_margin_results"],
      "source_direction":road["source_region_direction_support"]},ensure_ascii=False))
    rail=result["Kashmir_original_RAIL_within_source_group_species_frequency_fixed_reference_null"]
    for group in ("station","between"):
        for year in ("2014","2017"):
            x=rail["results"][group][year]["Interior Plot 2"]
            print(json.dumps({"version":"v27.0","source_study":"Kashmir_RAIL",
              "source_group":group,"year":year,
              "source_sites":x["source_original_spatial_sites"],
              "source_frequency_null":x["species_frequency_fixed_margin"],
              "accepted_swaps":x["accepted_swaps"]},ensure_ascii=False))


if __name__=="__main__":main()
