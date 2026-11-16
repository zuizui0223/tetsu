#!/usr/bin/env python3
"""v18.5: within-source taxonomic and lifeform null *jointly*.

The 359-taxon 15-station original Croatian/Bosnian railway flora is
compared to fixed-margin swaps of that SAME railway species pool.
This estimates no regional non-rail pool, transport, adaption, or
causal climate effect. Study n=1; station pairs n=105 are dependent.
"""
from __future__ import annotations
import argparse,json,random,statistics,sys
from pathlib import Path
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_Croatia_Bosnia_15station_lifeform_fixedmargin_null_v184 import (
    incidence,margins,climate_margins,swaps,STATION_MACROCLIMATE,SEED,
    N_CHAINS,DRAWS_PER_CHAIN,BURNIN_PROPOSALS,GAP_PROPOSALS
)
from extract_Croatia_Bosnia_2017_original_rotated_table_matrix_v183 import rows_from_original

def pairwise_metrics(sites,lifeforms):
    if len(sites)!=15:raise ValueError("Publisher 15 sites required")
    prof=[]
    for s in sites:
        if not s:raise ValueError("Empty publisher station")
        cc=Counter(lifeforms[k] for k in s)
        prof.append({k:cc[k]/len(s) for k in ("T","H","P","G","Ch")})
    tax=[];func=[];stratum_tax={"within":[],"between":[]};stratum_fun={"within":[],"between":[]}
    for i in range(15):
        for j in range(i+1,15):
            u=len(sites[i]|sites[j])
            jac=len(sites[i]&sites[j])/u if u else 0
            over=sum(min(prof[i][k],prof[j][k]) for k in ("T","H","P","G","Ch"))
            group="within" if STATION_MACROCLIMATE[i+1]==STATION_MACROCLIMATE[j+1] else "between"
            tax.append(jac);func.append(over)
            stratum_tax[group].append(jac)
            stratum_fun[group].append(over)
    if len(tax)!=105:raise AssertionError("15 choose 2 must be 105")
    return {
        "tax_jaccard_mean":statistics.mean(tax),
        "lifeform_overlap_mean":statistics.mean(func),
        "tax_jaccard_median":statistics.median(tax),
        "lifeform_overlap_median":statistics.median(func),
        "tax_within_climate_mean":statistics.mean(stratum_tax["within"]),
        "tax_between_climate_mean":statistics.mean(stratum_tax["between"]),
        "lifeform_within_climate_mean":statistics.mean(stratum_fun["within"]),
        "lifeform_between_climate_mean":statistics.mean(stratum_fun["between"]),
        "number_within_climate_pairs":len(stratum_tax["within"]),
        "number_between_climate_pairs":len(stratum_tax["between"]),
    }

STATS=("tax_jaccard_mean","lifeform_overlap_mean",
       "tax_within_climate_mean","tax_between_climate_mean",
       "lifeform_within_climate_mean","lifeform_between_climate_mean")

def summarize(real,samples):
    results={}
    for stat in STATS:
        vals=[d[stat] for d in samples];z=sorted(vals)
        results[stat]={
            "observed":round(real[stat],7),
            "internal_source_null_mean":round(statistics.mean(vals),7),
            "observed_minus_null":round(real[stat]-statistics.mean(vals),7),
            "null_q025":round(z[min(len(z)-1,int(0.025*len(z)))],7),
            "null_q975":round(z[min(len(z)-1,int(0.975*len(z)))],7),
            "technical_lower_tail_fraction_not_field_p":round((1+sum(x<=real[stat] for x in vals))/(1+len(vals)),7),
            "technical_upper_tail_fraction_not_field_p":round((1+sum(x>=real[stat] for x in vals))/(1+len(vals)),7),
            "chain_mean_gap":None,
        }
    return results

def generate(source_sites,lifeforms,climate_stratified,
             chains=N_CHAINS,draws=DRAWS_PER_CHAIN,
             burnin=BURNIN_PROPOSALS,gap=GAP_PROPOSALS):
    if chains<2 or draws<2 or gap<1:raise ValueError("Invalid diagnostic count")
    initial_margin=margins(source_sites,len(lifeforms))
    initial_climate=climate_margins(source_sites,len(lifeforms))
    observed=pairwise_metrics(source_sites,lifeforms)
    samples=[];by_chain=[];accepted=[]
    for chain in range(chains):
        rng=random.Random(SEED+100000*int(climate_stratified)+chain*7193)
        sites=[set(x) for x in source_sites]
        ac=swaps(sites,rng,burnin,climate_stratified)
        this=[]
        for _ in range(draws):
            ac+=swaps(sites,rng,gap,climate_stratified)
            if margins(sites,len(lifeforms))!=initial_margin:
                raise AssertionError("Swaps changed margins")
            if climate_stratified and climate_margins(sites,len(lifeforms))!=initial_climate:
                raise AssertionError("Within-climate swaps changed source local climate occupancy")
            row=pairwise_metrics(sites,lifeforms)
            samples.append(row);this.append(row)
        by_chain.append(this)
        accepted.append(ac)
    summaries=summarize(observed,samples)
    for stat in STATS:
        m=[statistics.mean(v[stat] for v in b) for b in by_chain]
        summaries[stat]["chain_mean_gap"]=round(max(m)-min(m),7)
    return {
        "null":"within_original_source_macroclimate_fixed_margins" if climate_stratified else "whole_original_source_fixed_margins",
        "original_lifeforms":["T","H","P","G","Ch"],
        "n_original_stations":15,
        "n_dependent_pairs":105,
        "original_study_N":1,
        "n_chains":chains,
        "draws_per_chain":draws,
        "draws_total_not_independent_field_replicates":len(samples),
        "burnin_proposals":burnin,
        "spacing_proposals":gap,
        "switches_accepted_by_chain":accepted,
        "observed_descriptive":{k:round(v,7) if isinstance(v,float) else v for k,v in observed.items()},
        "statistics":summaries,
        "source_taxon_occurrence_and_site_richness_FIXED":True,
        "macroclimate_taxon_occurrence_FIXED_if_stratified":climate_stratified,
        "this_is_source_internal_not_independent_nonrail_regional_species_pool":True,
        "cannot_claim_homogenization_relative_to_regional_species_pools":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-pdf",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    row,audit=rows_from_original(args.original_pdf)
    if not audit["full_matrix_verified"]:
        raise ValueError("Original PDF full taxon, station, family, life-form and frequency gates failed")
    sites,life=incidence(row)
    out={
        "version":"v18.5",
        "source":"Jasprica et al. 2017, 10.20302/NC.2017.26.21",
        "source_sha256":audit["publisher_original_pdf_sha256"],
        "source_count":len(row),
        "site_count":len(sites),
        "original_station_climate":STATION_MACROCLIMATE,
        "unrestricted":generate(sites,life,False),
        "within_publisher_climate":generate(sites,life,True),
        "interpretation":"Source-internal and not a comparison to independently surveyed local natural habitat pools. Cannot infer railway-specific global taxonomic or functional convergence, transport, or selection."
    }
    f=Path(args.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(out,ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":
    main()
