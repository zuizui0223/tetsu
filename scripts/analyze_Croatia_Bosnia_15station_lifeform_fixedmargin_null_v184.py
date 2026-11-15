#!/usr/bin/env python3
"""v18.4: SOURCE-VERIFIED within-railway community assembly diagnostic.

Jasprica et al. 2017, DOI 10.20302/NC.2017.26.21, 359 taxa x 15 stations.
Nulls: all-site 2x2 fixed-margin checkerboard switches, and switches only
among sites with the SAME original macroclimate. BOTH maintain every station
richness and every source taxon's occurrence frequency; the latter also fixes
within-macroclimate source taxon frequencies.

**NOT** an independent regional NON-rail vegetation species pool. A positive
or negative result cannot identify railway filtering, habitat analog
availability, train seed dispersal, convergence, adaptation, or global
biogeography; site pairs are not independent n=105 replicates.
"""
from __future__ import annotations
import argparse,copy,json,random,statistics,sys
from collections import Counter
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from extract_Croatia_Bosnia_2017_original_rotated_table_matrix_v183 import rows_from_original

SEED=20261115
N_CHAINS=2
DRAWS_PER_CHAIN=120
BURNIN_PROPOSALS=18000
GAP_PROPOSALS=1800
STATION_MACROCLIMATE={
    1:"TSS",2:"TWS",3:"TSS",4:"TWS",5:"TWS",6:"TWS",7:"TWS",8:"TWS",
    9:"TWS",10:"MED",11:"MED",12:"TSS",13:"TSS",14:"TSS",15:"TWS"
}
STATION_NAMES={
    1:"Zadar",2:"Benkovac",3:"Devr ske",4:"Kistanje",5:"Knin",
    6:"Drnis",7:"Unesic",8:"Koprno",9:"Perkovic",10:"Sibenik",
    11:"Split",12:"Ploce",13:"Metkovic",14:"Capljina",15:"Mostar"
}
# Source metadata transcribed independently from Tab. 1, p.274, not modeled climate.
# We do NOT use LAT/LON of individual sensitive plants; these are public stations.
assert len(STATION_MACROCLIMATE)==15 and set(STATION_MACROCLIMATE.values())=={"TWS","TSS","MED"}

def incidence(rows):
    if len(rows)!=359 or len(set(x["taxon_original"] for x in rows))!=359:
        raise ValueError("Use only verified true Jasprica original 359 taxa")
    sets=[set() for _ in range(15)]
    for k,r in enumerate(rows):
        for site in r["source_station_numbers"]:
            if not 1<=site<=15: raise ValueError("Invalid original station")
            sets[site-1].add(k)
    life=[r["life_form"] for r in rows]
    if any(g not in {"T","H","P","G","Ch"} for g in life):
        raise ValueError("Unknown publisher lifeform")
    return sets,life

def margins(sets,n_taxa):
    return (tuple(len(x) for x in sets),
            tuple(sum(k in x for x in sets) for k in range(n_taxa)))

def climate_margins(sets,n_taxa):
    return {g:tuple(sum(k in sets[site-1] for site in range(1,16)
                             if STATION_MACROCLIMATE[site]==g)
                    for k in range(n_taxa))
            for g in ("MED","TWS","TSS")}

def profile_overlap(sets, life):
    """Mean and median sum-min of five lifeform proportions across 105 pairs."""
    traits=sorted(set(life))
    probs=[]
    for site in sets:
        n=len(site)
        if n==0:raise ValueError("Unobserved species cannot be site zero under original")
        c=Counter(life[k] for k in site)
        probs.append({v:c[v]/n for v in traits})
    vals=[sum(min(probs[i][g],probs[j][g]) for g in traits)
          for i in range(15) for j in range(i+1,15)]
    return statistics.mean(vals),statistics.median(vals)

def swaps(sets,rng,n,climate_stratified):
    if climate_stratified:
        groups=[[site-1 for site in range(1,16) if STATION_MACROCLIMATE[site]==g]
                for g in ("MED","TWS","TSS")]
    else:
        groups=[list(range(15))]
    accepted=0
    for _ in range(n):
        indices=rng.choice(groups)
        i,j=rng.sample(indices,2)
        a=sets[i]-sets[j];b=sets[j]-sets[i]
        if not a or not b:continue
        x=rng.choice(tuple(a));y=rng.choice(tuple(b))
        sets[i].remove(x);sets[i].add(y)
        sets[j].remove(y);sets[j].add(x)
        accepted+=1
    return accepted

def null_model(source_sets,life,climate_stratified,chains=N_CHAINS,draws=DRAWS_PER_CHAIN,
               burnin=BURNIN_PROPOSALS,gap=GAP_PROPOSALS):
    original_margin=margins(source_sets,len(life))
    original_climate=climate_margins(source_sets,len(life))
    observed,_=profile_overlap(source_sets,life)
    samples=[];acceptances=[];chain_means=[]
    for ch in range(chains):
        rng=random.Random(SEED+100000*int(climate_stratified)+ch*7193)
        sets=[set(x) for x in source_sets]
        accepted=swaps(sets,rng,burnin,climate_stratified)
        one=[]
        for _ in range(draws):
            accepted+=swaps(sets,rng,gap,climate_stratified)
            if margins(sets,len(life)) != original_margin:
                raise AssertionError("Swap altered biological site or taxon margins")
            if climate_stratified and climate_margins(sets,len(life)) != original_climate:
                raise AssertionError("Climate-constrained source occurrences were altered")
            val,_=profile_overlap(sets,life)
            one.append(val);samples.append(val)
        chain_means.append(statistics.mean(one))
        acceptances.append({"chain":ch,"accepted_switches":accepted,
                            "proposal_count":burnin+draws*gap,
                            "acceptance_rate":round(accepted/(burnin+draws*gap),4)})
    out={
        "null":"same_source_climate_fixed_margins" if climate_stratified else "all_stations_fixed_margins",
        "number_climate_groups":3 if climate_stratified else 1,
        "n_chains":chains,"draws_per_chain":draws,
        "n_synthetic_replicates":len(samples),
        "seed":SEED,
        "burnin_proposals_per_chain":burnin,
        "proposals_between_replicates":gap,
        "observed_mean_site_pair_lifeform_profile_overlap":round(observed,6),
        "null_mean":round(statistics.mean(samples),6),
        "null_sd":round(statistics.stdev(samples),6) if len(samples)>1 else None,
        "null_q025":round(sorted(samples)[max(0,int(0.025*len(samples)))],6),
        "null_q975":round(sorted(samples)[min(len(samples)-1,int(0.975*len(samples)))],6),
        "observed_minus_null_mean":round(observed-statistics.mean(samples),6),
        "one_sided_tail_ge_observed":round((1+sum(x>=observed for x in samples))/(1+len(samples)),6),
        "one_sided_tail_le_observed":round((1+sum(x<=observed for x in samples))/(1+len(samples)),6),
        "chain_mean_values":[round(x,6) for x in chain_means],
        "chain_mean_max_minus_min":round(max(chain_means)-min(chain_means),6),
        "acceptance":acceptances,
        "pseudo_independent_station_pairs_not_samples":True,
        "technical_permutation_p_exploratory_samples_autocorrelated":True,
        "null_is_internal_rail_species_pool_NOT_natural_habitat_pool":True}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-pdf",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    rows,audit=rows_from_original(a.original_pdf)
    if not audit["full_matrix_verified"]:
        raise ValueError("Requires source-original PDF table2 exact 359×15 publisher audit")
    sites,life=incidence(rows)
    raw_count,raw_col=margins(sites,len(life))
    observed_mean,observed_median=profile_overlap(sites,life)
    data={
        "source_original_DOI":"10.20302/NC.2017.26.21",
        "source_original_pdf_sha256":audit["publisher_original_pdf_sha256"],
        "source_publisher_taxa":359,"source_publisher_stations":15,
        "field_one_study_across_two_countries":True,
        "source_survey_2013_2014_sampling_area_varies_200_to_800m2":True,
        "source_spring_only_surveys_absent_in_4_southern_stations":True,
        "source_publisher_original_table_anomaly_Canna":audit["source_frequency_anomaly_is_exactly_publisher_Canna_indica"],
        "source_publisher_singleton_table_vs_prose_mismatch":audit["publisher_table_singletons_prose_discrepancy"],
        "source_station_original_richness":dict(enumerate(raw_count,1)),
        "source_taxon_occupancy_frequency_sum":sum(raw_col),
        "source_station_macroclimate_from_publisher_table1":STATION_MACROCLIMATE,
        "observed_mean_pairwise_lifeform_overlap":round(observed_mean,6),
        "observed_median_pairwise_lifeform_overlap":round(observed_median,6),
        "null_1":null_model(sites,life,False),
        "null_2":null_model(sites,life,True),
        "CAUTION":"Internal railway species-pool fixed-margin nulls NOT regional non-rail species pool; no rail causal homogenization or selection.",
    }
    f=Path(a.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(data,ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":
    main()
