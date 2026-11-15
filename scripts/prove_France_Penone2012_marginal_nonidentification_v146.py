#!/usr/bin/env python3
"""v14.6: prove site-identity NON-identifiability from thesis species margins.

The two groups (35/36) are entirely hypothetical mathematical partitions,
NOT the actual French rail lines, station breaks or sampling strata.
All results are conditional on original 'observations' meaning unique sites.
"""
import argparse
import json
from math import comb
from pathlib import Path
from analyze_France_Penone2012_thesis_occupancy_v145 import load_verified_subset

def contrast_for_allocation(allocations, group_a, group_b):
    if min(group_a,group_b)<2:
        raise ValueError("Groups need at least two sites")
    within_pairs=comb(group_a,2)+comb(group_b,2)
    between_pairs=group_a*group_b
    w=b=0
    for name,(a,second) in allocations.items():
        if min(a,second)<0 or a>group_a or second>group_b:
            raise ValueError("Infeasible group allocation: "+name)
        w+=comb(a,2)+comb(second,2)
        b+=a*second
    return w/within_pairs-b/between_pairs

def exact_bounds(frequencies, group_a=35, group_b=36):
    if min(group_a,group_b)<2:
        raise ValueError("Groups must each have two or more sites")
    if not frequencies:
        raise ValueError("No source species")
    # Additive per-taxon contrast: group allocations can be chosen
    # independently when no per-site richness constraints are known.
    minimum_allocation={}
    maximum_allocation={}
    taxon_details=[]
    for name,f in sorted(frequencies.items()):
        if type(f) is not int or f<1 or f>group_a+group_b:
            raise ValueError("Invalid site frequency")
        feasible=[]
        for a in range(max(0,f-group_b),min(group_a,f)+1):
            b=f-a
            per=contrast_for_allocation({name:(a,b)},group_a,group_b)
            feasible.append((per,a,b))
        low_item=min(feasible)
        high_item=max(feasible)
        minimum_allocation[name]=(low_item[1],low_item[2])
        maximum_allocation[name]=(high_item[1],high_item[2])
        taxon_details.append({"species":name,"observations":f,
          "minimum_allocation_A_B":[low_item[1],low_item[2]],
          "maximum_allocation_A_B":[high_item[1],high_item[2]]})
    lower=contrast_for_allocation(minimum_allocation,group_a,group_b)
    upper=contrast_for_allocation(maximum_allocation,group_a,group_b)
    return {
        "version":"v14.6",
        "source":"Penone 2012 thesis Appendix B six audited original observation frequencies",
        "hypothetical_groups_not_actual_railway_lines":True,
        "group_A_n":group_a,"group_B_n":group_b,
        "conditional_original_observations_equal_unique_sites_UNVERIFIED":True,
        "complete_71_site_by_species_matrix_obtained":False,
        "original_six_taxon_frequency_vector":dict(sorted(frequencies.items())),
        "exact_minimum_within_minus_between_mean_six_species_shared_count":
            lower,
        "exact_maximum_within_minus_between_mean_six_species_shared_count":
            upper,
        "both_opposite_signs_feasible_with_identical_six_original_source_margins":
            lower<0<upper,
        "witness_minimum":{k:list(v) for k,v in minimum_allocation.items()},
        "witness_maximum":{k:list(v) for k,v in maximum_allocation.items()},
        "per_taxon_witness":taxon_details,
        "real_railway_connectivity_or_train_seed_dispersal_identified":False,
        "proof_scope":"Exact extrema conditional on six species unique-site frequency "
            "margins and a synthetic 35/36 partition with unconstrained "
            "site richness. Not true French 2 rail lines, no stations, no "
            "geography, no train wind, and not full-species Jaccard.",
    }

def construct_witness_sites(witness, a=35,b=36):
    sites=[set() for _ in range(a+b)]
    for species,(na,nb) in witness.items():
        for idx in range(na):
            sites[idx].add(species)
        for idx in range(nb):
            sites[a+idx].add(species)
    return sites

def direct_site_pair_contrast(sites, a):
    n=len(sites)
    if not 1<a<n-1:
        raise ValueError("Invalid partition")
    wi=[]; be=[]
    for i in range(n):
        for j in range(i+1,n):
            shared=len(sites[i]&sites[j])
            (wi if (i<a)==(j<a) else be).append(shared)
    return sum(wi)/len(wi)-sum(be)/len(be)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--six-csv",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    result=exact_bounds(load_verified_subset(args.six_csv))
    # Constructive witness: compare actual arbitrary binary arrangements.
    for label in ("minimum","maximum"):
        witness=result["witness_"+label]
        simulated=construct_witness_sites(witness)
        direct=direct_site_pair_contrast(simulated,35)
        score=result["exact_"+label+"_within_minus_between_mean_six_species_shared_count"]
        if abs(direct-score)>1e-12:
            raise AssertionError("Combinatorial proof witness and binary matrix mismatch")
    dest=Path(args.output); dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k.startswith("exact_") or
       k.startswith("both_") or k.startswith("real_")},ensure_ascii=False))

if __name__=="__main__":main()
