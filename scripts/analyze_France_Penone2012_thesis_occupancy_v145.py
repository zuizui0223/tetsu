#!/usr/bin/env python3
"""Source-audited six Penone 2012 thesis Appendix B observations.

Conditional combinatorial moments only: author Appendix B's
'Number of observations*' has not independently been proven to mean
distinct presences at 71 unique sites. No site incidence matrix obtained.
"""
import argparse
import csv
import json
from math import comb
from pathlib import Path

N = 71
EXACT = {
    "Arrhenatherum elatius": (69,175),
    "Brassica napus": (3,175),
    "Cirsium arvense": (17,175),
    "Papaver rhoeas": (39,177),
    "Rubus fruticosus": (53,178),
    "Senecio inaequidens": (39,178),
}

def load_verified_subset(filename):
    with open(filename, encoding="utf-8-sig", newline="") as h:
        r = csv.DictReader(h)
        if r.fieldnames != ["species","observations","page"]:
            raise ValueError("Penone original source transcription schema drift")
        observed = {}
        for row in r:
            name = row["species"]
            if name in observed or name not in EXACT:
                raise ValueError("Unverified or duplicated author taxon")
            source = (int(row["observations"]),int(row["page"]))
            if source != EXACT[name]:
                raise ValueError("Original species observation count/page mismatch")
            observed[name] = source[0]
    if set(observed) != set(EXACT):
        raise ValueError("Missing original verified author taxon")
    return observed

def marginal_moments(frequencies, n_sites=N):
    if type(n_sites) is not int or n_sites < 2 or not frequencies:
        raise ValueError("Invalid independent site count or empty source frequency")
    if any(type(v) is not int or v < 1 or v > n_sites for v in frequencies.values()):
        raise ValueError("An observation count is incompatible with site incidence")
    n_pairs = comb(n_sites,2)
    incidences = sum(frequencies.values())
    pair_intersections = sum(comb(v,2) for v in frequencies.values())
    return {
        "n_sites_conditionally_assumed":n_sites,
        "species_in_verified_subset":len(frequencies),
        "source_observations_sum_six":incidences,
        "unordered_site_pairs":n_pairs,
        "six_species_shared_presence_pair_total_conditional":pair_intersections,
        "six_species_contribution_mean_richness_conditional":incidences/n_sites,
        "six_species_contribution_mean_intersection_conditional":pair_intersections/n_pairs,
    }

def analyze(frequencies):
    result = marginal_moments(frequencies)
    result.update({
        "version":"v14.5",
        "original_primary_source":"Penone 2012 thesis Appendix B pp.175-178",
        "original_thesis_url":"https://www.ittecop.fr/doc/0THESES/PENONE.pdf",
        "same_original_field_cohort_paper_doi":"10.1016/j.biocon.2012.01.041",
        "paper_published_total_species":186,
        "thesis_appendix_C_all_species_N":187,
        "discrepancy_186_vs_187_unresolved":True,
        "appendix_B_label":"Number of observations*",
        "distinct_site_observation_grain_verified":False,
        "conditional_71_site_incidence_assumption":True,
        "complete_71_site_by_species_matrix":False,
        "can_calculate_distance_matched_Jaccard_or_route_break_effect":False,
        "can_infer_train_seed_dispersal_or_stress_adaptation":False,
        "added_independent_cohort_to_17_original_study_inventory":False,
        "original_taxon_observations":[
            {"species":key,"observations":value,"printed_page":EXACT[key][1],
             "fraction_of_distinct_site_pairs_both_present_conditional":
                comb(value,2)/comb(N,2)}
            for key,value in sorted(frequencies.items())
        ],
        "inference_limit":"Conditional mathematical six-taxon moment, NOT full flora "
            "pairwise Jaccard or source-confirmed site occupancy. The original "
            "author's 'Number of observations*' unit must be reconciled to "
            "unique sites before calling it a field lower bound. No spatial "
            "rail effect, train wind, mowing or propagule fitness identified.",
    })
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-csv",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    out=analyze(load_verified_subset(a.source_csv))
    p=Path(a.output)
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False))
