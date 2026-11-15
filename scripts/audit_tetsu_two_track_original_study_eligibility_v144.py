#!/usr/bin/env python3
"""v14.4: source-backed two-track railway macro study evidence eligibility.

This script is NOT a meta-analysis: each row describes a previously
published independent DOCUMENT (some papers use overlapping plant cohorts).
Stages are intentionally not numerically pooled. Only bona fide
field-site-by-species matrices can support common primary B Jaccard
reanalysis. Unretrieved French/Indian/Japanese raw plots are not
materialized, missing vegetation species never assumed absent.

Same site but overlapping years/locations in Swiss official monitoring
and published port/rail genetic papers are explicitly tagged.
The Mexico/Slovenia/Czech country list is NOT a global representative
random sample; one country is not one independent randomized train study.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
FIELDS=("study_key","source_identifier","country","track","observation_stage",
        "original_independent_units","year_range","original_field_unit",
        "true_site_taxon_matrix_obtained","true_site_coordinates_verified",
        "usable_for_common_species_beta","source_original_or_authorship",
        "exposure_or_comparator","principal_limit")
ELIGIBLE_TRACKS={"A","A_ROAD_CONTROL","A_GENETIC_TRIANGULATION",
                 "B","B_CANDIDATE","B_CONTEXT","B_ROAD_CONTROL"}
REAL_B_SITE_LIST={"SI_Sabic_2026","CZ_Kutlvasr_2024"}
SAME_GEO_REGION_OVERLAPS=(
  ("CH_Schoenenberger_2012","CH_Schulze_2014_port","CH_BaselCanton_2020"),
  ("FR_Penone_2012","FR_Blanchet_2015"),
)


def inventory(csvpath):
    with open(csvpath,encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f)
        if tuple(rd.fieldnames or ())!=FIELDS:
            raise ValueError("Two-track 2026 v14.4 evidence source schema drift")
        studies=list(rd)
    if len(studies)!=17 or len({r["study_key"] for r in studies})!=17:
        raise ValueError("Original v14.4 source evidence 17 publication records not preserved or duplicated")
    for r in studies:
        if (r["track"] not in ELIGIBLE_TRACKS or
            not r["source_identifier"] or not r["original_field_unit"]
            or not r["principal_limit"]):
            raise ValueError("Unspecified source study type or original sampling unit")
        if any(r[k] not in ("true","false")
               for k in ("true_site_taxon_matrix_obtained",
                         "true_site_coordinates_verified",
                         "usable_for_common_species_beta")):
            raise ValueError("Source true/false missing source availability")
        if r["usable_for_common_species_beta"]=="true":
            if r["track"]!="B" or r["true_site_taxon_matrix_obtained"]!="true":
                raise ValueError("Cannot qualify meta community analysis without authentic site species incidence")
        if int(r["original_independent_units"])<0:
            raise ValueError("Source authors independent unit count cannot be negative")
    seen={r["study_key"]:r for r in studies}
    if {k for k,r in seen.items() if r["track"]=="B"}!=REAL_B_SITE_LIST:
        raise ValueError("The only two currently verified primary field site-flora matrices are Slovenia/Czech")
    if any(seen[k]["usable_for_common_species_beta"]!="true"
           for k in REAL_B_SITE_LIST):
        raise ValueError("Primary original Slovenia/Czech source incidence lost")
    for group in SAME_GEO_REGION_OVERLAPS:
        if not all(k in seen for k in group):
            raise ValueError("Published overlapping Swiss/French source cohorts not logged")
    if seen["FR_Penone_2012"]["true_site_taxon_matrix_obtained"]!="false":
        raise ValueError("Do not claim 2012 original 71-site France flora raw matrix was downloaded")
    if seen["SI_Sabic_2026"]["true_site_coordinates_verified"]!="true":
        raise ValueError("Source 89 Slovenian actual start/end WGS84 GPS verification missing")
    if seen["CZ_Kutlvasr_2024"]["true_site_coordinates_verified"]!="false":
        raise ValueError("No verified original Czech 39 station-level geography in S1/S2")
    return studies


def result(studies):
    track=Counter(r["track"] for r in studies)
    country=defaultdict(set)
    for r in studies:country[r["track"]].add(r["country"])
    primaryB=[r for r in studies if r["usable_for_common_species_beta"]=="true"]
    stage=Counter(r["observation_stage"] for r in studies if r["track"].startswith("A"))
    return {
        "complete":True,"version":"v14.4",
        "literature_inventory_status":"TARGETED verified seed and railway floristic source evidence map, NOT systematic full PRISMA search or effect size meta",
        "n_unique_published_source_study_documents":len(studies),
        "n_original_studies_by_pipeline_track":dict(track),
        "n_countries_by_track_including_road_controls":{k:sorted(v) for k,v in country.items()},
        "stage_A_source_evidence_category_counts":dict(stage),
        "primary_B_source_true_plot_by_taxon_matrices_obtained":
            [r["study_key"] for r in primaryB],
        "source_2012_France_real_71site_break_vegetation_raw_matrix_acquired":False,
        "SAME_geography_or_overlapping_source_cohorts_MUST_not_count_as_independent_experiments":
            [list(g) for g in SAME_GEO_REGION_OVERLAPS],
        "n_effect_sizes_meta_harmonized_for_rail_vs_nonrail_or_rail_vs_road_with_same_trait_method":0,
        "n_actual_direct_rail_freight_cargo_grain_x_same_site_viable_germination_x_adult_matched_cohorts":0,
        "n_independent_countries_with_TRUE_ready_B_incidence_matrix":2,
        "no_cross_country_rail_dispersal_meta_effect_claimed":True,
        "next_source_gates":[
          "French Penone 2012 original 71-site by 186-species incidence and original railway break/geographic source exposure needed, data not yet acquired",
          "original French RER 394 equal-area plot data by author request",
          "author-verified original Czech station GPS and road/urban geographic source pool",
          "Original Slovenian same-country source non-rail controls and active freight transport OD needed",
          "Stage A viability and recruitment in SAME species/time/site cohorts, not product of unrelated marginal published estimates"
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-source-inventory",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=result(inventory(a.original_source_inventory))
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
