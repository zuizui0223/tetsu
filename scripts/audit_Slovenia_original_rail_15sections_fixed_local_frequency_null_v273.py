#!/usr/bin/env python3
"""v27.3 original Slovenia railway source BS/BT fixed species-frequency null.

True author-original Šabić et al 2026 CC BY supplement, publisher XLSX
401279 bytes SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f.
89 actual locations / 173 genuine same-side ballast SHOULDER (BS) versus
ballast TOE (BT) plot pairs nested in 15 ORIGINAL railway line sections.

Pool L/R matched side species sets WITHIN each genuine source location
before analyzing species identities; these are not independent new sites.
BS observed species sets FROZEN. BT source plot incidence checkerboard
swaps WITHIN each original railway section preserve (i) each source
original BT location's species richness and (ii) each source taxon's BT
occupancy within source original section EXACTLY. Source 14 sections with
>=2 sites/section (88 sites) form first comparator; stricter >=4 locations
per section sensitivity. Original 1-site source section not randomized.

Original v25.7 reference-label derangement on 14 sections gave
same-site Jaccard similarity 0.3632363548846349, genuine source verified.
Require exact reproduction before reporting source-fixed-margins result.

This isn't non-rail control: BOTH source microhabitats are railway
ballast zones, and BS vs BT surveyed areas cannot be assumed identical.
No traffic counts, train-seed carriage, genetic changes, wind-pressure
measurements or mowing dates. MCMC source ranks are technical diagnostics,
NOT randomized-field p-values. No original plant-by-site/EPPO/GPS exported.
"""
from __future__ import annotations
import argparse,json
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_rail_local_site_coupling_v257 import slovenia_location_sets
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    sampler,finish_null,METRICS
)

DRAWS=99
BURNIN=3000
SPACING=350
SEED=273118
EXPECTED_SIM_V257=.3632363548846349
SOURCE_MIN2={"sections":14,"sites":88}
MIN_SECTIONS_SENSITIVITY=4


def groups_from_original(sites:dict, min_sites:int)->dict:
    if min_sites<2:raise ValueError("Must preserve >1 true source site per original section")
    grouped=defaultdict(list)
    for location,entry in sites.items():
        if not entry["BS"] or not entry["BT"] or entry["matched_sides"] not in (1,2):
            raise ValueError("Original source shoulder/toe site pairing incomplete")
        grouped[entry["section"]].append({
            "original_location":location,
            "BS":set(entry["BS"]),"BT":set(entry["BT"]),
        })
    return {section:sorted(v,key=lambda r:str(r["original_location"]))
            for section,v in sorted(grouped.items()) if len(v)>=min_sites}


def analyze_group(grouped:dict,draws:int=DRAWS,
                  burnin:int=BURNIN,spacing:int=SPACING,
                  seed:int=SEED)->dict:
    source={}
    for idx,(section,items) in enumerate(sorted(grouped.items())):
        out=sampler([r["BS"] for r in items],[r["BT"] for r in items],
                    draws=draws,burnin=burnin,spacing=spacing,
                    seed=seed+1481*idx)
        summary=finish_null(out["observed"],out["null_draws"])
        source[section]={
            "source_original_locations":len(items),
            "source_original_BT_species_pool":out["original_reference_species_pool"],
            "source_accepted_checkerboard_swaps":out["accepted_swaps_total"],
            "original_site_and_species_margins_preserved_exact":True,
            "mean_source_section_frequency_fixed_null":summary,
            "technical_null_draws":out["null_draws"]
        }
    if not source:
        return {"identified":False,"source_rail_sections":0}
    raw={metric:mean(x["mean_source_section_frequency_fixed_null"][metric]["observed"] for x in source.values())
         for metric in METRICS}
    null=[{metric:mean(x["technical_null_draws"][k][metric] for x in source.values())
           for metric in METRICS} for k in range(draws)]
    overall=finish_null(raw,null)
    for row in source.values():
        del row["technical_null_draws"]
    return {
        "identified":True,
        "source_original_rail_sections":len(source),
        "source_original_spatial_locations":sum(x["source_original_locations"] for x in source.values()),
        "n_source_sections_with_no_accepted_swaps":sum(x["source_accepted_checkerboard_swaps"]==0 for x in source.values()),
        "source_equal_section_global_metric":overall,
        "n_original_sections_observed_jaccard_below_frequency_null":sum(
            x["mean_source_section_frequency_fixed_null"]["beta_jac"]["observed_minus_frequency_null"]<0
            for x in source.values()),
        "n_original_sections_observed_shared_taxa_above_frequency_null":sum(
            x["mean_source_section_frequency_fixed_null"]["shared_taxa"]["observed_minus_frequency_null"]>0
            for x in source.values()),
        "original_author_section_aggregate_support":source,
        "site_species_occupancy_fixed_inside_original_section":True,
        "source_sections_not_independent_countries_or_random_field_trials":True
    }


def run(xlsx:Path,audit:Path,draws:int=DRAWS,
        burnin:int=BURNIN,spacing:int=SPACING)->dict:
    original=extract_original_xlsx(xlsx,audit)
    sites,source=slovenia_location_sets(original)
    panel2=groups_from_original(sites,2)
    qa={"sections":len(panel2),"sites":sum(len(v) for v in panel2.values())}
    if qa!=SOURCE_MIN2:
        raise ValueError(f"Frozen 14 genuine original sections with 88 locations changed: {qa}")
    main=analyze_group(panel2,draws,burnin,spacing,SEED)
    observed_similarity=1.-main["source_equal_section_global_metric"]["beta_jac"]["observed"]
    if abs(observed_similarity-EXPECTED_SIM_V257)>1e-10:
        raise ValueError(f"Earlier independent genuine source v25.7 similarity not reproduced: {observed_similarity:.12f} != {EXPECTED_SIM_V257:.12f}")
    panel4=groups_from_original(sites,4)
    narrow=analyze_group(panel4,draws,burnin,spacing,SEED+47001)
    return {
        "version":"v27.3",
        "genuine_original_railway_study_DOI":"10.1007/s44473-026-00155-z",
        "source_original_author_QA":source,
        "source_original_study_sites":89,
        "source_original_study_matched_sides":173,
        "source_original_study_rail_sections":15,
        "source_14_multisite_88location_panel":main,
        "source_four_original_locations_per_section_sensitivity":narrow,
        "previous_v25_7_source_same_site_jaccard_similarity_reproduced_exact_1e_minus_10":True,
        "original_source_BS_BT_both_true_RAIL_microhabitats":True,
        "original_BS_BT_field_survey_area_equality_UNVERIFIED":True,
        "original_BS_BT_pairwise_comparison_not_nonrail_reference":True,
        "MCMC_draws_not_independent_species_or_field_experiments":True,
        "no_train_movement_transport_mowing_wind_or_genetic_adaptation_identified":True,
        "no_species_site_original_rows_or_locality_export":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-xlsx",type=Path,required=True)
    p.add_argument("--source-audit",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    a=p.parse_args()
    out=run(a.source_xlsx,a.source_audit,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for label,key in (("min2_sections","source_14_multisite_88location_panel"),
                      ("min4_sections","source_four_original_locations_per_section_sensitivity")):
        s=out[key]
        print(json.dumps({
            "version":"v27.3","source_group":label,
            "identified":s["identified"],
            "sections":s.get("source_original_rail_sections"),
            "sites":s.get("source_original_spatial_locations"),
            "zero_acceptance_original_sections":s.get("n_source_sections_with_no_accepted_swaps"),
            "beta":s.get("source_equal_section_global_metric",{}).get("beta_jac"),
            "shared_taxa":s.get("source_equal_section_global_metric",{}).get("shared_taxa"),
            "n_sections_with_observed_below_freqnull":
                s.get("n_original_sections_observed_jaccard_below_frequency_null"),
        },ensure_ascii=False))


if __name__=="__main__":main()
