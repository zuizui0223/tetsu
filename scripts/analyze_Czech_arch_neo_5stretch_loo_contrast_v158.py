#!/usr/bin/env python3
"""v15.8: posthoc paired ARCHAEOPHYTE minus NEOPHYTE same-stretch identity excess.

Author-original Czech botanical records: Kutlvašr et al. (2024),
Preslia DOI 10.23855/preslia.2024.247. Use 35/39 real original station
floras exactly matched to independent 2026 SR70 station centroid proxies.
Subgroups are explicitly author tagged archaeophyte ("arch") and
neophyte ("neo") taxa; untagged plants NEVER declared native.
Previous v15.7 estimated the subgroups separately. This study
POSTHOC contrasts the two groups ON THE SAME ORIGINAL PAIR WINDOW
and checks all five original route-stretch leave-one-out subsets.
The groups have different prevalence, richness and arrival histories:
difference of their raw Jaccard gaps is NOT a dispersal-rate difference.
No independent train-seed carriage or mowing/wind observed.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from statistics import mean
from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import (
    original_35_station_sites, real_site_pairs, real_site_pairs_without_n_gate,
    conditional_source_gap
)
from analyze_Czech_2024_real_archaeophyte_neophyte_geographic_identity_v157 import (
    read_source_status, group_sets, source_class_pairs, source_subgroup_gap
)

MODES=("geography_only","geography_and_function","geography_function_and_area")
MIN_STRETCHES=5
MIN_SITES=30

def paired_group_stats(sites,guild,status,methods=MODES):
    if len({x["section"] for x in sites})<4:
        return {"identified":False,"reason":"Too few original Czech railway route stretches"}
    source = real_site_pairs(sites,guild) if len(sites)==35 else real_site_pairs_without_n_gate(sites,guild)
    classes=group_sets(sites,status)
    if not all(classes[x]["n_actual_original_site_stations_with_class"]==len(sites)
               for x in ("archaeophyte","neophyte")):
        return {"identified":False,
                "reason":"At least one original source station lacks explicit arch or neo taxon; do not compare double-empty Jaccard"}
    results={}
    for label in ("archaeophyte","neophyte"):
        d=classes[label]
        sets=d["source_site_sets"]
        rows=source_class_pairs(source,sets)
        results[label]={
          "source_original_taxa":d["n_source_distinct_observed_class_taxa"],
          "original_station_taxon_incidences":d["n_source_site_taxon_incidence"],
          "mean_plant_taxa_per_study_station":d["n_source_site_taxon_incidence"]/len(sites),
          "raw_source_gap":source_subgroup_gap(rows,sets)["same_minus_other"],
          "conditioned":{
            mode:conditional_source_gap(rows,sets=sets,mode=mode)
            for mode in methods
          }
        }
    paired={}
    for mode in methods:
        a=results["archaeophyte"]["conditioned"][mode]
        n=results["neophyte"]["conditioned"][mode]
        if not a["identified"] or not n["identified"]:
            paired[mode]={"identified":False,
                          "reason":"At least one source group fails common geography/guild/area support"}
            continue
        # Both have exactly same pair IDs and same base original source
        # functional/geographical (and optionally area) strata, because
        # both class sets appear in all study sites.
        if a["source_support"]!=n["source_support"]:
            raise ValueError("Subgroups have inconsistent original matched pair support")
        paired[mode]={
          "identified":True,
          "archaeophyte_original_matched_source_gap":a["weighted_same_minus_different_jaccard"],
          "neophyte_original_matched_source_gap":n["weighted_same_minus_different_jaccard"],
          "arch_minus_neo_source_gap_contrast":(
              a["weighted_same_minus_different_jaccard"]-
              n["weighted_same_minus_different_jaccard"]),
          "identical_original_source_pair_and_distance_guild_strata":True,
          "shared_original_pair_support":a["source_support"],
        }
    return {
      "identified":True,
      "n_original_source_stations":len(sites),
      "n_original_source_stretches":len(set(s["section"] for s in sites)),
      "n_dependent_source_station_pairs":len(source),
      "subgroup":results,
      "raw_arch_minus_neo_gap":(
          results["archaeophyte"]["raw_source_gap"]-
          results["neophyte"]["raw_source_gap"]),
      "matched_group_contrast":paired,
      "not_a_paired_taxon_experiment":True
    }

def run(s1,s2,disp,life,sr70):
    sites,guild,src=original_35_station_sites(s1,s2,disp,life,sr70)
    status,catprov=read_source_status(s1)
    if len(sites)!=35 or len({s["section"] for s in sites})!=5:
        raise ValueError("Original 35 properly matched Czech stations from all 5 study stretches needed")
    main=paired_group_stats(sites,guild,status)
    if not main["identified"]:
        raise ValueError("Original two source class contrasts not identifiable")
    out=[]
    for stretch in sorted({s["section"] for s in sites}):
        subset=[x for x in sites if x["section"]!=stretch]
        x=paired_group_stats(subset,guild,status)
        out.append({
          "source_stretch_excluded":stretch,
          "n_original_remaining_stations":len(subset),
          "n_independent_countries_still_just_one":1,
          "both_groups_common_support_identified":x["identified"] and x["matched_group_contrast"]["geography_and_function"]["identified"],
          "matched_arch_minus_neo_original_source_gap":(
              x["matched_group_contrast"]["geography_and_function"]["arch_minus_neo_source_gap_contrast"]
              if x["identified"] and x["matched_group_contrast"]["geography_and_function"]["identified"] else None),
          "separate_group_covered_original_pair_counts": (
              x["matched_group_contrast"]["geography_and_function"]["shared_original_pair_support"]
              if x["identified"] and x["matched_group_contrast"]["geography_and_function"]["identified"] else None),
        })
    values=[x["matched_arch_minus_neo_original_source_gap"] for x in out
            if x["matched_arch_minus_neo_original_source_gap"] is not None]
    return {
      "version":"v15.8_POSTHOC_author_original_Czech_arch_minus_neo_matched_shared_source",
      "source_preslia_DOI":"10.23855/preslia.2024.247",
      "publisher_and_independent_station_geometry_provenance":src,
      "explicit_original_source_taxon_category_provenance":catprov,
      "fixed_same_pair_window_arch_vs_neo_contrast":main,
      "all_original_stretch_leave_one_out_posthoc_sensitivity":{
        "n_original_stretches":5,"n_estimable":len(values),
        "n_positive_arch_minus_neo":sum(v>0 for v in values),
        "median_posthoc_loo_arch_minus_neo":sorted(values)[len(values)//2] if values else None,
        "minimum":min(values) if values else None,
        "maximum":max(values) if values else None,
        "individual_original_stretch_leaveouts":out,
      },
      "limits":[
        "POSTHOC comparison after v15.7 subgroup gaps were seen, not a preregistered causal test.",
        "No true author 2020-22 plot GPS: 2026 SR70 railway station centroids are proxies.",
        "5 route stretches/35 railway stations in only one original Czech field study; 595 station pairs share individual station endpoints.",
        "Archaeophyte versus neophyte original taxon prevalence and per-station richness differ; difference in Jaccard is NOT a difference in per-seed railway dispersal probabilities.",
        "Original categories are author S1 explicitly arch/neo tagged only. Source unlabelled plants never imputed as natives.",
        "The same source geographical/functional matching bins inherited from Slovenia/Czech v15.4 are held fixed; no new adaptive bin tuning.",
        "One original field matrix cannot identify freight train operational route, on-train propagules, grass-cutting pressure, wind-adaptation or genetic evolution.",
      ],
      "no_railway_transport_or_evolution_causal_claim":True
    }

def main():
    p=argparse.ArgumentParser()
    for n in ("s1","s2","pladias-diaspore","pladias-lifeform","sr70","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    result=run(a.s1,a.s2,a.pladias_diaspore,a.pladias_lifeform,a.sr70)
    o=Path(a.output);o.parent.mkdir(parents=True,exist_ok=True)
    o.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    main=result["fixed_same_pair_window_arch_vs_neo_contrast"]
    print(json.dumps({
       "version":result["version"],
       "original_source_subgroup_archaeophyte":main["subgroup"]["archaeophyte"],
       "original_source_subgroup_neophyte":main["subgroup"]["neophyte"],
       "source_original_matched_group_contrasts":main["matched_group_contrast"],
       "source_leave_stretch_out":result["all_original_stretch_leave_one_out_posthoc_sensitivity"],
       "posthoc_NOT_randomized_causal_train_test":True,
    },ensure_ascii=False))

if __name__=="__main__":main()
