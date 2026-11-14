#!/usr/bin/env python3
"""v11.6: original Slovenia RAIL dispersal-strategy placebo panel.

The 2026 source study codes dispersal.strategy values by original
named exemplar strategies (Epilobium, Phragmites, Allium, Bidens, etc.).
Question: Is frequency-sensitive BS−BT habitat preference a special
property of source wind dispersal classes or do OTHER source dispersal
categories yield a similar weight-scheme gap?

All original source 89 matched BS/BT sites/15 sections/261 known
life-history+strategy taxa from Šabić et al. 2026
DOI 10.1007/s44473-026-00155-z (published external field source).
Reuses exact site-side/section tests from verified v11.5.
Runs ONLY groups with >=3 source labeled taxa in BOTH strict-annual and
other-known-lifespan strata; otherwise reports NOT IDENTIFIABLE.
Tests every original strategy and two overlap-defined wind guilds;
no favorable subset hand picked.
Bonferroni is an exploratory family-wise statistical adjustment,
not a causal randomized ecological inference. Site areas unequal.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_DOI,SOURCE_SHA256,extract_original_xlsx,
)
from analyze_Slovenia_rail_species_occupancy_permutation_v109 import (
    STRATA,WIND,original_location_matched_sets,
)
from analyze_Slovenia_rail_side_aggregation_section_flip_v115 import (
    extract_side_pairs,preference_site_matrix,crossing,exact_section_label_swaps,
)

SOURCE_CATEGORIES=("Allium","Epilobium","Cornus","Bidens",
                   "Sparganium","Phragmites","Lycopodium","Zea")
COMPARISONS={
    "strict_wind_Epilobium_Phragmites":("Epilobium","Phragmites"),
    "broad_wind_plus_Lycopodium_spore":("Epilobium","Phragmites","Lycopodium"),
    **{"source_strategy_"+c:(c,) for c in SOURCE_CATEGORIES},
}


def label_strategy(rows,traits,selected_categories):
    if not set(selected_categories).issubset(SOURCE_CATEGORIES):
        raise ValueError("A specified dispersal type is not in original study source codebook")
    alt=[]
    for original in rows:
        t=dict(original)
        code=original["taxon_code"]
        value=traits[code].get("Diaspore.strategy")
        if value not in SOURCE_CATEGORIES:
            raise ValueError("Missing or unexpected original strategy in previously source-verified trait-known cohort")
        t["wind"]=value in selected_categories
        alt.append(t)
    return alt


def support(rows):
    bins={s:{"focus":0,"other":0} for s in STRATA}
    for r in rows:
        bins[r["stratum"]]["focus" if r["wind"] else "other"]+=1
    eligible=all(min(bins[s].values())>=3 for s in STRATA)
    return bins,eligible


def run(original,prev115):
    if (prev115.get("complete") is not True or prev115.get("version")!="v11.5"
        or prev115.get("original_source_SHA256")!=SOURCE_SHA256
        or prev115.get("n_original_true_RAIL_field_sites")!=89):
        raise ValueError("Require authenticated original v11.5 paired-side site source")
    sites,traits,integrity=original_location_matched_sets(original)
    side,trait2,_=extract_side_pairs(original)
    if len(sites)!=89 or len(side)!=89 or traits!=trait2:
        raise ValueError("Source taxon and true field location data drifted")
    rows,section=preference_site_matrix(side,traits,mode="union")
    if len(rows)!=261 or len(set(section.values()))!=15:
        raise ValueError("Source 261 taxa/15 original rail sections not recovered")
    outputs={}
    eligible=[]
    for label,cats in COMPARISONS.items():
        classified=label_strategy(rows,traits,cats)
        n,ok=support(classified)
        result={
            "defined_original_diaspore_strategy_codes":list(cats),
            "counts_by_original_lifespan":n,
            "estimable_two_strata_min_three_source_taxa":ok,
        }
        if ok:
            curve=crossing(classified)
            shuffle=exact_section_label_swaps(classified,section)
            result.update({
                "species_equal_original_source":curve["species_equal"],
                "occupancy_weighted_original_source":curve["site_opportunity_weighted"],
                "weight_scheme_gap":curve["species_equal"]-curve["site_opportunity_weighted"],
                "source_frequency_exponent_zero_crossing":curve["source_frequency_exponent_zero_crossing_0_to_1"],
                "source_opposite_signs":curve["species_equal"]*curve["site_opportunity_weighted"]<0,
                "original_exact_section_label_flip_p_gap":
                    shuffle["two_sided_full_15section_flip_p_weighting_gap"],
                "original_exact_section_label_flip_p_equal":
                    shuffle["two_sided_full_15section_flip_p_species_equal"],
                "original_exact_section_label_flip_p_weighted":
                    shuffle["two_sided_full_15section_flip_p_occ_weighted"],
                "n_original_full_section_swaps":shuffle["n_possible_section_exchangeability_flips"],
            })
            eligible.append(label)
        outputs[label]=result
    for label in eligible:
        p=outputs[label]["original_exact_section_label_flip_p_gap"]
        outputs[label]["exploratory_Bonferroni_across_all_estimable_strategy_gaps"]=min(1.0,p*len(eligible))
        outputs[label]["n_original_strategy_contrasts_in_multiplicity_family"]=len(eligible)
    strict=outputs["strict_wind_Epilobium_Phragmites"]
    if not strict["estimable_two_strata_min_three_source_taxa"]:
        raise ValueError("Source strict-wind minimum diagnostic support unexpectedly missing")
    orig=prev115["exact_32768_rail_section_BS_BT_label_swaps_original_union"]
    if not math.isclose(strict["weight_scheme_gap"],orig["original_gap_equal_minus_occupancy"],abs_tol=1e-12):
        raise ValueError("Original source strict-wind 261 taxon primary gap not reproduced")
    return {
        "complete":True,"version":"v11.6",
        "original_source_DOI":SOURCE_DOI,"original_supplement_SHA256":SOURCE_SHA256,
        "n_original_RAIL_field_locations":89,"n_original_RAIL_sections":15,
        "n_original_strict_same_side_BS_BT_pairs":173,
        "n_source_taxa_with_original_known_strategy_and_lifespan":261,
        "original_all_source_diaspore_strategies_included":list(SOURCE_CATEGORIES),
        "n_fully_identified_source_trait_group_vs_rest_comparisons":len(eligible),
        "all_original_strategy_controls_including_missing_strata":outputs,
        "focus":strict,
        "limits":[
            "Exploratory multiple original named categories share species and are not independent biological treatment arms.",
            "Plant-specific strategy labels are fixed traits and may be phylogenetically autocorrelated.",
            "Exact section label flips are not an experimental randomization; track shoulder/toe plot geometry may be different.",
            "Source plant occupancy is not seed influx, functional fitness, phenotypic plasticity, train wind or maintenance intensity.",
            "No original biological data outside Slovenia added in v11.6."
        ],
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("original-xlsx","source-audit102","verified115","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    data=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    parent=json.loads(Path(a.verified115).read_text(encoding="utf-8"))
    out=run(data,parent)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,"version":"v11.6",
        "n_testable_original_classes":out["n_fully_identified_source_trait_group_vs_rest_comparisons"],
        "all_source_strategy_results":{
            name:{"counts":v["counts_by_original_lifespan"],
                  "estimable":v["estimable_two_strata_min_three_source_taxa"],
                  "equal":v.get("species_equal_original_source"),
                  "occupancy":v.get("occupancy_weighted_original_source"),
                  "gap":v.get("weight_scheme_gap"),
                  "flip_p_gap":v.get("original_exact_section_label_flip_p_gap"),
                  "bonf":v.get("exploratory_Bonferroni_across_all_estimable_strategy_gaps")}
            for name,v in out["all_original_strategy_controls_including_missing_strata"].items()}
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
