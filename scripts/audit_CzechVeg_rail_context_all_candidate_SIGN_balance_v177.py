#!/usr/bin/env python3
"""v17.7 source-backed sign/robustness audit of ALL 39-station CzechVeg taxa.

Input MUST be the completed original v17.4 GitHub Actions full-source
CzechVeg+Preslia report, not the 20 hand-picked positive taxa in Git.
This script does NOT reconstruct plot records or independently validate the
source; it uses the already independently verified 117,739 plot-header /
2,510,340 source-species-row report. No species are removed on q sign before
the signed denominator is formed. It separates source original-study
jackknife from replication and does not infer rail selection or p values.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

ORIGINAL = "BROAD_original_same_STUDY_region_DECADE_area_class_all_years"
RECENT = "BROAD_original_same_STUDY_region_DECADE_area_class_2010_onward"
STRICT = "STRICT_HABITAT_original_same_STUDY_region_DECADE_area_class_all_years"
Q = "CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence"
NR = "CzechVeg_same_study_region_decade_area_source_N_RAIL_text_plot_with_this_taxon"
NC = "CzechVeg_same_study_region_decade_area_source_N_ROAD_text_plot_with_this_taxon"
NS = "CzechVeg_same_study_source_N_distinct_publications_with_taxon"
STATION = "station_original_2020_22_n_of_39"
LLO = "CzechVeg_original_AUTHOR_STUDY_scored_loo_min_delta"
LHI = "CzechVeg_original_AUTHOR_STUDY_scored_loo_max_delta"
LOO_N = "CzechVeg_all_original_source_AUTHOR_STUDY_scored_loo_full_source_studies"
VERSION = "v17.7_original_CzechVeg_39station_SIGN_unfiltered_candidates_and_sensitivity"

def sign(x):
    if x is None:
        return "unestimated"
    if x > 1e-12:
        return "positive"
    if x < -1e-12:
        return "negative"
    return "zero"

def eligible_without_q_sign(r):
    """Frozen v17.4 opportunity/station gates, with q sign NOT used."""
    return (
        r.get(Q) is not None
        and r[STATION] >= 5
        and r[NR] + r[NC] >= 5
        and r[NS] >= 2
    )

def audit_panel(panel, modern, strict, ranked_names):
    eligible = []
    for species, row in sorted(panel.items()):
        if not eligible_without_q_sign(row):
            continue
        q = row[Q]
        m = modern.get(species, {}).get(Q)
        s = strict.get(species, {}).get(Q)
        ll, hh = row.get(LLO), row.get(LHI)
        loo = ("always_positive_after_one_study_delete"
               if ll is not None and sign(ll) == "positive" and sign(q) == "positive"
               else "always_negative_after_one_study_delete"
               if hh is not None and sign(hh) == "negative" and sign(q) == "negative"
               else "one_study_deletion_can_remove_or_reverse_effect")
        eligible.append({
            "species": species, "q": q, "q_sign": sign(q),
            "actual_Preslia_station_occupancy_n39": row[STATION],
            "original_CzechVeg_source_studies_with_species": row[NS],
            "original_CzechVeg_rail_context_species_occurrences": row[NR],
            "original_CzechVeg_road_context_species_occurrences": row[NC],
            "n_all_studies_in_leave_one_original_paper_out": row.get(LOO_N),
            "loo_min_q": ll, "loo_max_q": hh,
            "loo_description_NOT_confidence_interval": loo,
            "2010plus_q": m, "2010plus_sign": sign(m),
            "strict_HABITAT_q": s, "strict_HABITAT_sign": sign(s),
            "was_in_original_q_positive_top20": species in ranked_names,
        })
    counts = Counter(x["q_sign"] for x in eligible)
    if not eligible:
        raise ValueError("No source-sign-independent species meet original support gates")
    # Modern and strict panels are SMALL observational strata, not independent
    # validation of railway effects; report the denominator explicitly.
    def sensitivity(key):
        supported = [r for r in eligible if r[key] is not None]
        directional = [r for r in supported if sign(r[key]) in ("positive", "negative")]
        return {
            "n_source_signindependent_species_with_estimable_q": len(supported),
            "n_signed_nonzero_both_original_and_this_panel": len([
                r for r in directional if r["q_sign"] in ("positive", "negative")]),
            "n_same_nonzero_direction_as_original_full_panel": sum(
                r["q_sign"] == sign(r[key])
                for r in directional
                if r["q_sign"] in ("positive", "negative")),
            "n_direction_reversed": sum(
                r["q_sign"] != sign(r[key])
                for r in directional
                if r["q_sign"] in ("positive", "negative")),
            "observational_not_independent_replication": True,
        }
    with_no_sign = sorted(eligible, key=lambda r: (-abs(r["q"]), r["species"]))
    prior = [r for r in eligible if r["was_in_original_q_positive_top20"]]
    return eligible, {
        "n_original_rail_station_taxa_passing_station_plot_study_SUPPORT_NO_SIGN_SELECTION":len(eligible),
        "n_q_positive_UNSELECTED_directional_denominator":counts["positive"],
        "n_q_negative_UNSELECTED_directional_denominator":counts["negative"],
        "n_q_exactly_zero":counts["zero"],
        "n_original_q_positive_top20_recognized":len(prior),
        "n_eligible_source_taxa_positive_loo_for_every_single_original_study_deletion":sum(
            r["loo_description_NOT_confidence_interval"] == "always_positive_after_one_study_delete"
            for r in eligible),
        "n_eligible_source_taxa_negative_loo_for_every_single_original_study_deletion":sum(
            r["loo_description_NOT_confidence_interval"] == "always_negative_after_one_study_delete"
            for r in eligible),
        "n_original_top20_positive_loo_under_single_original_study_deletion":sum(
            r["loo_description_NOT_confidence_interval"] == "always_positive_after_one_study_delete"
            for r in prior),
        "modern_2010plus_sparse_same_original_study_sensitivity":sensitivity("2010plus_q"),
        "STRICT_label_only_very_small_panel_sensitivity":sensitivity("strict_HABITAT_q"),
        "original_top20_species_with_negative_or_zero_2010plus_score":sorted(
            r["species"] for r in prior if r["2010plus_q"] is not None and
            r["2010plus_sign"] != "positive"),
        "original_top20_species_with_negative_or_zero_STRICT_HABITAT_score":sorted(
            r["species"] for r in prior if r["strict_HABITAT_q"] is not None and
            r["strict_HABITAT_sign"] != "positive"),
        "top10_BY_ABSOLUTE_signed_q_UNSELECTED_FOR_SIGN":with_no_sign[:10],
        "no_species_level_p_values_and_no_phylogenetic_or_study_cluster_CI":True,
        "station_origin_and_context_are_NOT_replicated_random_sample_of_rail_ecosystems":True,
        "single_study_deletion_IS_NOT_independent_study_replication":True,
    }

def audit(report):
    if (report.get("complete") is not True or
        report.get("version") !=
            "v17.4_CzechVeg_original_plot_baselined_39Preslia_station_taxa_railroad_source_proxy" or
        report.get("n_station_source_taxa") != 761 or
        report.get("n_unique_publisher_stations") != 39):
        raise ValueError("Not the actual completed v17.4 CzechVeg+Preslia39 original source report")
    orig = report["source_original_observational_rail_vs_road_comparisons"][ORIGINAL]
    if (orig["n_authentic_same_study_plot_balance_blocks"] != 44 or
        orig["n_original_rail_text_plots"] != 229 or
        orig["n_original_road_text_plots"] != 462 or
        orig["n_original_source_field_studies"] != 36):
        raise ValueError("True original broad CzechVeg 44 matched blocks / 36 original studies / 229:462 plots not reproduced")
    data = report["original_CzechVeg_full_country_plot_and_selected_taxa_source_counts"]
    if (data["source_original_CzechVeg_site_x_taxon_CSV_rows_scanned"] != 2510340 or
        data["source_original_plot_headers_with_recorded_plot_area_1_to_1000m2"] != 106670):
        raise ValueError("Source Zenodo 2.51m / 106670 original area-valid plot count not verified")
    frames = report["source_original_station_taxa_independent_CzechVeg_reference_context_per_species"]
    base = frames[ORIGINAL]
    if len(base) != 761:
        raise ValueError("Original 761 exact normalized Preslia species panel incomplete")
    ranked_names = {
        r["scientific_name_original_exact_normalized"]
        for r in report["independent_SOURCE_TAXON_RAIL_TEXT_CONTEXT_ranked_TOP20_NOT_true_railroad_affinity"]
    }
    rows, stat = audit_panel(base, frames[RECENT], frames[STRICT], ranked_names)
    expected = report["source_eligibility_criteria_AND_no_forced_rank20"][
        "n_eligible_station_taxa_with_positive_independent_CzechVeg_railroad_proxy_and_source_support"]
    if stat["n_q_positive_UNSELECTED_directional_denominator"] != expected or expected != 104:
        raise ValueError("All-source sign-independently eligible taxa do not reconstruct 104 original q-positive candidates")
    if stat["n_original_q_positive_top20_recognized"] != 20:
        raise ValueError("Prior original top20 names were not recovered")
    echium = next((r for r in rows if r["species"] == "Echium vulgare"), None)
    if echium is None or abs(echium["q"]-0.06909200397868588) > 1e-9:
        raise ValueError("Prespecified Echium original same-study q missing/drifted")
    info = {
        "complete": True, "version": VERSION,
        "source_genuine_v174_prior_successful_run_ID": 43700792828,
        "source_input_source_MD5s_preserved_from_original_v174_report": {
            "CzechVeg_headers": report["CzechVeg_full_original_source_HEADER_MD5_verified"],
            "CzechVeg_species": report["CzechVeg_full_original_source_2510340_species_MD5_verified"]
        },
        "original_source_CzechVeg_and_Preslia_papers": report["independent_original_botanist_sources"],
        "source_44_matched_blocks_36_studies_229_rail_462_road": {
            "blocks": 44, "studies": 36, "rail_releves": 229,
            "road_releves": 462
        },
        "all_source_SIGN_UNFILTERED_support_eligible": stat,
        "prespecified_Echium_vulgare_source_diagnostic": {
            **echium,
            "unadjusted_overall_RAIL_minus_ROAD_plot_prevalence_NOT_same_block_adjusted": (
                echium["original_CzechVeg_rail_context_species_occurrences"]/229
                -echium["original_CzechVeg_road_context_species_occurrences"]/462
            ),
            "q_is_within_author_study_region_decade_area_matched_prevalence_difference":True,
        },
        "scope": [
            "ALL eligible species are counted before looking at the sign of q; the earlier 104 positive taxa were necessarily q>0 by selection, not evidence of a 104/761 global positive signal.",
            "The signs are correlated among taxa; the source study and plots are not exchangeable random rail-versus-road samples. No binomial sign-test p, no ecological causal estimates or cross-country generalization.",
            "LOO min/max is deletion sensitivity across original source papers, NOT 36 independent taxon effect replications, 95% CI, random effects, or proof of selection.",
            "Even source strict HABITAT and 2010+ comparison use only 20/25 and 32/9 original source plots respectively, and may be missing/biased.",
            "Original 39 source stations aggregate multiple microhabitats; comparators are smaller fixed-area original botanical relevés; 2023 OSM line proximity is separate from actual rail management.",
            "This is a second-stage arithmetic reanalysis of the genuine v17.4 report, not a fresh ingestion of the original full 2.51m records. Source verified there.",
            "No photo phenotype was read, and no plant plasticity, morphological selection, slipstream stress, dispersal or genetic adaptation was estimated."
        ],
    }
    return rows, info

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source-v174", required=True)
    p.add_argument("--output-json", required=True)
    p.add_argument("--output-csv", required=True)
    a = p.parse_args()
    report = json.loads(Path(a.source_v174).read_text(encoding="utf-8"))
    rows, info = audit(report)
    out = Path(a.output_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    dest = Path(a.output_csv)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "version": info["version"],
        "support_sign_unfiltered":info["all_source_SIGN_UNFILTERED_support_eligible"],
        "prespecified_Echium":info["prespecified_Echium_vulgare_source_diagnostic"]
    }, ensure_ascii=False))

if __name__ == "__main__":
    main()
