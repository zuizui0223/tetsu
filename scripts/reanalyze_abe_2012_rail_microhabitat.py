#!/usr/bin/env python3
"""Audited descriptive reanalysis of Abe et al. (2012) JR East weeds.

The 2011 railway questionnaire reported occurrence *mentions* by microhabitat
and broad region, not individuals and not independent experimental units.
No treatment causal effects or p-values. Read:
doi:10.3719/weed.57.159, Table 4 (journal p.161).
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = (
    "embankment_slope", "formation_base", "trackbed",
    "off_track", "boundary", "other"
)
REGION_TOTALS = {
    "warm": {
        "total_reported": 592,
        "embankment_slope": 247, "formation_base": 229,
        "trackbed": 16, "off_track": 30, "boundary": 10,
        "other": 60,
    },
    "cold": {
        "total_reported": 742,
        "embankment_slope": 395, "formation_base": 162,
        "trackbed": 54, "off_track": 16, "boundary": 5,
        "other": 110,
    },
}
SOURCE = "https://doi.org/10.3719/weed.57.159"
SOURCE_TABLE = "Abe et al. 2012 Table 4, JR East 2011 questionnaire: 109 rail sections"


def read_and_validate(path: str | Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "region", "source_name_ja", "taxon_group", "life_history",
            "total_reported", *FIELDS
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing columns: {sorted(required - set(reader.fieldnames or []))}")
        records = list(reader)
    if len(records) != 22:
        raise ValueError("Table 4 must have 11 named/other categories per region")
    keys = set()
    audit = defaultdict(Counter)
    for item in records:
        region = item["region"]
        if region not in REGION_TOTALS:
            raise ValueError(f"Unknown source region {region}")
        tag = (region, item["taxon_group"])
        if tag in keys:
            raise ValueError(f"Duplicate source category {tag}")
        keys.add(tag)
        if item["life_history"] not in ("annual", "perennial", "mixed_unknown"):
            raise ValueError("Unrecognized life history from transcribed source")
        for key in (*FIELDS, "total_reported"):
            try:
                value = int(item[key])
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Noninteger count {tag}:{key}") from exc
            if value < 0:
                raise ValueError(f"Negative count {tag}:{key}")
            item[key] = value
            audit[region][key] += value
        if sum(item[f] for f in FIELDS) != item["total_reported"]:
            raise ValueError(f"Habitat counts do not sum to species total for {tag}")
    for region, expected in REGION_TOTALS.items():
        if audit[region] != expected:
            raise ValueError(
                f"Published region {region} Table 4 subtotal mismatch: "
                f"{dict(audit[region])} != {expected}"
            )
        if (region, "OTHER_TOP10_RESIDUAL") not in keys:
            raise ValueError(f"Published 'other species' residual missing from {region}")
    return records


def row_by(records: list[dict], region: str, taxon: str) -> dict:
    found = [
        x for x in records
        if x["region"] == region and x["taxon_group"] == taxon
    ]
    if len(found) != 1:
        raise ValueError(f"Expected exactly one category {region}/{taxon}")
    return found[0]


def habitat_enrichment(records, taxon: str, region: str | None = None,
                       habitat: str = "trackbed") -> dict:
    """Reported occurrence fraction / all-habitat occurrence fraction.

    Enrichment is an observational composition ratio only. It is NOT
    a per-site odds ratio or survival/adaptation effect.
    """
    if habitat not in FIELDS:
        raise ValueError("Unknown source microhabitat")
    if region is not None and region not in REGION_TOTALS:
        raise ValueError("Unknown region")
    selected = [
        x for x in records if region is None or x["region"] == region
    ]
    subset = [x for x in selected if x["taxon_group"] == taxon]
    if not subset:
        raise ValueError(f"No published taxon category {taxon}")
    taxon_habitat = sum(x[habitat] for x in subset)
    n_habitat = sum(x[habitat] for x in selected)
    n_taxon = sum(x["total_reported"] for x in subset)
    all_reports = sum(x["total_reported"] for x in selected)
    if n_habitat == 0 or all_reports == 0 or n_taxon == 0:
        raise ValueError("Undefined source composition")
    fraction_habitat = taxon_habitat / n_habitat
    baseline = n_taxon / all_reports
    enrichment = fraction_habitat / baseline
    return {
        "taxon_group": taxon,
        "region": region or "pooled_warm_cold",
        "habitat": habitat,
        "taxon_occurrence_reports_in_habitat": taxon_habitat,
        "all_taxa_occurrence_reports_in_habitat": n_habitat,
        "taxon_occurrence_reports_all_habitats": n_taxon,
        "all_taxa_occurrence_reports_all_habitats": all_reports,
        "taxon_share_of_habitat_reports": round(fraction_habitat, 6),
        "taxon_share_of_all_reports": round(baseline, 6),
        "descriptive_enrichment": round(enrichment, 5),
        "not_independent_replicates": True
    }


def summarize(records):
    taxa = sorted({x["taxon_group"] for x in records
                   if x["taxon_group"] != "OTHER_TOP10_RESIDUAL"})
    region_reports = []
    for region in ("warm", "cold", None):
        for taxon in taxa:
            # Region top10 differ: do NOT treat taxon missing from this
            # reported top10 list as a biological zero in that region.
            if not any(x["taxon_group"] == taxon and
                       (region is None or x["region"] == region)
                       for x in records):
                continue
            for zone in ("trackbed", "embankment_slope", "formation_base"):
                region_reports.append(
                    habitat_enrichment(records, taxon, region, zone)
                )
    focus_taxa = (
        "Equisetum arvense", "Solidago altissima", "Setaria spp.",
        "Miscanthus sinensis", "Imperata cylindrica"
    )
    focus = {}
    for taxon in focus_taxa:
        focus[taxon] = {}
        for region in ("warm", "cold", None):
            if any(x["taxon_group"] == taxon and
                   (region is None or x["region"] == region)
                   for x in records):
                focus[taxon][region or "pooled"] = {
                    zone: habitat_enrichment(records, taxon, region, zone)
                    for zone in ("trackbed", "embankment_slope")
                }
    horsetail = focus["Equisetum arvense"]
    condition_ratios = {}
    for region in ("warm", "cold", "pooled"):
        track = horsetail[region]["trackbed"]
        slope = horsetail[region]["embankment_slope"]
        denom = slope["taxon_share_of_habitat_reports"]
        condition_ratios[region] = round(
            track["taxon_share_of_habitat_reports"] / denom, 5
        ) if denom else None
    return {
        "type": "independently computed, source-audited descriptive reanalysis",
        "source": SOURCE, "source_table": SOURCE_TABLE,
        "n_source_rail_sections": 109,
        "n_source_published_species_and_other_rows": len(records),
        "sample_sizes_are": (
            "reported weed occurrence mentions in microhabitat categories; "
            "NOT unique plants, survey sites, independent treatments, "
            "or observation effort per habitat"
        ),
        "warm_region_total_reports": REGION_TOTALS["warm"]["total_reported"],
        "cold_region_total_reports": REGION_TOTALS["cold"]["total_reported"],
        "all_region_occurrence_reports": 1334,
        "trackbed_reports": 70,
        "embankment_slope_reports": 642,
        "formation_base_reports": 391,
        "focus": focus,
        "horsetail_trackbed_vs_slope_ratio_of_report_shares": condition_ratios,
        "audits": {
            "row_marginals_matched": True,
            "warm_region_table_subtotals_matched": True,
            "cold_region_table_subtotals_matched": True
        },
        "limitations": [
            "Raw Table 4 is an aggregated 2011 railway operator questionnaire.",
            "One railway section can contribute more than one species and/or microhabitat count; no independent statistical test is licensed by these data.",
            "The warm/cold region top-10 species lists differ; absence from top 10 does not mean species absent. OTHER is an unresolved residual, not biological diversity.",
            "Different habitats have different maintenance and reporting/survey opportunities; trackbed enrichment is a composition association, not causally induced stress tolerance.",
            "No cutting dates, herbicide product/dose, slipstream velocities, photo ROIs, flowering dates or genetics are in these source counts.",
            "This is a NEW reanalysis of 2012 published aggregates, NOT new field measurements or evidence of train-carried seeds."
        ]
    }, region_reports


def save_csv(path, rows):
    if not rows:
        raise ValueError("empty analysis")
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default="data/literature/abe_2012_railway_microhabitat_occurrences.csv"
    )
    parser.add_argument("--output-dir", default="results/literature_rail_microhabitat")
    args = parser.parse_args()
    records = read_and_validate(args.input)
    summary, categories = summarize(records)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )
    save_csv(out / "species_habitat_composition.csv", categories)
    h = summary["focus"]["Equisetum arvense"]["pooled"]
    print(json.dumps({
        "all_source_tables_reconcile": summary["audits"],
        "rail_sections": summary["n_source_rail_sections"],
        "all_habitat_occurrence_reports": summary["all_region_occurrence_reports"],
        "trackbed_occurrences": summary["trackbed_reports"],
        "horsetail_trackbed_24_of_70": h["trackbed"]["taxon_share_of_habitat_reports"],
        "horsetail_all_habitats_98_of_1334": h["trackbed"]["taxon_share_of_all_reports"],
        "horsetail_trackbed_enrichment": h["trackbed"]["descriptive_enrichment"],
        "horsetail_track_vs_slope_share_ratio":
            summary["horsetail_trackbed_vs_slope_ratio_of_report_shares"],
        "report": str(out / "summary.json")
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
