#!/usr/bin/env python3
"""Descriptive year-on-year reference sensitivity from audited Kashmir JSON.

No new field resampling or evidence of a train-mediated causal effect.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

SOURCE_BLOB_SHA1 = "75dd78dd2b66a67543c248ff686c3a9dc016ea73"
ZONE_PAIRS = (
    ("Track", "Interior Plot 2"),
    ("Track", "Interior Plot 1"),
    ("Embankment", "Interior Plot 2"),
)
GROUPS = ("station", "between")


def git_blob_sha1(raw: bytes) -> str:
    """Compute Git's blob SHA-1 for a source byte string."""
    header = b"blob " + str(len(raw)).encode("ascii") + b"\0"
    return hashlib.sha1(header + raw).hexdigest()


def extract(data: dict) -> dict:
    """Aggregate preexisting source-derived means, not paired site changes."""
    x = data["all_species"]
    index: dict[tuple[str, tuple[str, str], int], dict] = {}
    for record in x:
        key = (record["group"], tuple(record["zones"]), int(record["year"]))
        if key in index:
            raise ValueError(f"Duplicate source stratum: {key}")
        index[key] = record
    expected = {
        (g, zones, year)
        for g in GROUPS for zones in ZONE_PAIRS for year in (2014, 2017)
    }
    if set(index) != expected:
        raise ValueError(f"Unexpected source strata: missing={expected-set(index)} extra={set(index)-expected}")
    out = []
    for group in GROUPS:
        for zones in ZONE_PAIRS:
            a, b = index[(group, zones, 2014)], index[(group, zones, 2017)]
            if a["n"] < 2 or b["n"] < 2:
                raise ValueError("Insufficient source site count")
            out.append({
                "group": group,
                "rail_zone": zones[0],
                "reference_zone": zones[1],
                "site_count_2014": a["n"],
                "site_count_2017": b["n"],
                "raw_delta_2014": a["all_species_distance_delta"],
                "raw_delta_2017": b["all_species_distance_delta"],
                "raw_temporal_change": b["all_species_distance_delta"] - a["all_species_distance_delta"],
                "equalized_delta_2014": a["richness_equalized_delta"],
                "equalized_delta_2017": b["richness_equalized_delta"],
                "equalized_temporal_change": b["richness_equalized_delta"] - a["richness_equalized_delta"],
            })
    station_vs_between = []
    for rail, reference in ZONE_PAIRS:
        s = next(row for row in out if row["group"] == "station" and row["rail_zone"] == rail and row["reference_zone"] == reference)
        b = next(row for row in out if row["group"] == "between" and row["rail_zone"] == rail and row["reference_zone"] == reference)
        station_vs_between.append({
            "rail_zone": rail,
            "reference_zone": reference,
            "raw_station_minus_between_temporal_change": s["raw_temporal_change"] - b["raw_temporal_change"],
            "equalized_station_minus_between_temporal_change": s["equalized_temporal_change"] - b["equalized_temporal_change"],
        })
    return {
        "method": "descriptive difference of preaggregated source-site Jaccard distance contrasts",
        "input_is_derived_not_original_field_matrix": True,
        "different_year_site_membership_not_verified_from_aggregate": True,
        "paired_site_level_temporal_change_identified": False,
        "traffic_or_genetic_causal_effect_identified": False,
        "strata": out,
        "reference_sensitivity": station_vs_between,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1] / "research/validation_20261117/kashmir_results.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw = args.source.read_bytes()
    found = git_blob_sha1(raw)
    if found != SOURCE_BLOB_SHA1:
        raise SystemExit(f"Source Git blob SHA mismatch: expected {SOURCE_BLOB_SHA1}; actual {found}")
    summary = extract(json.loads(raw))
    summary["source_git_blob_sha1"] = found
    final = json.dumps(summary, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.write_text(final, encoding="utf-8")
    else:
        print(final, end="")


if __name__ == "__main__":
    main()
