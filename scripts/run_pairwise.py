#!/usr/bin/env python3
"""Descriptive overlaps only. There is NO causal or inferential p-value output."""
import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.pairwise import prepare_pairwise, read_incidence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--incidence", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    incidence = read_incidence(args.incidence)
    records = prepare_pairwise(incidence, args.pairs)
    if not records:
        raise ValueError("no pairs")
    with (out / "pairwise.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    valid = [r for r in records if r["alien_minus_native"] is not None]
    summary = {
        "disclaimer": "Exploratory descriptive overlaps, not evidence of train transport.",
        "n_regions": len(incidence), "n_pairs": len(records),
        "n_valid_contrasts": len(valid),
        "n_rail_connected": sum(r["rail_connected"] for r in records)
    }
    for label, connected in (("rail_connected", 1), ("rail_unconnected", 0)):
        values = [r["alien_minus_native"] for r in valid if r["rail_connected"] == connected]
        summary[label + "_n"] = len(values)
        summary[label + "_mean_alien_minus_native_jaccard"] = (
            statistics.mean(values) if values else None
        )
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
