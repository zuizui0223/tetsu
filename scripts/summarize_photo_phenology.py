#!/usr/bin/env python3
"""Audit photo-based phenophases from structured survey visits.

For opportunistic photos, missing nonflowering plants are not observable and
this script MUST NOT be used to claim true flowering onset/offset.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.phenology import load_events, normalize_visits, weekly_summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--visits", required=True)
    p.add_argument("--events", default=None)
    p.add_argument("--output-dir", default="results/phenology")
    args = p.parse_args()
    cut_events = load_events(args.events) if args.events else {}
    rows, skipped = normalize_visits(args.visits, cut_events)
    if not rows:
        raise SystemExit("no eligible complete visits with reviewed phenophases")
    weekly = weekly_summary(rows)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for fname, data in (
        ("eligible_visits.csv", rows),
        ("weekly_observed_flowering.csv", weekly)
    ):
        with (output / fname).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    manifest = {
        "type": "descriptive phenophase audit; NOT flowering onset, causal mowing effect or train adaptation",
        "eligible_complete_visits": len(rows),
        "distinct_sites": len({x["site_id"] for x in rows}),
        "distinct_taxa": len({x["taxon"] for x in rows}),
        "excluded_rows": skipped,
        "n_weekly_strata": len(weekly),
        "n_screening_only_weekly_strata": sum(x["screening_only"] for x in weekly),
        "warning": "Only repeated, complete survey visits include a valid observed nonflowering denominator."
    }
    (output / "audit.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
