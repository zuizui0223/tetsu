#!/usr/bin/env python3
"""Small deterministic GBIF taxon/occurrence coverage check for railway-study taxa.

Checks availability of recent georeferenced observations in two independently
defined Japanese alien-plant pathway groups. No dispersal-effect estimates.
Not a representative GBIF download: results serve ONLY as feasibility screening.
"""
from __future__ import annotations

import argparse
import io
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
from openpyxl import load_workbook
from audit_public_data import ZENODO, fetch_bytes
from summarize_alien_checklist import flag, ACCIDENTAL_COLS, INTENTIONAL_COLS

MATCH = "https://api.gbif.org/v1/species/match"
OCC = "https://api.gbif.org/v1/occurrence/search"


def candidate_groups(workbook_bytes: bytes) -> dict[str, list[str]]:
    wb = load_workbook(io.BytesIO(workbook_bytes), read_only=True, data_only=True)
    ws = wb["v1_FirstRecordPathway"]
    iterator = ws.iter_rows(values_only=True)
    headers = [str(x) for x in next(iterator)]
    groups = {"accidental_only": [], "ornamental_no_accidental": []}
    for cells in iterator:
        row = dict(zip(headers, cells))
        name = str(row.get("ScientificName_w.oAuthor") or "").strip()
        if not name or not row.get("Family"):
            continue
        accidental = any(flag(row.get(k)) for k in ACCIDENTAL_COLS)
        intentional = any(flag(row.get(k)) for k in INTENTIONAL_COLS)
        if accidental and not intentional:
            groups["accidental_only"].append(name)
        if flag(row.get("EscapeOrnamental")) and not accidental:
            groups["ornamental_no_accidental"].append(name)
    wb.close()
    return {k: sorted(set(v)) for k, v in groups.items()}


def spaced_sample(items, n):
    if len(items) <= n:
        return items
    return [items[int((i + 0.5) * len(items) / n)] for i in range(n)]


def coverage(groups: dict[str, list[str]], n: int) -> dict:
    sess = requests.Session()
    sess.headers.update({"User-Agent": "tetsu-feasibility-check/0.2 (https://github.com/zuizui0223/tetsu)"})
    result = {}
    for group, names in groups.items():
        rows = []
        for name in spaced_sample(names, n):
            item = {"canonical_name": name, "matching_status": "unmatched", "observations": None}
            try:
                match = sess.get(MATCH, params={"name": name, "strict": "true"}, timeout=(15, 45))
                match.raise_for_status()
                data = match.json()
                mtype = data.get("matchType")
                conf = int(data.get("confidence") or 0)
                rank = data.get("rank")
                item.update({"matchType": mtype, "confidence": conf, "rank": rank,
                             "accepted_canonical": data.get("canonicalName")})
                if mtype == "EXACT" and conf >= 90 and rank in (
                    "SPECIES", "SUBSPECIES", "VARIETY", "FORM", "HYBRID", "CULTIVAR"
                ):
                    key = data.get("usageKey")
                    item["usageKey"] = key
                    count = sess.get(OCC, params={
                        "country": "JP", "taxonKey": key, "hasCoordinate": "true",
                        "hasGeospatialIssue": "false", "year": "2018,2025", "limit": 0
                    }, timeout=(15, 45))
                    count.raise_for_status()
                    item["observations"] = int(count.json().get("count") or 0)
                    item["matching_status"] = "matched_exact"
                else:
                    item["matching_status"] = "excluded_low_confidence_or_rank"
            except (requests.RequestException, ValueError) as exc:
                item["matching_status"] = "error"
                item["error"] = str(exc)[:180]
            rows.append(item)
            time.sleep(0.10)
        vals = [x["observations"] for x in rows if x["observations"] is not None]
        result[group] = {
            "group_population_size": len(names),
            "n_sampled_names": len(rows),
            "n_exact_high_confidence": len(vals),
            "n_with_at_least_5_records": sum(v >= 5 for v in vals),
            "n_with_at_least_20_records": sum(v >= 20 for v in vals),
            "n_with_at_least_50_records": sum(v >= 50 for v in vals),
            "sample": rows
        }
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names-per-group", type=int, default=12)
    ap.add_argument("--output", default="results/audit/gbif_taxon_coverage.json")
    args = ap.parse_args()
    groups = candidate_groups(fetch_bytes(ZENODO))
    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "description": "Deterministic evenly spaced species-name sampling; NOT a floristic analysis or train transport test.",
        "gbif_filter": "JP, taxon match EXACT confidence>=90, georeferenced no geospatial issue, years 2018-2025",
        "result": coverage(groups, args.names_per_group),
        "warnings": [
            "Full species data and independent observation-effort denominators are needed for any biogeographic model.",
            "Zero count in GBIF does not imply biological absence.",
            "Japanese national introduction pathway cannot be interpreted as a rail-specific pathway.",
            "Any sampling and taxonomy-match failure is explicitly retained."
        ]
    }
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
