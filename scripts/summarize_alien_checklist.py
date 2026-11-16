#!/usr/bin/env python3
"""Summarize Egawa & Koyama's real Japanese alien plant checklist.

Counts pathway tags as non-mutually exclusive. National introduction years
must not be interpreted as the date of arrival at a railway station.
"""
from __future__ import annotations

import argparse
import io
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from openpyxl import load_workbook

from audit_public_data import ZENODO, ZENODO_MD5, checksum, fetch_bytes


PATHWAY_COLS = (
    "Stowaway", "Contaminant", "Release", "EscapeAgriculture",
    "EscapeMedicinal", "EscapeOrnamental", "EscapeResearch",
    "EscapePetFood", "EscapeBotanicalgarden"
)
INTENTIONAL_COLS = ("Release", "EscapeAgriculture", "EscapeMedicinal",
                    "EscapeOrnamental", "EscapeResearch", "EscapePetFood",
                    "EscapeBotanicalgarden")
ACCIDENTAL_COLS = ("Stowaway", "Contaminant")


def flag(value) -> bool:
    """Count only an explicit 1 as presence, never classify NA as false."""
    return value == 1 or str(value).strip() == "1"


def valid_year(value) -> bool:
    try:
        value = int(value)
        return 1600 <= value <= 2026
    except (TypeError, ValueError):
        return False


def parse_checklist(raw: bytes) -> dict:
    md5 = checksum(raw, "md5")
    if md5 != ZENODO_MD5:
        raise ValueError("Zenodo workbook MD5 mismatch; refuse data")
    wb = load_workbook(io.BytesIO(raw), read_only=True, data_only=True)
    ws = wb["v1_FirstRecordPathway"]
    rows = ws.iter_rows(values_only=True)
    keys = [str(cell).strip() if cell is not None else "" for cell in next(rows)]
    missing = set(PATHWAY_COLS + ("Family", "ScientificName_w.oAuthor",
                                  "FirstRecordYear", "IntroductionYear", "PathwayUnknown")) - set(keys)
    if missing:
        raise ValueError(f"Missing expected fields: {sorted(missing)}")
    counts = Counter()
    by_pathway = Counter()
    seen = set()
    missing_name_rows = 0
    skipped_samples = []
    examples = {"accidental_only": [], "ornamental_only": []}
    for record in rows:
        row = dict(zip(keys, record))
        name = str(row.get("ScientificName_w.oAuthor") or "").strip()
        # The dataset has non-taxon explanatory rows; don't silently count them.
        family = str(row.get("Family") or "").strip()
        if not family or not name:
            missing_name_rows += 1
            if len(skipped_samples) < 25:
                skipped_samples.append({"family": family, "canonical": name})
            continue
        counts["n_valid_taxon_rows"] += 1
        if name in seen:
            counts["n_duplicate_canonical_names"] += 1
        seen.add(name)
        if valid_year(row["FirstRecordYear"]):
            counts["n_first_record_years"] += 1
        if valid_year(row["IntroductionYear"]):
            counts["n_introduction_years"] += 1
        if str(row["PathwayUnknown"] or "").strip().lower() == "y":
            counts["n_pathway_unknown_y"] += 1
        for key in PATHWAY_COLS:
            if flag(row[key]):
                by_pathway[key] += 1
        accidental = any(flag(row[k]) for k in ACCIDENTAL_COLS)
        intentional = any(flag(row[k]) for k in INTENTIONAL_COLS)
        if accidental and intentional:
            counts["accidental_and_intentional"] += 1
        elif accidental:
            counts["accidental_only"] += 1
            if len(examples["accidental_only"]) < 5:
                examples["accidental_only"].append(name)
        elif intentional:
            counts["intentional_only"] += 1
        else:
            counts["neither_group"] += 1
        if flag(row["EscapeOrnamental"]) and not accidental:
            counts["ornamental_without_accidental"] += 1
            if len(examples["ornamental_only"]) < 5:
                examples["ornamental_only"].append(name)
    wb.close()
    return {
        "source": ZENODO, "md5": md5,
        "n_canonical_names": len(seen), "non_taxon_or_blank_rows_skipped": missing_name_rows, "skipped_examples": skipped_samples,
        "summary": dict(counts), "pathway_flag_positive": dict(by_pathway),
        "examples_not_an_effect_estimate": examples,
        "notes": [
            "Pathway flags can overlap and therefore their counts are not mutually exclusive.",
            "An introduction pathway to Japan is not necessarily an introduction pathway to any specific railway station.",
            "GBIF taxon matching and local presence/absence screening remain to be performed.",
            "Known first national record is not the local arrival year."
        ]
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/audit/alien_pathway_summary.json")
    args = parser.parse_args()
    summary = parse_checklist(fetch_bytes(ZENODO))
    summary["audit_time_utc"] = datetime.now(timezone.utc).isoformat()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
