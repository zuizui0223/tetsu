#!/usr/bin/env python3
"""Audit independently sourced public data BEFORE attempting biological inference.

Does not claim seed dispersal, assemblage similarity or railway effects.
Outputs JSON with archive checksum, sizes, feature counts and API counts.
Raw downloaded third-party files are not committed to the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import requests
from openpyxl import load_workbook

ZENODO = "https://zenodo.org/records/7597598/files/AlienPlantsJapan_v1_FirstRecordPathway.xlsx?download=1"
ZENODO_MD5 = "312eb8548df4a94946d026450d975f04"
MLIT = "https://nlftp.mlit.go.jp/ksj/gml/data/N02/N02-25/N02-25_GML.zip"
GBIF = "https://api.gbif.org/v1/occurrence/search"

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "tetsu-reproducible-research-pilot/0.2 (https://github.com/zuizui0223/tetsu)"})
TIMEOUT = (15, 100)


def fetch_bytes(url: str) -> bytes:
    response = SESSION.get(url, timeout=TIMEOUT)
    response.raise_for_status()
    return response.content


def checksum(raw: bytes, algo: str = "sha256") -> str:
    return hashlib.new(algo, raw).hexdigest()


def audit_alien_workbook(raw: bytes) -> dict:
    md5 = checksum(raw, "md5")
    if md5 != ZENODO_MD5:
        raise ValueError(f"Zenodo workbook md5 mismatch: got {md5}; expected {ZENODO_MD5}")
    wb = load_workbook(io.BytesIO(raw), read_only=True, data_only=True)
    out = {
        "source": ZENODO,
        "doi": "10.5281/zenodo.7597598",
        "bytes": len(raw), "md5": md5, "sha256": checksum(raw),
        "sheets": []
    }
    for ws in wb.worksheets:
        head = list(ws.iter_rows(min_row=1, max_row=min(12, ws.max_row), max_col=min(30, ws.max_column), values_only=True))
        candidate_headers = [
            {"row": i + 1, "values": [str(v)[:72] if v is not None else "" for v in row]}
            for i, row in enumerate(head)
            if any(v is not None for v in row)
        ]
        out["sheets"].append({
            "name": ws.title, "max_row": ws.max_row, "max_col": ws.max_column,
            "header_candidates": candidate_headers[:4]
        })
    wb.close()
    return out


def audit_rail_archive(raw: bytes) -> dict:
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("N02 source did not return a valid zip archive")
    with zipfile.ZipFile(io.BytesIO(raw)) as bundle:
        members = bundle.namelist()
        shp_members = [p for p in members if p.lower().endswith(".shp")]
        section_shp = next((p for p in shp_members if "railroadsection" in p.lower()), None)
        station_shp = next((p for p in shp_members if "station" in p.lower()), None)
        result = {
            "source": MLIT, "reference_year": 2025, "bytes": len(raw),
            "sha256": checksum(raw), "archive_members": len(members),
            "shapefiles": shp_members,
            "n_sections": None, "n_station_shapes": None,
            "route_categories": None,
            "limitations": ["Passenger railway and tram routes only; freight-only lines and yards are not included.",
                            "Station shapes are rail line features; do not assume unique station points."]
        }
        if section_shp:
            import shapefile
            base = section_shp[:-4]
            names = {p.lower(): p for p in members}
            def get(ext):
                return io.BytesIO(bundle.read(names[(base + ext).lower()]))
            reader = shapefile.Reader(shp=get(".shp"), shx=get(".shx"), dbf=get(".dbf"), encoding="cp932")
            result["n_sections"] = len(reader)
            field_names = [field[0] for field in reader.fields[1:]]
            result["fields"] = field_names
            if "N02_001" in field_names:
                ix = field_names.index("N02_001")
                counts = Counter(str(record[ix]) for record in reader.iterRecords())
                result["route_categories"] = dict(counts.most_common())
            reader.close()
        if station_shp:
            import shapefile
            base = station_shp[:-4]
            names = {p.lower(): p for p in members}
            def get(ext):
                return io.BytesIO(bundle.read(names[(base + ext).lower()]))
            reader = shapefile.Reader(shp=get(".shp"), shx=get(".shx"), dbf=get(".dbf"), encoding="cp932")
            result["n_station_shapes"] = len(reader)
            reader.close()
        return result


def audit_gbif() -> dict:
    # Aggregate API count only, no records/download DOI. A count is NOT a
    # usable or citation-complete species incidence table.
    filters = [
        {"country": "JP", "taxonKey": 6, "hasCoordinate": "true", "limit": 0},
        {"country": "JP", "taxonKey": 6, "hasCoordinate": "true",
         "hasGeospatialIssue": "false", "year": "2018,2025", "limit": 0},
    ]
    queries = []
    for params in filters:
        response = SESSION.get(GBIF, params=params, timeout=TIMEOUT)
        response.raise_for_status()
        payload = response.json()
        queries.append({"parameters": params, "count": payload.get("count"), "api_url": response.url})
    return {
        "source": GBIF, "queries": queries,
        "limitations": ["taxonKey 6 = Plantae, not restricted to vascular plants.",
                        "The API counts do not indicate recorder balance or species coverage.",
                        "A DOI-bearing GBIF bulk download is required for a full inferential study."]
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/audit/public_sources.json")
    parser.add_argument("--strict", action="store_true",
                        help="Return nonzero when a required source fails")
    args = parser.parse_args()

    report = {
        "audit_time_utc": datetime.now(timezone.utc).isoformat(),
        "type": "source feasibility audit, NOT a train dispersal test",
        "datasets": {}, "errors": {}
    }
    for name, fn in (
        ("japan_alien_plant_checklist", lambda: audit_alien_workbook(fetch_bytes(ZENODO))),
        ("japan_rail_2025", lambda: audit_rail_archive(fetch_bytes(MLIT))),
        ("gbif_counts", audit_gbif)
    ):
        try:
            report["datasets"][name] = fn()
        except Exception as exc:
            report["errors"][name] = f"{type(exc).__name__}: {str(exc)[:350]}"
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if args.strict and any(x in report["errors"] for x in ("japan_alien_plant_checklist", "japan_rail_2025")):
        raise SystemExit("Required public-source audit failed; see JSON for details")


if __name__ == "__main__":
    main()
