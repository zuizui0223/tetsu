#!/usr/bin/env python3
"""Extract a *candidate* JR Freight container-handling station list.

Source: https://www.jrfreight.co.jp/service/area (public JR Freight page).
The list includes ORS (off-rail stations): these MUST NOT be treated as
rail-served freight depots. Coordinates and actual train routes are NOT inferred.
This is a source/coverage audit, not a test of seed dispersal.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import requests
from bs4 import BeautifulSoup

URL = "https://www.jrfreight.co.jp/service/area"
PREFECTURES = ("北海道", "青森県", "岩手県", "宮城県", "秋田県", "山形県",
              "福島県", "茨城県", "栃木県", "群馬県", "埼玉県", "千葉県",
              "東京都", "神奈川県", "新潟県", "富山県", "石川県", "福井県",
              "山梨県", "長野県", "岐阜県", "静岡県", "愛知県", "三重県",
              "滋賀県", "京都府", "大阪府", "兵庫県", "奈良県", "和歌山県",
              "鳥取県", "島根県", "岡山県", "広島県", "山口県", "徳島県",
              "香川県", "愛媛県", "高知県", "福岡県", "佐賀県", "長崎県",
              "熊本県", "大分県", "宮崎県", "鹿児島県", "沖縄県")


def clean(cell):
    return " ".join(cell.get_text(" ", strip=True).replace("\u3000", " ").split())


def extract_station_rows(html: str) -> tuple[list[dict], dict]:
    soup = BeautifulSoup(html, "html.parser")
    accepted = []
    stats = Counter()
    for table in soup.find_all("table"):
        text = clean(table)[:800]
        if "駅名" not in text or "県名" not in text:
            continue
        stats["candidate_tables"] += 1
        last_pref = ""
        for row in table.find_all("tr"):
            cells = [clean(x) for x in row.find_all(["th", "td"], recursive=False)]
            # The official table uses a postal-code cell for station rows.
            postal = next((i for i, x in enumerate(cells) if "〒" in x), -1)
            if postal < 1:
                continue
            station = cells[postal - 1]
            if any(x in station for x in ("駅名", "支店名", "営業部")) or not station:
                continue
            # Some rows omit prefecture due to HTML rowspan.
            if postal >= 2:
                candidate_pref = cells[postal - 2]
                if candidate_pref in PREFECTURES:
                    last_pref = candidate_pref
            pref = last_pref
            address = cells[postal]
            if not pref or not re.search(r"〒\s*\d{3}-\d{4}", address):
                stats["unresolved_candidate_rows"] += 1
                continue
            is_ors = bool(re.search(r"\bORS\b|オフレール", station, re.IGNORECASE))
            accepted.append({
                "prefecture": pref,
                "station": station.replace(" ", ""),
                "address": address,
                "off_rail_service_ors": int(is_ors),
                "port_name_candidate": int("港" in station or "埠頭" in station or "ふ頭" in station),
                "source_url": URL
            })
    # Detect ambiguous repeated facilities rather than silently dropping them.
    keys = [(x["prefecture"], x["station"]) for x in accepted]
    stats["duplicate_prefecture_station_rows"] = len(keys) - len(set(keys))
    return accepted, dict(stats)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="results/audit")
    args = parser.parse_args()
    response = requests.get(URL, timeout=(15, 90), headers={
        "User-Agent": "tetsu-source-audit/0.2 (https://github.com/zuizui0223/tetsu)"
    })
    response.raise_for_status()
    records, stats = extract_station_rows(response.text)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if records:
        with (out / "freight_station_candidates.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    report = {
        "audit_time_utc": datetime.now(timezone.utc).isoformat(),
        "source": URL,
        "sha256_of_source_html": hashlib.sha256(response.content).hexdigest(),
        "rows_extracted": len(records),
        "off_rail_service_ors_rows": sum(x["off_rail_service_ors"] for x in records),
        "non_ors_rows": sum(not x["off_rail_service_ors"] for x in records),
        "port_name_candidates": sum(x["port_name_candidate"] for x in records),
        "prefectures_covered": len(set(x["prefecture"] for x in records)),
        "sample_station_names": [x["station"] for x in records[:10]],
        "parse_diagnostics": stats,
        "caveats": [
            "Container handling site is not equivalent to directly rail-served freight station.",
            "ORS explicitly denotes off-rail service. Station rows may not equal JR Freight's corporate station count.",
            "Station name contains 港/埠頭 does not establish physical port adjacency.",
            "Independent geocoding, active route topology, year reconciliation and license checks required before geospatial use."
        ]
    }
    (out / "freight_station_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not records:
        raise SystemExit("No station rows were extracted: HTML changed, refusing silent success")


if __name__ == "__main__":
    main()
