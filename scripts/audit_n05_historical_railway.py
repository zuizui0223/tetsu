#!/usr/bin/env python3
"""Audit official MLIT N05-25 nationwide passenger railway HISTORICAL GIS.

The 2025 N05 contains opening and closing/change attributes back to 1950,
unlike the N02 present-day snapshots. This is a provenance and cohort
screen, NOT a seed-dispersal or plant flowering outcome.

IMPORTANT: N05 is NONCOMMERCIAL-ONLY due to third-party JTB source rights.
Never upload original N05 shapefiles/line coordinates or assume freight-only
lines are included. An N05 record ending can be an ATTRIBUTE change, not
necessarily a railway closure. See docs/n05_historical_rail_protocol.md.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from datetime import datetime,timezone
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path
import requests
import shapefile

SOURCE="https://nlftp.mlit.go.jp/ksj/gml/data/N05/N05-25/N05-25_GML.zip"
INFO="https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N05-2025.html"
PUBLICATION_YEAR=2025
MAX_ZIP_BYTES=110_000_000
TREATMENT_START=2005


def known_year(val):
    """Normalize N05 year tokens; 999=unknown, 9999=current."""
    if val is None:
        return None
    # N05 attribute often strings "2020" or "2020-01-01".
    s=str(val).strip()
    if s in ("", "999", "9999", "0", "None", "NULL"):
        return None
    m=re.fullmatch(r"(\d{4})(?:[-/]\d{2}[-/]\d{2})?",s)
    if not m:
        return None
    year=int(m.group(1))
    return year if 1950<=year<=PUBLICATION_YEAR else None


def is_current_end(value):
    return str(value).strip() in ("9999","9999-12-31")


def bounded_zip_download():
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-rail-history-ecology/1.1 (https://github.com/zuizui0223/tetsu)"})
    with client.get(SOURCE,timeout=(25,180),stream=True) as response:
        response.raise_for_status()
        content=bytearray()
        for chunk in response.iter_content(1024*1024):
            content.extend(chunk)
            if len(content)>MAX_ZIP_BYTES:
                raise ValueError("N05 ZIP exceeded safe maximum")
    result=bytes(content)
    if not zipfile.is_zipfile(io.BytesIO(result)):
        raise ValueError("Official N05 source did not return a ZIP")
    return result


def shape_reader(archive, token):
    """Read SHP/DBF from ZIP without saving MLIT noncommercial files."""
    items={p.lower():p for p in archive.namelist()}
    candidates=[p for p in archive.namelist()
                if p.lower().endswith(".shp") and
                token.lower() in p.lower()]
    if not candidates:
        raise ValueError(f"N05 expected {token} .shp not found; shapefiles={list(items)[:25]}")
    # Prefer UTF-8 variant when multiple geometry formats included.
    candidates.sort(key=lambda x:(0 if "utf-8" in x.lower() else 1,len(x)))
    shp=candidates[0]
    stem=shp[:-4]
    def buf(ext):
        key=(stem+ext).lower()
        if key not in items:
            raise ValueError(f"Missing shapefile sidecar {key}")
        return io.BytesIO(archive.read(items[key]))
    reader=shapefile.Reader(shp=buf(".shp"),shx=buf(".shx"),
                            dbf=buf(".dbf"),encoding="utf-8",
                            encodingErrors="replace")
    return reader, shp


def audit_records(reader, label):
    names=[x[0] for x in reader.fields[1:]]
    must={"N05_004","N05_005b","N05_005e","N05_006"}
    if not must.issubset(names):
        raise ValueError(f"{label} missing expected history fields {sorted(must-set(names))}")
    counters=Counter()
    openings=Counter()
    final_ends=Counter()
    by_group=defaultdict(list)
    unknown_open_samples=Counter()
    by_corridor_open_year=defaultdict(set)
    for raw in reader.iterRecords():
        row=dict(zip(names,raw))
        open_raw=row["N05_004"]
        b_raw=row["N05_005b"]
        e_raw=row["N05_005e"]
        opening=known_year(open_raw)
        begin=known_year(b_raw)
        ending=known_year(e_raw)
        current=is_current_end(e_raw)
        group=str(row.get("N05_006") or "").strip()
        if not group:
            counters["missing_relation_group_id"]+=1
        if opening is not None:
            openings[opening]+=1
            if group:
                by_corridor_open_year[group].add(opening)
        else:
            counters["unknown_or_pre1950_opening_year"]+=1
            unknown_open_samples[str(open_raw)[:20]]+=1
        if begin is not None:
            counters["known_start_year_row"]+=1
        if ending is not None:
            final_ends[ending]+=1
            counters["ended_rows_including_attribute_versions"]+=1
        if current:
            counters["current_end_sentinel_9999"]+=1
        else:
            counters["not_current_end_sentinel"]+=1
        if begin is not None and ending is not None and begin>ending:
            counters["inconsistent_start_after_end"]+=1
        if group:
            by_group[group].append((opening,begin,ending,current))
        counters["n_rows"]+=1

    # First-year counts are deduped per official relation ID, not per
    # SHP row: multiple versions can share one logical railway segment.
    first_by_group={}
    current_groups=set()
    potential_terminal=Counter()
    for group, versions in by_group.items():
        known={year for opening,b,e,current in versions
               for year in ([opening] if opening is not None else [])}
        if len(known)>1:
            counters["inconsistent_opening_year_within_relation_group"]+=1
        if known:
            first_by_group[group]=min(known)
        if any(v[3] for v in versions):
            current_groups.add(group)
        terminal=[v[2] for v in versions if v[2] is not None]
        if terminal and not any(v[3] for v in versions):
            potential_terminal[max(terminal)]+=1

    out={
        "shape_type":label,
        "n_feature_rows":counters["n_rows"],
        "n_distinct_relation_groups":len(by_group),
        "n_relation_groups_current_2025":len(current_groups),
        "n_relation_groups_with_known_opening":len(first_by_group),
        "first_opening_year_groups_by_decade":dict(sorted(Counter(
            (year//10)*10 for year in first_by_group.values()
        ).items())),
        "candidate_relation_groups_newly_opened_2005_2025":{
            str(y):sum(year==y for year in first_by_group.values())
            for y in range(2005,PUBLICATION_YEAR+1)
        },
        "noncurrent_group_terminal_end_years_2005_2024":{
            str(y):potential_terminal[y] for y in range(2005,2025)
        },
        "raw_row_opening_year_2005_2025":{
            str(y):openings[y] for y in range(2005,2026)
        },
        "raw_row_ending_year_2005_2024_INCLUDES_ATTRIBUTE_CHANGES":{
            str(y):final_ends[y] for y in range(2005,2025)
        },
        "source_diagnostics":dict(counters),
        "unknown_opening_tokens_max5":unknown_open_samples.most_common(5),
        "warning":(
            "Relation group IDs are GIS metadata, not independent railway "
            "lines. Temporal end can be an attribute change, not physical "
            "closure. Candidate terminal end needs route/service verification. "
            "2025 and 1950/opening bounds not necessarily historically complete."
        )
    }
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="results/audit/n05_2025_historical_rail_summary.json")
    args=parser.parse_args()
    zip_raw=bounded_zip_download()
    result={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "source":SOURCE,"official_data_page":INFO,
        "sha256":hashlib.sha256(zip_raw).hexdigest(),
        "downloaded_zip_bytes":len(zip_raw),
        "copyright_license":"NONCOMMERCIAL ONLY; original licensed railway geometry not republished",
        "min_valid_history_year":1950,"max_valid_history_year":2025,
        "freight_only_and_yards_are_excluded":True,
        "outputs":{},"failures":{}
    }
    with zipfile.ZipFile(io.BytesIO(zip_raw)) as archive:
        for shape in ("RailroadSection2","Station2"):
            try:
                reader, member=shape_reader(archive,shape)
                out=audit_records(reader,shape)
                out["source_shapefile_basename"]=Path(member).name
                reader.close()
                result["outputs"][shape]=out
            except Exception as exc:
                result["failures"][shape]=f"{type(exc).__name__}: {str(exc)[:320]}"
    output=Path(args.output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    for shape,values in result["outputs"].items():
        print(json.dumps({
            "type":shape,
            "rows":values["n_feature_rows"],
            "groups":values["n_distinct_relation_groups"],
            "new_opened_group_2010_2025":sum(
                n for yr,n in values["candidate_relation_groups_newly_opened_2005_2025"].items()
                if int(yr)>=2010
            ),
            "possible_terminal_groups_2010_2024":sum(
                n for yr,n in values["noncurrent_group_terminal_end_years_2005_2024"].items()
                if int(yr)>=2010
            ),
            "data_errors":values["source_diagnostics"].get(
                "inconsistent_start_after_end",0
            ),
        },ensure_ascii=False),flush=True)
    if not result["outputs"]:
        raise SystemExit("N05 historical GIS source audit failed entirely")
    if "RailroadSection2" not in result["outputs"]:
        raise SystemExit("No verified N05 historical railway line geometries")


if __name__=="__main__":
    main()
