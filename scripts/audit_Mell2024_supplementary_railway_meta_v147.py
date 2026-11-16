#!/usr/bin/env python3
"""v14.7: audit ORIGINAL Mell et al. (2024) systematic-review supplement.

DOI: 10.1186/s13750-024-00328-3.
This source has 316 included studies, 205 meta-analysis cases, and reports
1% of its meta-analysis cases from railways. The study examines verge
HABITAT contrasts, not proof of vehicle-borne seed dispersal.

Read publisher original XLSX files 5 (study metadata) and 8 (meta cases)
with the Python standard library OOXML XML parser, preserving row and
sheet labels. No original PDF or publisher workbook is redistributed.
Diagnostic rail-token row matches are CANDIDATES, not yet validated rail
effect estimates, and overlap across two supplements must not be counted
as independent original studies. Store only short source row descriptors.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import posixpath
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
DOC_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
RAIL = re.compile(r"\brail(?:ways?|roads?|tracks?)?\b|\btrains?\b", re.I)
DOI = "10.1186/s13750-024-00328-3"
SOURCE = "https://link.springer.com/article/" + DOI
FILES = {
    "study_metadata_5": "13750_2024_328_MOESM5_ESM.xlsx",
    "effect_cases_8": "13750_2024_328_MOESM8_ESM.xlsx",
}


def xlsx_cell_index(address):
    match = re.match(r"^([A-Z]+)[1-9]\d*$", address.upper())
    if not match:
        raise ValueError("Unsafe or invalid original Excel A1 cell reference: " + address[:50])
    n = 0
    for c in match.group(1):
        n = 26 * n + ord(c) - ord("A") + 1
    if n > 512:
        raise ValueError("Publisher spreadsheet >512-column safety limit")
    return n - 1


def get_shared_strings(z):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in node.findall(".//{" + MAIN + "}t"))
            for node in root.findall("{" + MAIN + "}si")]


def sheet_files(z):
    if "xl/workbook.xml" not in z.namelist():
        raise ValueError("Not a valid publisher OOXML XLSX")
    book = ET.fromstring(z.read("xl/workbook.xml"))
    rels_path = "xl/_rels/workbook.xml.rels"
    rels = {}
    if rels_path in z.namelist():
        xr = ET.fromstring(z.read(rels_path))
        for r in xr.findall("{" + DOC_REL + "}Relationship"):
            rels[r.get("Id")] = r.get("Target")
    sheets = []
    for i, sh in enumerate(book.findall(".//{" + MAIN + "}sheet")):
        name = sh.get("name") or f"sheet_{i + 1}"
        target = rels.get(sh.get("{" + REL + "}id"))
        if target:
            path = (target.lstrip("/") if target.startswith("/")
                    else posixpath.normpath(posixpath.join("xl", target)))
        else:
            path = f"xl/worksheets/sheet{i+1}.xml"
        if path not in z.namelist() or not path.startswith("xl/worksheets/"):
            raise ValueError("Missing or unsafe workbook sheet path")
        sheets.append((name, path))
    if not sheets:
        raise ValueError("Publisher workbook contains no worksheets")
    return sheets


def worksheet_rows(z, sheet_path, shared):
    root = ET.fromstring(z.read(sheet_path))
    for row in root.findall(".//{" + MAIN + "}sheetData/{" + MAIN + "}row"):
        vals = {}
        for cell in row.findall("{" + MAIN + "}c"):
            index = xlsx_cell_index(cell.get("r", ""))
            kind = cell.get("t")
            if kind == "inlineStr":
                inline = cell.find("{" + MAIN + "}is")
                val = "" if inline is None else "".join(
                    x.text or "" for x in inline.findall(".//{" + MAIN + "}t"))
            else:
                v = cell.find("{" + MAIN + "}v")
                val = "" if v is None or v.text is None else v.text
                if kind == "s" and val:
                    n = int(val)
                    if n < 0 or n >= len(shared):
                        raise ValueError("Invalid original publisher shared-string index")
                    val = shared[n]
            val = str(val).strip()
            if val:
                vals[index] = val
        yield (int(row.get("r")), vals)


def short_cells(vals, max_cols=42, max_chars=108):
    return [[n + 1, value[:max_chars]] for n, value in sorted(vals.items())[:max_cols]]


def examine(path, kind):
    raw = Path(path).read_bytes()
    if len(raw) < 2500 or not zipfile.is_zipfile(path):
        raise ValueError("Not an original publisher downloadable Excel workbook")
    sha = hashlib.sha256(raw).hexdigest()
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:
            raise ValueError("Original Springer appendix XLSX ZIP CRC failure")
        shared = get_shared_strings(z)
        sheets = []
        for name, xml in sheet_files(z):
            matches, previews, n_rows = [], [], 0
            for row_id, vals in worksheet_rows(z, xml, shared):
                n_rows += 1
                if len(previews) < 8 and vals:
                    previews.append({"source_excel_row":row_id,
                                     "nonempty_cells":short_cells(vals, 30, 90)})
                matching = [[n + 1, value[:115]]
                            for n, value in vals.items() if RAIL.search(value)]
                if matching and len(matches) < 75:
                    matches.append({
                        "source_excel_row": row_id,
                        "rail_keyword_columns": matching[:12],
                        "other_nonempty_cells": short_cells(vals, 38, 110)
                    })
            sheets.append({"publisher_sheet": name,
                           "total_original_nonempty_or_empty_rows": n_rows,
                           "raw_rail_keyword_row_candidates_shown":len(matches),
                           "first_8_nonempty_source_rows":previews,
                           "candidate_rail_rows_for_manual_column_validation":matches})
    return {
        "supplement_kind":kind,
        "original_publisher_supplement_file":FILES[kind],
        "source_SHA256":sha,
        "source_download_bytes":len(raw),
        "source_worksheet_count":len(sheets),
        "source_sheets":sheets
    }


def run(paths):
    result = {
        "version":"v14.7_source_diagnostic_only",
        "source_article_DOI":DOI,
        "source_article":SOURCE,
        "source_316_included_studies_authors_report":316,
        "source_47_studies_205_meta_cases_authors_report":205,
        "source_rail_fraction_approx_one_percent_authors_report":True,
        "source_original_meta_comparator":"transportation verge versus away from infrastructure",
        "NO_railway_dispersal_effect_from_these_meta_cases":True,
        "NO_unvalidated_rail_token_row_count_as_study_count":True,
        "source_supplements":[examine(path,kind) for kind,path in paths.items()]
    }
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--study-metadata-5",required=True)
    p.add_argument("--effect-cases-8",required=True)
    p.add_argument("--output",required=True)
    args = p.parse_args()
    result = run({"study_metadata_5":args.study_metadata_5,
                  "effect_cases_8":args.effect_cases_8})
    dest = Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    brief = [{
      "kind": x["supplement_kind"], "sha256":x["source_SHA256"],
      "sheets":[{"name":s["publisher_sheet"],
                 "rows":s["total_original_nonempty_or_empty_rows"],
                 "candidate_keyword_rows":s["raw_rail_keyword_row_candidates_shown"],
                 "previews":s["first_8_nonempty_source_rows"][:5],
                 "candidate_rows":s["candidate_rail_rows_for_manual_column_validation"][:12]}
                for s in x["source_sheets"]]
    } for x in result["source_supplements"]]
    print(json.dumps({"complete":True,
                      "study_supplement_diagnostic_only":brief,
                      "no_rail_effect_claim_yet":True},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
