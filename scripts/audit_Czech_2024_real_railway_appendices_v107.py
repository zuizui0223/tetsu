#!/usr/bin/env python3
"""v10.7 source audit for actual Czech railway station botanical survey appendices.

Kutlvašr et al. 2024, Preslia 96:247–266, doi 10.23855/preslia.2024.247.
No claim of raw Czech–Slovenia plot-area or habitat-zone equivalence.
We do NOT republish the authors' original worksheets: only safe metadata,
sheet headers, small schema examples, original file SHA256, and data-quality
counts in the JSON report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from collections import Counter
from pathlib import Path

SOURCE_DOI="10.23855/preslia.2024.247"
S1_URL="https://www.preslia.cz/appendices/P24302KutlvasrSupplTableS1.xlsx"
S2_URL="https://preslia.cz/appendix/file?id=540"


def clean_cell(x):
    if x is None:
        return None
    if isinstance(x,(int,float,bool,str)):
        return str(x)[:130]
    return f"<{type(x).__name__}>"


def audit_sheet(ws):
    # This is for discovering source schema, NOT an ecological model.
    rows=list(ws.iter_rows(min_row=1,max_row=min(12,ws.max_row),values_only=True))
    headers=[[clean_cell(v) for v in row[:55]] for row in rows]
    column_preview=[]
    for c in range(1,min(41,ws.max_column+1)):
        cells=[r[c-1] for r in rows if c-1 < len(r)]
        x=[clean_cell(v) for v in cells if v is not None]
        column_preview.append({
            "one_based_column_index":c,
            "n_nonempty_first_12":len(x),
            "first_three_nonempty":x[:3]})
    return {
        "sheet_name":ws.title,
        "reported_max_rows_including_possible_formatting":ws.max_row,
        "reported_max_columns_including_possible_formatting":ws.max_column,
        "first_eight_rows_first_55_cells":headers[:8],
        "first_40_columns_preview":column_preview,
        "actual_trait_coded_species_versus_presence_matrix_is_not_assumed":True,
    }


def audit_workbook(path):
    from openpyxl import load_workbook
    data=Path(path).read_bytes()
    if len(data)<500 or not zipfile.is_zipfile(path):
        raise ValueError("Publisher source not a valid nonempty OOXML XLSX")
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:
            raise ValueError("Original Czech appendix XLSX ZIP CRC mismatch")
    wb=load_workbook(path,read_only=True,data_only=True)
    sheets=[audit_sheet(wb[n]) for n in wb.sheetnames]
    wb.close()
    return {"original_source_n_bytes":len(data),
            "original_source_sha256":hashlib.sha256(data).hexdigest(),
            "original_source_sheet_count":len(sheets),
            "original_source_sheet_schemas":sheets,
            "NOT_original_republication_or_new_rail_plant_survey":True}


def run(s1_path,s2_path):
    s1=audit_workbook(s1_path)
    s2=audit_workbook(s2_path)
    return {
        "complete":True,
        "version":"v10.7",
        "article_DOI":SOURCE_DOI,
        "record_type":"original_czech_station_survey_appendix_schema_audit_ONLY",
        "first_appendix":{"publisher_url":S1_URL,**s1},
        "second_appendix":{"publisher_url":S2_URL,**s2},
        "no_Czech_vs_Slovenia_biological_effect_fitted_yet":True,
        "authors": "Kutlvašr, Turková, Hejda, Vojík, Kadlecová, Berchová Bímová, Pyšek & Pergl (2024)",
        "publication_reported_design":"39 railway localities in Czech Republic during 2020–2022; tracks/embankments vs wider surroundings vs buildings; published 11139 species records and 764 taxa (article metadata, not independently verified workbook count).",
        "limitations":"Independent Czech source zones and sampling units must be validated from supplement before comparison. We cannot pool Czech station floristics with Slovenia BS-BT or MIREN road quadrats as identical design. Original source files downloaded ephemerally, not redistributed.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--s1",required=True)
    p.add_argument("--s2",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    result=run(args.s1,args.s2)
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
