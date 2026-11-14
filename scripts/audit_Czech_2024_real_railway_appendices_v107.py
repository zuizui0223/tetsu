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



def audit_presence_matrix_header_s1(s1_path,s2_path):
    """Source-only diagnostics for the unusual merged locality headers.

    Row 5 is the locality code with merged cells, row 6 is the source
    rec/hist type; DO NOT assume 39 paired historical stations.
    """
    from openpyxl import load_workbook
    a=load_workbook(s1_path,read_only=True,data_only=True)
    b=load_workbook(s2_path,read_only=True,data_only=True)
    sheet=next(iter(a.worksheets))
    rows=list(sheet.values)
    if len(rows)<800 or len(rows[4])<50:
        raise ValueError("Czech original S1 presence sheet unexpectedly changed")
    stations=[]
    last=None
    for col in range(2,len(rows[4])):
        code=rows[4][col]
        if code is not None and str(code).strip():
            last=str(code).strip()
        rec=rows[5][col]
        rec=str(rec).strip().lower() if rec is not None else None
        if last is not None and rec in ("rec","hist"):
            stations.append({"col_one_based":col+1,"locality_code":last,
                             "source_record_type":rec})
    categories=Counter()
    values=Counter()
    marker_by_type=Counter()
    presence_by_station=Counter()
    recent_alien_by_station=Counter()
    taxa=[]
    for row in rows[7:]:
        species=row[1] if len(row)>1 else None
        if species is None or not str(species).strip():
            continue
        taxa.append(str(species).strip())
        cat=str(row[0]).strip() if row[0] is not None else "(blank)"
        categories[cat]+=1
        for desc in stations:
            col=desc["col_one_based"]-1
            if col<len(row):
                v=row[col]
                if v is not None and str(v).strip():
                    values[str(v).strip().lower()]+=1
                    typ=desc["source_record_type"]
                    marker_by_type[(typ,str(v).strip().lower())]+=1
                    if typ=="rec":
                        presence_by_station[desc["locality_code"]]+=1
                        if "arch" in cat.lower() or "neo" in cat.lower():
                            recent_alien_by_station[desc["locality_code"]]+=1
    w2=next(iter(b.worksheets))
    w2rows=list(w2.values)
    s2codes=[str(x[1]).strip() for x in w2rows[5:]
             if len(x)>1 and x[1] is not None and str(x[1]).strip()]
    s2recent={}
    for x in w2rows[5:]:
        if len(x)<6 or x[1] is None or not str(x[1]).strip():
            continue
        code=str(x[1]).strip()
        s2recent[code]={"reported_taxa":int(x[3]),"reported_aliens":int(x[4]),
                        "reported_invasives":int(x[5]),
                        "source_total_area_ha":float(x[10]) if x[10] is not None else None}
    a.close()
    b.close()
    combos=sorted(set((x["locality_code"],x["source_record_type"]) for x in stations))
    return {
        "source_S1_row5_locality_by_merged_cells_row6_hist_rec":stations,
        "n_S1_station_x_record_type_columns":len(stations),
        "n_unique_S1_locality_codes":len(set(x["locality_code"] for x in stations)),
        "n_S1_union_taxon_rows_not_recent_only":len(taxa),
        "n_unique_S1_taxon_labels":len(set(taxa)),
        "n_repeated_taxon_labels":len(taxa)-len(set(taxa)),
        "original_source_S1_category_col_values":dict(categories.most_common()),
        "original_S1_presence_marker_types":dict(values.most_common()),
        "original_presence_marker_by_rec_hist":{typ+"__"+mark:n for (typ,mark),n in marker_by_type.most_common()},
        "n_recent_station_columns":sum(d["source_record_type"]=="rec" for d in stations),
        "n_historical_station_columns":sum(d["source_record_type"]=="hist" for d in stations),
        "station_s1_rec_any_marker_count_vs_s2_richness":{k:{"recent_marker_taxa":presence_by_station[k],"published_s2_taxa":s2recent[k]["reported_taxa"]} for k in sorted(s2recent)},
        "station_s1_rec_explicit_arch_neo_vs_s2_aliens":{k:{"recent_arch_neo_tagged":recent_alien_by_station[k],"published_s2_aliens":s2recent[k]["reported_aliens"]} for k in sorted(s2recent)},
        "n_original_S2_station_code_rows":len(s2codes),
        "n_original_S2_unique_station_codes":len(set(s2codes)),
        "S1_locality_codes_not_in_S2":sorted(set(x["locality_code"] for x in stations)-set(s2codes)),
        "S2_locality_codes_not_in_S1":sorted(set(s2codes)-set(x["locality_code"] for x in stations)),
        "source_observer_effort_historical_not_assumed_comparable":True,
    }


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
        "original_presence_table_semantics_audit":audit_presence_matrix_header_s1(s1_path,s2_path),
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
