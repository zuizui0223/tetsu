#!/usr/bin/env python3
"""v11.2: exactly verified Dryad 2021 Kashmir Himalayan railway vegetation matrix.

Original: Rashid et al. J Appl Ecol 2021 10.1111/1365-2664.13961;
actual correct Dryad DOI 10.5061/dryad.0p2ngf20h verified with
NLM catalog AND Dryad API in real v11.1, version 99735, file 546003.
This audit only describes original workbook schema, not a fitted
within-rail-transect biological effect or linked wind-diaspore trait.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from pathlib import Path

SOURCE_DOI="10.5061/dryad.0p2ngf20h"
SOURCE_FILE_URL="https://datadryad.org/api/v2/files/546003/download"
SOURCE_SHA256="529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474"
SOURCE_N_BYTES=122187


def verify_original_file(path):
    data=Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest()!=SOURCE_SHA256 or len(data)!=SOURCE_N_BYTES:
        raise ValueError("Original 2021 Dryad railway vegetation Matrix_plot XLSX SHA/bytes differ")
    if not zipfile.is_zipfile(path):
        raise ValueError("Original Dryad matrix not OOXML")
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:
            raise ValueError("Original Dryad XLSX CRC failed")
    return True


def audit_original_source_xlsx(path):
    from openpyxl import load_workbook
    verify_original_file(path)
    w=load_workbook(path,read_only=True,data_only=True)
    sheets=[]
    for ws in w.worksheets:
        first=list(ws.iter_rows(min_row=1,max_row=min(ws.max_row,12),values_only=True))
        mincols=min(ws.max_column,35)
        def clean(x):
            if x is None:return None
            a=str(x).strip()
            return a[:85]
        sheets.append({
          "sheet_name":ws.title,
          "publisher_max_rows":ws.max_row,
          "publisher_max_columns":ws.max_column,
          "first_10_source_rows_first_35_columns":[
              [clean(v) for v in row[:mincols]] for row in first[:10]
          ],
          "first_row_nonempty_header_values_first_60":[
              clean(x) for x in first[0][:60] if x is not None] if first else [],
          "first_data_row_nonempty_values_first_60":[
              clean(x) for x in first[1][:60] if x is not None] if len(first)>1 else [],
        })
    w.close()
    return {
      "complete":True,"version":"v11.2",
      "original_published_paper_DOI":"10.1111/1365-2664.13961",
      "original_public_Dryad_DOI":SOURCE_DOI,
      "original_data_author_title":"Species abundance along the railway of Kashmir Himalaya",
      "source_Dryad_version_id":99735,"source_Dryad_file_id":546003,
      "download_url":SOURCE_FILE_URL,
      "original_source_SHA256":SOURCE_SHA256,
      "source_file_bytes":SOURCE_N_BYTES,
      "n_original_xlsx_sheets":len(sheets),
      "source_workbook_sheets":sheets,
      "no_species_wind_trait_category_yet":True,
      "no_imputed_2014_or_2017_repeat_sites":True,
      "not_a_train_seed_or_mowing_causality_estimate":True,
      "original_field_data_belongs_to_Rashid_and_colleagues":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    data=audit_original_source_xlsx(a.original_xlsx)
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(data,ensure_ascii=False),flush=True)


if __name__=="__main__":main()
