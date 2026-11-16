#!/usr/bin/env python3
"""v11.8 source-only audit of public Pladias Czech flora dispersal and life-form traits.

The published Slovenia railway 2026 field workbook uses original
Sádlo et al. 2018 eight/few strategy classes (Epilobium, Phragmites,
Allium ...); Pladias public page exposes the same source strategy
archetype feature, plus Raunkiaer life-form. Public website at
https://pladias.cz/en/download/features, original author citation:
Sádlo J., Chytrý M., Pergl J., Pyšek P. (2018) Preslia 90:1–22.
Official download URLs feature/141 (strategy) and /136 (life form).
Pladias (2021) DOI 10.23855/preslia.2021.001.

This is metadata/schema/coverage AUDIT only. Do NOT claim the original
railway survey has already been matched to these trait records. Do
not archive or republish raw provider trait files. Substantial use for
a publication entails contacting the Pladias Governing Board.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,zipfile
from collections import Counter
from pathlib import Path

DISP_URL="https://pladias.ibot.cas.cz/public/traits/downloadTraitData/feature/141/lang/en"
LIFE_URL="https://pladias.ibot.cas.cz/public/traits/downloadTraitData/feature/136/lang/en"
DISP_CODES=("Allium","Bidens","Cornus","Epilobium","Lycopodium",
            "Phragmites","Sparganium","Wolffia","Zea")

def csv_rows(raw):
    txt=None;enc=None
    for codepage in ("utf-8-sig","utf-16","cp1250"):
        try:
            value=raw.decode(codepage)
        except UnicodeDecodeError:continue
        # utf-16 should not be selected for ordinary ascii from arbitrary pairs.
        if "\x00" in value[:100]:continue
        txt=value;enc=codepage;break
    if txt is None:raise ValueError("Unrecognized authentic public source character encoding")
    first=txt[:1000].lower().lstrip()
    if first.startswith(("<!doctype","<html","<?xml")):
        raise ValueError("Publisher endpoint returned HTML/XML not trait table")
    lines=txt.splitlines()
    if len(lines)<10:raise ValueError("Public trait download has too few rows")
    first_line=lines[0]
    try: dialect=csv.Sniffer().sniff("\n".join(lines[:15]),delimiters=",;\t|")
    except csv.Error:
        dialect=csv.excel_tab if "\t" in first_line else csv.excel
    reader=csv.DictReader(io.StringIO(txt),dialect=dialect)
    fields=reader.fieldnames or []
    data=[dict(r) for r in reader]
    if len(fields)<2 or len(data)<50:
        raise ValueError("Public Pladias source not a useful species trait table")
    return enc,fields,data


def xlsx_rows(raw):
    from openpyxl import load_workbook
    if not zipfile.is_zipfile(io.BytesIO(raw)):raise ValueError("Source ZIP not XLSX")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip() is not None:raise ValueError("Pladias office spreadsheet ZIP CRC failed")
        if not any(s.startswith("xl/workbook.xml") for s in z.namelist()):
            raise ValueError("Archive is not a source XLSX spreadsheet")
    wb=load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    ws=wb.worksheets[0];it=ws.values
    headers=next(it)
    fields=[str(x).strip() if x is not None else "column_"+str(i) for i,x in enumerate(headers)]
    rows=[dict(zip(fields,row)) for row in it]
    wb.close()
    if len(rows)<50:raise ValueError("Published Pladias source XLSX has too few taxa")
    return "xlsx",fields,rows


def inspect_source(path,trait):
    raw=Path(path).read_bytes()
    if len(raw)<200:raise ValueError("Pladias source too small, probably source access error")
    kind="zip-office" if raw[:2]==b"PK" else "csv-text"
    enc,headers,rows=xlsx_rows(raw) if kind=="zip-office" else csv_rows(raw)
    cats=Counter()
    nonempty=Counter()
    for row in rows:
        for k,v in row.items():
            if v is None:continue
            s=str(v).strip()
            if not s or s.lower() in ("na","nan","null"):continue
            nonempty[k]+=1
            for cat in DISP_CODES:
                if re.search(r"\b"+re.escape(cat)+r"\b",s,re.I):
                    cats[cat]+=1
    return {
        "trait_source":trait,"original_source_bytes":len(raw),
        "source_file_sha256":hashlib.sha256(raw).hexdigest(),
        "file_kind":kind,"parsed_encoding":enc,
        "fieldnames":headers,
        "n_original_species_trait_source_rows":len(rows),
        "nonempty_columns":dict(nonempty),
        "strategy_archetype_token_counts_all_fields":dict(cats),
        "first_four_original_field_rows_for_schema_only":[{
            k:(str(v)[:90] if v is not None else None) for k,v in row.items()}
            for row in rows[:4]],
        "original_download_was_not_CzechVeg_railway_species_join":True,
    }


def run(disp_file,life_file):
    d=inspect_source(disp_file,"published_Pladias_diaspore_strategy")
    l=inspect_source(life_file,"published_Pladias_life_form")
    if d["source_file_sha256"]==l["source_file_sha256"]:
        raise ValueError("Same published traits file returned for two distinct source features")
    return {
        "complete":True,"version":"v11.8",
        "Pladias_database_DOI":"10.23855/preslia.2021.001",
        "Pladias_publisher_feature_page":"https://pladias.cz/en/download/features",
        "SADLO2018_strategies_source":"Sádlo et al. (2018) Preslia 90:1-22",
        "strategy_endpoint":DISP_URL,"lifeform_endpoint":LIFE_URL,
        "strategy_original_source_metadata":d,
        "lifeform_original_source_metadata":l,
        "no_original_CzechVeg_species_names_merged":True,
        "no_Czech_rail_trait_effect_fitted":True,
        "source_rights":"Public Pladias data: source citation required. Contact database Governing Board to discuss collaboration if substantial part of a publication.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--dispersal-source",required=True)
    p.add_argument("--life-form-source",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    report=run(a.dispersal_source,a.life_form_source)
    o=Path(a.output);o.parent.mkdir(parents=True,exist_ok=True)
    o.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False),flush=True)


if __name__=="__main__":main()
