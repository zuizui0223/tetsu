#!/usr/bin/env python3
"""Inspect genuine Everingham et al 2019 raw Australian vegetation Appendix S4.

Original ecological source access only, no imputation of taxon nativeness
or display of original individual plot positions. Read via archive ZIP/XML.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
from collections import Counter
from audit_Mell_2024_railway_meta_original_xlsx_v317 import source_original_workbook

SOURCE_DOI="10.1111/aec.12719"
SOURCE_FILE="aec12719-sup-0002-AppendixS4.xlsx"
GENUINE_PUBLISHED_RAW_KIND="Appendix S4: Raw data of species collected during fieldwork"
ORIGINAL_META_CASE_IDS=["Everingham_2019_3","Everingham_2019_6"]

def audit(path):
    if path.name!=SOURCE_FILE:
        raise ValueError("Expected only the original official paper Appendix S4 name")
    sources,qa=source_original_workbook(path)
    if not 130000<=qa["bytes"]<=450000:
        raise ValueError("Original published ~191.2 kB XLSX not in plausible original byte range")
    sheets=[]
    for name,rows in sources.items():
        sample=[]
        all_text=Counter()
        for r in rows:
            for col,val in r.items():
                if isinstance(val,str):
                    v=val.strip()
                    if v and len(v)<80 and re.search(r"rail|fire trail|woodland|shoulder|lawns|roadside|native|exotic",v,re.I):
                        all_text[v]+=1
        for r in rows[:5]:
            sample.append([[col,str(value)[:85]] for col,value in sorted(r.items())[:35]])
        sheets.append({
            "name":name,"n_source_rows":len(rows),
            "max_source_column_index":max((max(r) for r in rows if r),default=-1),
            "first_5_source_rows_bounded_preview":sample,
            "field_label_matches_about_rail_fire_native_etc":dict(all_text.most_common(24))
        })
    if not sheets:raise ValueError("Source workbook has no plant survey worksheets")
    return {
        "version":"v318","original_publication_doi":SOURCE_DOI,
        "original_publisher_supplement_file":SOURCE_FILE,
        "original_publisher_supplement_public_description":GENUINE_PUBLISHED_RAW_KIND,
        "source_original_hash":qa,
        "worksheets":sheets,
        "original_meta_two_cases":ORIGINAL_META_CASE_IDS,
        "original_meta_railway_vs_firetrail_native_richness_Hedges_d":-1.63326905465261,
        "original_meta_railway_vs_firetrail_native_cover_Hedges_d":-0.516972129076364,
        "n_native_source_rail_firetrail_original_plots_expected_per_group":24,
        "original_species_by_railway_plot_reconstructed":False,
        "not_a_confirmatory_rail_habitat_selection_result":True,
        "no_raw_original_site_species_matrix_exported":True,
        "status":"ORIGINAL_2019_FIELDWORK_SOURCE_STRUCTURE_ONLY"
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    j=audit(a.source)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(j,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":j["status"],"original_sha":j["source_original_hash"]["sha256"],
                      "source_sheets":[{"name":x["name"],"rows":x["n_source_rows"],
                                       "columns":x["max_source_column_index"]+1}
                                       for x in j["worksheets"]]},ensure_ascii=False,indent=2))
if __name__=="__main__":main()
