#!/usr/bin/env python3
"""Find independent railway field studies in Mell 2024 all-inclusion metadata.

Additional file 5 carries review-screening and validity tables beyond the
47 primary studies that supplied 205 numerical Q2 meta cases. No raw plant
data are reconstructed or invented here. Original source XLSX stays
ephemeral. Count source studies rather than effect rows.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from collections import Counter
from pathlib import Path
from audit_Mell_2024_railway_meta_original_xlsx_v317 import source_original_workbook,MELL5
EXACT_PUBLICATION_DOI="10.1186/s13750-024-00328-3"

def audit(path):
    if path.name!=MELL5:raise ValueError("Not the original additional file 5")
    sheet,qa=source_original_workbook(path)
    if qa["bytes"]<250000:raise ValueError("Source publication metadata XLSX unexpectedly small")
    results=[]
    for name,rows in sheet.items():
        hits=[]
        for i,row in enumerate(rows):
            cells=[[col,str(v)[:130]] for col,v in sorted(row.items()) if
                   isinstance(v,str) and re.search(r"\brail(?:way|road|track|ways)?\b",v,re.I)]
            if cells:
                hits.append({"row_1based":i+1,
                             "literal_source_cells":cells[:8],
                             "other_original_first_10_cells":[[col,str(v)[:75]] for col,v in sorted(row.items())[:10]]})
        results.append({
            "original_sheet":name,
            "n_source_rows":len(rows),"n_source_columns":max((max(r) for r in rows if r),default=-1)+1,
            "first_three_publisher_rows_bounded":[[[k,str(v)[:85]] for k,v in sorted(r.items())[:32]] for r in rows[:3]],
            "n_row_mentions_rail_any_field":len(hits),
            "bounded_railway_mentioned_source_rows":hits[:65]
        })
    return {
        "version":"v317_source_all_studies_additional_file5",
        "official_supplementary_filename":MELL5,
        "original_publication_doi":EXACT_PUBLICATION_DOI,
        "original_XLSX_sha_and_structure":qa,
        "source_all_study_sheet_structure":results,
        "number_of_all_included_studies_reported":316,
        "source_screening_eligibility_not_meta_effects":True,
        "original_railway_abstract_or_metadata_mentions_are_candidates_not_railway_meta_cases":True,
        "source_taxon_site_or_GPS_exported":False,
    }
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--xlsx",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    x=a.parse_args()
    j=audit(x.xlsx)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(j,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"sha256":j["original_XLSX_sha_and_structure"]["sha256"],
                       "sheets":[{"name":z["original_sheet"],"rows":z["n_source_rows"],
                                  "cols":z["n_source_columns"],"rail_keyword_rows":z["n_row_mentions_rail_any_field"]}
                                  for z in j["source_all_study_sheet_structure"]]},indent=2,ensure_ascii=False))
if __name__=="__main__":main()
