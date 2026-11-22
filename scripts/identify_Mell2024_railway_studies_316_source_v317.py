#!/usr/bin/env python3
"""Exact original systematic-review Railway study identities from Mell2024 S5.

The 316 studies are not interchangeable with 205 effect rows; this source
encodes all Q1–Q6 designs including railway floristics, management,
connectivity and invasive-plant comparisons.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from collections import Counter
from audit_Mell_2024_railway_meta_original_xlsx_v317 import source_original_workbook
SOURCE_SHA="85042f37984a6f8f316f1f578e7f3093334ded057ebdc0c72437ac8b34579c90"
REQUIRED=("ID_doc","Publication_author","Publication_year","Publication_title","DOI",
          "Question","Risk_of_bias","Country","LTI","LTI verge","Comparison",
          "Outcomes","Included_in_metaanalysis")

def inspect(path):
    original,qa=source_original_workbook(path)
    if qa["sha256"]!=SOURCE_SHA:
        raise ValueError("Original review file5 SHA mismatch; no unsupported railway cohort accepted")
    if "D - Systematic review table" not in original:
        raise ValueError("Missing original review's 316-study full inclusion sheet D")
    source=original["D - Systematic review table"]
    header=source[0]
    keys={str(value):i for i,value in header.items()}
    if any(k not in keys for k in REQUIRED):
        raise ValueError("Original review source LTI/outcomes/risk header changed")
    if len(source)!=317:
        raise ValueError("316 included source study rows not reproduced")
    lookup=[]
    for row in source[1:]:
        kind=str(row.get(keys["LTI"],""))
        habitat=str(row.get(keys["LTI verge"],""))
        by_LTI=any(token.strip().lower()=="railway" for token in kind.split("|"))
        by_verge=("railway" in habitat.casefold() or "railroad" in habitat.casefold())
        if by_LTI or by_verge:
            out={key:row.get(index) for key,index in keys.items()
                 if key not in ("GPS coordinates","Regions")}
            out["identified_in_explicit_source_LTI"]=by_LTI
            out["identified_in_source_verge_description"]=by_verge
            lookup.append(out)
    idlist=sorted({str(r["ID_doc"]) for r in lookup})
    base=[r for r in lookup if r["identified_in_explicit_source_LTI"]]
    explicit_unique=sorted({str(r["ID_doc"]) for r in base})
    quality=Counter(str(r.get("Risk_of_bias")) for r in lookup)
    questions=Counter(str(r.get("Question")) for r in lookup)
    return {
        "version":"v317_original_316_study_rail_id_map",
        "source_original_publisher_Additional5_SHA256":SOURCE_SHA,
        "source_file5_original_study_rows":len(source)-1,
        "original_review_railway_studies_reported":10,
        "n_original_study_rows_with_explicit_railway_LTI":len(base),
        "n_unique_original_review_document_IDs_railway_LTI":len(explicit_unique),
        "n_rows_with_railway_mentioned_in_LTI_or_verge":len(lookup),
        "n_unique_review_document_IDs_with_railway_anywhere":len(idlist),
        "question_fields_count":dict(questions),
        "original_source_risk_by_rail_document_rows":dict(quality),
        "source_study_rows":lookup,
        "author_publication_distinct_not_assumed_independent_plots":True,
        "comparability_caveat":"Original railway source Q1–Q6 studies use different comparators, biological response and controls; two are in 205-case quantitative Q2 meta",
        "no_railway_meta_effect_pooled":True,
        "no_raw_species_by_field_plot_or_coordinate_data_exported":True
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--xlsx",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    obj=inspect(a.xlsx)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
       "genuine_original_source_316_studies":obj["source_file5_original_study_rows"],
       "rail_study_unique_IDs":obj["n_unique_original_review_document_IDs_railway_LTI"],
       "rail_study_entries":obj["n_original_study_rows_with_explicit_railway_LTI"],
       "qualifying_source_cases":[{k:str(x.get(k,""))[:150] for k in ("ID_doc","Publication_title",
                        "Country","Question","Risk_of_bias","LTI","Comparison","Outcomes")}
                     for x in obj["source_study_rows"]]
    },indent=2,ensure_ascii=False))
if __name__=="__main__":main()
