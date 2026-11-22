#!/usr/bin/env python3
"""Original report 2015 actual railway site 17–21,24 flora-section feasibility.

Read complete original 420-page site chapters, not incidental names across
the territory or the management summary. A printed narrative with some Latin
binomials does NOT create a complete source species-by-station inventory.
"""
from __future__ import annotations
import argparse,hashlib,json,re,unicodedata
from pathlib import Path
import fitz

SOURCE_SHA="25a1a38512156bf1264b697efd78984a1abf9217f1da21e6ef2c5623da509405"
FLORA_PATTERNS=("la flore","synthèse floristique","relevé floristique",
                "relevés botaniques","placettes du botaniste",
                "infrastructure / matrice","nombre d'espèces","richesse floristique")
LATIN=re.compile(r"\b[A-Z][a-z]{3,}\s+[a-z]{3,}\b")
# These are printed SOURCE chapter starts, not invented sampling transect
# coordinates. 24 is the last site before synthesis at printed p345.
RAIL_BLOCKS=[(17,263,269),(18,269,276),(19,276,284),
             (20,284,292),(21,292,300),(24,324,345)]

def fold(s):
    return "".join(c for c in unicodedata.normalize("NFKD",s).lower()
                   if not unicodedata.combining(c))
def audit(path):
    buf=path.read_bytes()
    if hashlib.sha256(buf).hexdigest()!=SOURCE_SHA:
        raise ValueError("Publisher genuine 2015 full original PDF bytes changed")
    d=fitz.open(stream=buf,filetype="pdf")
    if len(d)!=420:raise ValueError("Original report site chapter structure changed")
    out={}
    for code,lo,hi in RAIL_BLOCKS:
        if not (100<lo<hi<350):
            raise ValueError("Original chapter range not within botanical site part")
        pages=[d[i-1].get_text("text",sort=True) for i in range(lo,hi)]
        counts={}
        match_pages={}
        for key in FLORA_PATTERNS:
            xs=[lo+j for j,t in enumerate(pages) if fold(key) in fold(t)]
            counts[key]=len(xs)
            match_pages[key]=xs
        # Genuine complete inventories would require a source species list
        # with one uniquely keyed botanical-relevé site and a trusted end,
        # which mere Latin two-word tokens never demonstrate.
        plant_names=sorted({m.group(0) for t in pages for m in LATIN.finditer(t)})
        botanical_heading_page=next((lo+i for i,t in enumerate(pages)
                 if re.search(r"(?im)^\s*(?:la flore|flore|synthèse floristique)\s*$",fold(t))),None)
        out[str(code)]={
            "source_report_railway_site_code":code,
            "original_pdf_chapter_page_interval_1based":[lo,hi-1],
            "n_pdf_pages_inspected":len(pages),
            "flora_and_plant_reference_page_counts":counts,
            "source_report_page_refs_by_flora_topic":match_pages,
            "has_standalone_floristic_narrative_heading":botanical_heading_page is not None,
            "source_original_floristic_narrative_heading_pageno_1based":botanical_heading_page,
            "possible_binomial_pattern_distinct_text_strings_UNVALIDATED":len(plant_names),
            "is_full_original_species_by_site_list_proven":False,
            "same_source_local_nature_reference_plant_species_matrix_proven":False
        }
    d.close()
    return {
        "version":"v320_biology_section_original_full_report",
        "source_report_full_pdf_SHA256":SOURCE_SHA,
        "source_report_original_railway_site_chapter_count_checked":6,
        "source_site_chapter_plant_sections":out,
        "author_methods_describe_infrastructure_matrix_botanist_plots_WHERE_POSSIBLE":True,
        "source_limitation":"Maps have legend for botanical plots infrastructure/matrix and narrative floristics; this doesn't prove accessible complete paired site×taxa matrix and equal plot grain per railway site.",
        "eligible_original_rail_vs_local_nature_plant_plot_pairs":0,
        "new_species_traits_or_effect_estimates_approved":0,
        "do_not_call_0_verified_pairs_0_actual_pairs_in_source":True
    }
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pdf",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    x=audit(a.pdf)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
       "source_rail_sites_audited":x["source_report_original_railway_site_chapter_count_checked"],
       "flora_sections":{k:{
          "pages":v["source_report_page_refs_by_flora_topic"]["la flore"],
          "matched_plant_site_matrix_verified":False,
          "main_chapter_pages":v["original_pdf_chapter_page_interval_1based"]}
          for k,v in x["source_site_chapter_plant_sections"].items()},
       "all_true_rail_vs_matched_nature_incidence_columns_verified":0
    },indent=2,ensure_ascii=False))
if __name__=="__main__":main()
