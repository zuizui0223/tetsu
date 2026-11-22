#!/usr/bin/env python3
"""v320 precise inspection of original FR/BE full-report six railway site chapters.

This identifies chapter page pointers, measured flora columns and long
Latin botanical names embedded in source page text. It deliberately does
NOT assert a true site-by-species vegetation matrix from incidental
mentions in tables, citations or landscape descriptions.

Publisher original 420-page report is acquired transiently in CI; only
source QA page pointers and counts are emitted.
"""
from __future__ import annotations
import argparse,hashlib,json,re,unicodedata
from pathlib import Path
from collections import Counter
import fitz

EXPECTED_SHA="25a1a38512156bf1264b697efd78984a1abf9217f1da21e6ef2c5623da509405"
SITES={
"17":("Esplechin","Paris-Bruxelles","TGV"),
"18":("Fretin","Paris-Lille","TGV"),
"19":("Chapelle d'Armentières","Lille-Londres","TGV"),
"20":("Phalempin","Lille-Douai","conventional"),
"21":("Lauwe","Lille-Kortrijk","conventional"),
"24":("Lompret","Lille-Dunkerque","conventional"),
}
RE_SCIENTIFIC=re.compile(r"(?<![A-Za-zÀ-ž])([A-Z][a-z]{3,})\s+([a-z]{3,})(?![A-Za-zÀ-ž])")
LEX={
"botanique":"botani","floristique":"florist","plant_species":"espèces",
"vegetation":"végéta","vascular":"vascul","species_listing":"espèces relevées",
"railway":"ferrovia","sample":"relevé","transect":"transect",
"cover":"recouvrement",
}
def normalize(s):
    return "".join(a for a in unicodedata.normalize("NFKD",s).lower()
                   if not unicodedata.combining(a))
def source_audit(path):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=EXPECTED_SHA:
        raise ValueError("Full genuine 2015 original PDF source changed")
    d=fitz.open(stream=raw,filetype="pdf")
    if len(d)!=420:raise ValueError("Original complete report no longer 420 pages")
    pages=[p.get_text("text",sort=True) for p in d]
    normals=[normalize(t) for t in pages]
    out={}
    for code,(place,line,rail_type) in SITES.items():
        token=normalize(place)
        hits=[k for k,p in enumerate(normals) if token in p]
        # Avoid table of contents/infrastructure overview. Main report's
        # indexed site chapters begin about printed page 100, while appendices
        # begin after 350. We do not infer plot site lists from all mentions.
        main=[k for k in hits if 85<=k<=355]
        report_keyword_pages=[]
        for k in main:
            observed={t:len(re.findall(normalize(lex),normals[k])) for t,lex in LEX.items()}
            if sum(observed.values())>=2:
                report_keyword_pages.append({
                    "pdf_page_1based":k+1,"flora_source_keyword_counts":observed,
                    "possible_Latin_binomial_text_glyph_count_UNVALIDATED":len(RE_SCIENTIFIC.findall(pages[k])),
                })
        out[code]={
            "original_physical_infrastructure_type":rail_type,
            "source_line_designation":line,
            "original_report_site_label_suffix":place,
            "n_full_report_pages_where_location_name_appears":len(hits),
            "main_site_chapter_name_found_on_pdf_pages_1based":[i+1 for i in main],
            "bounded_pages_with_flora_keywords_and_same_railway_site_name":report_keyword_pages[:24],
            "raw_source_full_species_by_site_directly_confirmed":False,
        }
    summary={
      "version":"v320_rail_original_chapter_pointer_audit",
      "original_complete_source_sha256":EXPECTED_SHA,
      "source_report_n_pages":420,
      "original_report_chapter_site_candidate_numbers":[17,18,19,20,21,24],
      "original_railway_site_names_chapter_source_support":out,
      "n_original_rail_candidate_sites_in_toc":6,
      "not_six_matched_natural_reference_pairs":True,
      "n_verified_railway_site_by_taxon_incidence":0,
      "no_true_matched_railway_vs_local_nonrail_effect_estimated":True,
      "not_inferred_thermal_niche_from_document_text":True,
      "safe_note":"Location-name matches and Latin binomial-like tokens DO NOT license per-site species incidences; verify actual author subplot field inventories and equal non-rail matched plot design before biology."
    }
    d.close()
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pdf",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    r=source_audit(a.pdf)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
      "source_sha":r["original_complete_source_sha256"],
      "chapter_site_page_reference":{k:{
        "where_site_name_main_report":v["main_site_chapter_name_found_on_pdf_pages_1based"],
        "source_biology_hits":v["bounded_pages_with_flora_keywords_and_same_railway_site_name"][:3]
      } for k,v in r["original_railway_site_names_chapter_source_support"].items()},
      "real_matched_nonrail_pairs_verified":0,
    },indent=2,ensure_ascii=False))
if __name__=="__main__":main()
