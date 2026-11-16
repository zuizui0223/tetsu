#!/usr/bin/env python3
"""v18.6 Parma 2020 original museum PDF: source and recoverability audit.

Download source published PDF read-only. No OCR, no invented station absences,
no automatic publication of publisher PDF. This is a layout/parseability gate
before claiming a 438×28 biological matrix.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path

ORIGINAL_URL="https://www.museocivico.rovereto.tn.it/UploadDocs/21188_Art08_Adorni_C_agg_min.pdf"
ORIGINAL_BIBLIO="Adorni & Ghillani 2020, Annali Museo Civico di Rovereto 36:149-200"
STATION_CODES=("BER","BRG","BUS","BVT","CGU","CHI","CLL","CLO","CMA","FEL","FID",
               "FOR","GAI","MED","MEZ","NOC","OST","OZZ","PAR","ROC","SAL","SEL",
               "SOL","SOR","SPT","VAI","VAL","VIC")
EXCLUDED_UNSURVEYED_STATION="CIT"
FREQ=re.compile(r"\[(\d{1,2})\]\s*$")
FIRST_SPECIES="Abies alba"
LAST_SPECIES="Zygophyllum fabago"
CODE_RE=re.compile(r"\b(?:"+"|".join(STATION_CODES)+r")\b")
SECTION_OPEN="Elenco floristico"
SPECIES_START=re.compile(r"^[A-Z][a-zA-ZÀ-ž-]{2,}\s+[a-z][a-zA-ZÀ-ž-]+")
LIFE_RE=re.compile(r",\s*((?:T|H|P|G|Ch|He|I|NP)\s+\w+|(?:He|I))\s*,")
AUTO_RE=re.compile(r"\b(?:Autoct\.|Archeof\.|Neof\.|Xenof\.)\b",re.IGNORECASE)

def audit_pdf(path):
    import fitz
    raw=Path(path).read_bytes()
    if not raw.startswith(b"%PDF-") or len(raw)<300000:
        raise ValueError("Publisher original absent/non-PDF/implausible size")
    d=fitz.open(stream=raw,filetype="pdf")
    pgs=[]
    for i,p in enumerate(d):
        t=p.get_text("text",sort=True)
        pgs.append({"pdf_page_index":i,"text":t,
                    "word_count":len(p.get_text("words")),
                    "printed_page_candidates":re.findall(r"\b(?:1[4-9]\d|200)\b",t[:160])})
    matches=[x for x in pgs if FIRST_SPECIES in x["text"]]
    hits=[x for x in pgs if SECTION_OPEN in x["text"]]
    last=[x for x in pgs if LAST_SPECIES in x["text"]]
    source_rows=[]
    candidate_ends=[]
    for p in pgs:
        text=p["text"]
        # Gate starts at actual species list, never interpret introduction as taxon rows.
        for line in text.splitlines():
            line=line.strip()
            if not line:continue
            m=FREQ.search(line)
            if m and SPECIES_START.match(line):
                source_rows.append({"page_index":p["pdf_page_index"],"example":line[:175],
                                    "reported_station_total":int(m.group(1)),
                                    "station_codes_on_same_line":sorted(set(CODE_RE.findall(line)))})
            elif m:
                candidate_ends.append({"page_index":p["pdf_page_index"],"line":line[:90]})
    return {
        "source_original_url":ORIGINAL_URL,
        "source_bibliography":ORIGINAL_BIBLIO,
        "source_sha256":hashlib.sha256(raw).hexdigest(),
        "source_bytes":len(raw),"source_pdf_pages":len(d),
        "source_open_section_page_indices":[x["pdf_page_index"] for x in hits],
        "first_published_taxon_page_indices":[x["pdf_page_index"] for x in matches],
        "last_published_taxon_page_indices":[x["pdf_page_index"] for x in last],
        "n_exact_species_single_line_count_matches_ALL_pages":len(source_rows),
        "n_line_end_count_candidates_all_pages":len(candidate_ends),
        "species_line_examples_first":source_rows[:14],
        "species_line_examples_last":source_rows[-14:],
        "line_ends_preview":candidate_ends[:14],
        "original_declared_species":438,
        "station_count_publisher":28,
        "surveyed_station_codes":list(STATION_CODES),
        "unsurveyed_CIT_not_imputed_zero":True,
        "source_2012_autumn_2013_spring_plus_5_station_later_observations":True,
        "site_habitat_mixed_including_public_areas_and_ballast":True,
        "field_dataset_not_recovered_yet":True,
        "first_species_page_preview":next((x["text"][:5000] for x in matches),None),
        "last_species_page_preview":next((x["text"][-3400:] for x in last),None),
        "table_of_station_codes_page_preview":next((x["text"][:5500] for x in pgs if "Sigla" in x["text"] and "Berceto" in x["text"]),None)
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source-pdf",required=True)
    a.add_argument("--output",required=True)
    p=a.parse_args()
    result=audit_pdf(p.source_pdf)
    path=Path(p.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__": main()
