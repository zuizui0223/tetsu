#!/usr/bin/env python3
"""v12.8.1 DIRECT original-Austrian-government-PDF Table 14 machine audit.

The v12.8 original Pascher, Hainz-Renetzeder & Schneeweiss (2016)
Austrian Federal Ministry of Health 107-page PDF report was actually
downloaded successfully on GitHub Actions source run 43588505410:
raw ORIGINAL PDF 18,158,388 bytes, SHA256
52a1910fb24cc9a38da841322ce28e9c97b51758a632146551981edac696125d.
This script reads the original PDF itself via pdftotext -layout,
isolates published Table 14 (60 original survey sites × two years),
re-extracts EACH source site code, 2014/2015 number of sampled FERAL
Brassica napus PLANTS and checks all exact 60 original code×year
numerical pairs against separately manually transcribed v12.8 input.

If PDF extraction fails or differs on ANY site we must explicitly
fail; no invented "independently verified" original PDF results.
Caution: Table 14 numeric counts are plants sampled for leaf/DNA
analysis, NOT caught propagule seeds nor source population census.
This is published AUTHORS' original field study, not new tetsu data.
"""
from __future__ import annotations
import argparse,hashlib,json,re,subprocess
from collections import defaultdict
from pathlib import Path
from analyze_Austria_ferreral_osr_actual_sites_v128 import (
    SOURCE_PDF,canonical_source_codes,validate_original_table,
    observed_count_summary
)

ORIGINAL_GOV_PDF_BYTES=18158388
ORIGINAL_GOV_PDF_SHA256="52a1910fb24cc9a38da841322ce28e9c97b51758a632146551981edac696125d"
SITE_CODE=re.compile(r"^\s*((?:BHR|BHK|GBH|FB|DH|SR|SK|OEM)0?\d{1,2})\b",re.I)
NUMERIC_PAIR=re.compile(r"(?<!\d)(\d{1,5})\s+(\d{1,5})\s*$")
ORIGINAL_SITE_KEYS=canonical_source_codes()


def audit_actual_gov_pdf(original_pdf):
    raw=Path(original_pdf).read_bytes()
    h=hashlib.sha256(raw).hexdigest()
    if (len(raw)!=ORIGINAL_GOV_PDF_BYTES or h!=ORIGINAL_GOV_PDF_SHA256
        or not raw.startswith(b"%PDF-")):
        raise ValueError("Government ORIGINAL Austrian 2016 PDF bytes/SHA256 changed; no fake provenance")
    p=subprocess.run(["pdftotext","-layout",str(original_pdf),"-"],
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                     check=True,timeout=120)
    if not p.stdout.strip():
        raise ValueError("Original Austrian PDF yielded no machine-readable text")
    txt=p.stdout.decode("utf-8",errors="strict")
    return txt,{"original_gov_source_report_URL":SOURCE_PDF,
                "original_gov_pdf_SHA256":h,
                "original_gov_pdf_size_bytes":len(raw),
                "native_pdf_text_extracted_without_OCR":True,
                "original_full_report_sha_pinned_from_real_v12_8_gitHub_source_run":True}


def isolate_table14(txt):
    # Original publication page 53/54. Use published table heading,
    # not other mentions in TOC/Results and not site codes in Figure maps.
    a=re.search(r"Table\s*14\s*:\s*Sample\s+sites\s+and\s+number\s+of\s+sampled\s+feral",txt,re.I)
    if not a:
        # Punctuation may differ in text layout; but do not accept arbitrary context.
        a=re.search(r"Table\s*14\s*:\s*Sample\s+sites",txt,re.I)
    if not a:raise ValueError("Original 2016 government PDF Table14 anchored heading not found")
    rest=txt[a.end():]
    end_candidates=[]
    for pattern in (r"In\s+total\s*,?\s*2[,]?113\s+feral",
                    r"Figure\s*24\s*:",r"9\.2\.\s*Occurrence\s+of\s+hybridisation"):
        q=re.search(pattern,rest,re.I)
        if q:end_candidates.append(q.start())
    end=min(end_candidates) if end_candidates else len(rest)
    segment=rest[:end]
    if len(segment)<1500:
        raise ValueError("Original published Table14 isolated segment unexpectedly short")
    return segment


def extract_original_pdf_table_pairs(table14_text):
    records={}
    lines=table14_text.splitlines()
    duplicates=[]
    for i,line in enumerate(lines):
        code_search=SITE_CODE.search(line)
        if not code_search:continue
        code=code_search.group(1).upper()
        if code not in ORIGINAL_SITE_KEYS:continue
        # All Table14 source observations are printed as one row of site
        # code + name + regional label + two numeric year columns.
        # Retain up to two consecutive physical text lines only if
        # original PDF rendering wrapped the site name.
        joined=line
        numbers=NUMERIC_PAIR.search(joined)
        if numbers is None:
            for j in (i+1,i+2):
                if j>=len(lines):break
                if SITE_CODE.search(lines[j]):break
                joined=joined.rstrip()+" "+lines[j].strip()
                numbers=NUMERIC_PAIR.search(joined)
                if numbers is not None:break
        if numbers is None:
            raise ValueError("Actual original government Table14 source site "+code+
                             " did not contain two source year numeric columns")
        current=(int(numbers.group(1)),int(numbers.group(2)))
        if code in records:
            duplicates.append(code)
        records[code]=current
    if duplicates:
        raise ValueError("Original PDF table includes repeated site IDs: "+str(duplicates))
    if set(records)!=ORIGINAL_SITE_KEYS:
        missing=sorted(ORIGINAL_SITE_KEYS-set(records))
        extra=sorted(set(records)-ORIGINAL_SITE_KEYS)
        raise ValueError("Cannot faithfully recover all original 60 PDF Table14 site-year rows: "+
                         "missing="+str(missing[:15])+",extra="+str(extra[:15]))
    if sum(x[0] for x in records.values())!=1001 or sum(x[1] for x in records.values())!=1112:
        raise ValueError("Source PDF derived original 2014+2015 sampled plant totals disagree")
    if sum(a+b>0 for a,b in records.values())!=44:
        raise ValueError("Source PDF 44/60 feral sampled plant source-site detection mismatch")
    return records


def compare_original_pdf_to_manual_csv(author_pdf,transcribed_csv):
    text,source=audit_actual_gov_pdf(author_pdf)
    table=isolate_table14(text)
    original=extract_original_pdf_table_pairs(table)
    manual=validate_original_table(transcribed_csv)
    mismatches=[]
    for site in manual:
        key=site["source_site_code"]
        expected=(site["sampled_2014"],site["sampled_2015"])
        if original[key]!=expected:
            mismatches.append({"original_site":key,
                               "gov_report_pdf_year_counts":original[key],
                               "manually_transcribed_year_counts":expected})
    if mismatches:
        raise ValueError("Manual original table source rows differ from government original PDF; "+
                         json.dumps(mismatches[:10]))
    return {
        "complete":True,"version":"v12.8.1",
        "original_authentic_government_pdf_source":source,
        "n_site_codes_verified_individually_against_original_PDF":len(original),
        "n_original_PDF_year_by_site_numeric_cells_verified":len(original)*2,
        "n_manual_transcription_mismatches":len(mismatches),
        "original_PDF_author_table14_exact_source_sums":{
             "original_sites":60,"original_2014_sampled_plants":1001,
             "original_2015_sampled_plants":1112,
             "original_site_sampled_OSR_in_at_least_one_year":44},
        "original_report_measured_sAMPLED_FERAL_PLANTS_NOT_ACTUAL_SEED_SPILLAGE":True,
        "not_rail_freight_causal_hypothesis_test":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gov-pdf",required=True)
    p.add_argument("--manually-transcribed-v128",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=compare_original_pdf_to_manual_csv(a.gov_pdf,a.manually_transcribed_v128)
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)


if __name__=="__main__":main()
