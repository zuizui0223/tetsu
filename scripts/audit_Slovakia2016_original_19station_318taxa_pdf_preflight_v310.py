#!/usr/bin/env python3
"""v310 publisher original Slovak 2016 19-station FULL numeric incidence source audit.

Non-OCR use original author's published scientific PDF; no raw source PDF,
full list or original station×taxon incidence to GitHub. Source page semantics
and candidate entries only; no ecological train-number data is invented.
"""
from __future__ import annotations
import argparse,hashlib,json,re,unicodedata
from pathlib import Path
import fitz

SOURCE_URL="https://cbrb.sav.sk/usr/Janka/docs/BulletinSBS_2016_Majekova-Limanek.pdf"
DOI=None
ORIGINAL_PAPER="Májeková & Limánek (2016) Bulletin Slovenskej botanickej spoločnosti 38(1):37–45"
ORIGINAL_STATIONS=19
ORIGINAL_TAXA=318

def fold(s):
    v=unicodedata.normalize("NFKD",s)
    return "".join(c for c in v if not unicodedata.combining(c)).casefold()

def clean_text(t):
    t=t.replace("\u00ad","").replace("\xa0"," ").replace("\u200b","")
    t=re.sub(r"(?<=[a-zA-ZÀ-ž])-\s*\n\s*(?=[a-zA-ZÀ-ž])","",t)
    return t

def preflight(source:Path)->dict:
    buf=source.read_bytes()
    if not buf.startswith(b"%PDF") or len(buf)<25000:
        raise ValueError("Original author paper missing or non-PDF")
    doc=fitz.open(stream=buf,filetype="pdf")
    if doc.needs_pass or len(doc)!=9:raise ValueError("Publisher reported 9 printed pages 37-45")
    pages=[clean_text(p.get_text("text",sort=True)) for p in doc]
    doc.close()
    txt="\n".join(pages)
    flat=fold(txt)
    if "poprad" not in flat or "plavec" not in flat or "majekova" not in flat:
        raise ValueError("Downloaded PDF is not the verified 2016 Poprad–Plaveč paper")
    if not ("zistenych taxonov" in flat or "zoznam zistenych" in flat):
        raise ValueError("The raw numbered source species×station list cannot be found")
    start=flat.find("zoznam zistenych")
    headerpage=next((i+1 for i,p in enumerate(pages) if "zoznam zistenych" in fold(p)),None)
    pg=[]
    for i,t in enumerate(pages):
        # Raw page text has the original taxon list only from published p.41 onward.
        no=t.count(":")
        sem=t.count(";")
        r=re.findall(r"[A-Za-zÀ-ž\.\s×-]{3,100}:\s*(?:\d{1,2}(?:\s*[,–-]\s*\d{1,2})*)",t)
        pg.append({"pageno":i+1,"chars":len(t),"colons":no,"semicolons":sem,
                   "numeric_species_entry_candidates":len(r)})
    return {
       "author_original_paper":ORIGINAL_PAPER,
       "source_url":SOURCE_URL,"pdf_sha256":hashlib.sha256(buf).hexdigest(),
       "source_bytes":len(buf),"n_pdf_pages":len(pages),
       "n_original_recorded_taxa_according_to_publisher":ORIGINAL_TAXA,
       "n_numbered_original_station_sites_according_to_publisher":ORIGINAL_STATIONS,
       "source_station_numbers_are_not_train_frequencies":True,
       "numbered_taxon_incidence_list_found":start>=0,
       "source_incidence_list_first_page_index":headerpage,
       "all_original_page_lexical_profiles":pg,
       "source_true_2012_2014_train_volume_per_station":None,
       "source_original_2012_2014_station_x_taxon_matrix_NOT_YET_parsed":True,
       "this_stage_no_new_traffic_or_ecological_beta_claim":True
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--pdf",required=True,type=Path)
    a.add_argument("--out",required=True,type=Path)
    args=a.parse_args()
    r=preflight(args.pdf)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(r,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
