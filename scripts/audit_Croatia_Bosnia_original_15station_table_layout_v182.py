#!/usr/bin/env python3
"""Source-first page geometry audit for Jasprica et al. 2017 railway flora.
Publisher PDF: https://hrcak.srce.hr/file/282429 ; DOI 10.20302/NC.2017.26.21.
Taxon/rail-station incidence resides in an unusually ROTATED table. This
audit does NOT assume sorted text is equivalent to a taxon×15-site matrix.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

def audit(original_pdf):
    import fitz
    raw=Path(original_pdf).read_bytes()
    if not raw.startswith(b"%PDF-"): raise ValueError("Missing original PDF")
    doc=fitz.open(stream=raw,filetype="pdf")
    out={"original_DOI":"10.20302/NC.2017.26.21",
         "original_publication":"Jasprica et al. Natura Croatica 26:271-303 (2017)",
         "original_pdf_url":"https://hrcak.srce.hr/file/282429",
         "original_pdf_sha256":hashlib.sha256(raw).hexdigest(),
         "original_pdf_bytes":len(raw),"n_pages":len(doc),"expected_n_station":15,
         "expected_n_taxa":359,"matrix_extracted":False,"pages":[]}
    for pi in (6,7,13,20,25,28):
        if pi>=len(doc):continue
        pg=doc[pi]
        blocks=pg.get_text("dict")["blocks"]
        specimens=[]
        directions={}
        for block in blocks:
            for line in block.get("lines",[]):
                di=tuple(round(float(v),2) for v in line.get("dir",[1,0]))
                key=str(di)
                directions[key]=directions.get(key,0)+1
                st="".join(sp.get("text","") for sp in line.get("spans",[]))
                if st.strip() and len(specimens)<22:
                    specimens.append({"text":st[:130],
                                      "bbox":[round(x,1) for x in line["bbox"]],
                                      "dir":di})
        out["pages"].append({"page_index":pi,"page_number":pi+1,
                             "width":round(pg.rect.width,1),"height":round(pg.rect.height,1),
                             "directions":directions,
                             "plain_text_start":pg.get_text("text")[:1600],
                             "sample_lines":specimens})
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--pdf",required=True)
    p.add_argument("--output",required=True)
    x=p.parse_args()
    out=audit(x.pdf)
    dest=Path(x.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(out,ensure_ascii=False,indent=2))
