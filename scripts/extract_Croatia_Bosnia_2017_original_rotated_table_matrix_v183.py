#!/usr/bin/env python3
"""Original-source 2017 Croatia/Bosnia railway stations: rotated vector Tab2.

Jasprica et al., Natura Croatica 26:271-303; doi:10.20302/NC.2017.26.21.
Recover only PDF-native source text (NOT OCR). Publisher survey of live track
areas at 15 stations. Source includes station numbers 1-15, life form and 
regional status; this is not experimental ballast equivalence across nations.
"""
import argparse,csv,hashlib,json,re
from pathlib import Path
PDF_URL="https://hrcak.srce.hr/file/282429"
SPECIES_RE=re.compile(r"^\*?[A-Z][A-Za-zÀ-ž\-]+\s+[a-z][a-zA-ZÀ-ž\-]+")
RECORD_RE=re.compile(r"^\s*\d{1,2}(?:\s*,\s*\d{1,2})*\s*,?\s*$")
LF={"T","H","P","Ch","G"}
FREQUENCIES={round(n*100/15) for n in range(1,16)}

def rows_from_original(pdf):
    import fitz
    raw=Path(pdf).read_bytes()
    if not raw.startswith(b"%PDF-"):raise ValueError("Not a PDF")
    doc=fitz.open(stream=raw,filetype="pdf")
    if len(doc)!=33:raise ValueError("Publisher PDF page count changed")
    out=[];page_reports=[]
    for pi in range(6,24):
        page=doc[pi]
        lines=[]
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines",[]):
                di=line.get("dir",(1,0))
                if abs(di[0])>.05 or abs(di[1]+1)>.05:continue
                txt="".join(sp.get("text","") for sp in line.get("spans",[])).strip()
                if not txt: continue
                x0,y0,x1,y1=line["bbox"]
                if x0<102 or x0>440:continue
                lines.append(dict(text=txt,x=(x0+x1)/2,y0=y0,y1=y1,x0=x0))
        taxon_anchors=[l for l in lines if l["y1"]>=560 and SPECIES_RE.match(l["text"])]
        page_reports.append({"page":pi+1,"n_taxon_anchor":len(taxon_anchors),
              "anchor_examples":[{"name":x["text"][:80],"x":round(x["x"],1)} for x in taxon_anchors[:4]]})
        for cand in taxon_anchors:
            # A row is a constant PDF x: all columns are vertically rotated.
            # Long species authority lines may use one extra glyph line.
            same=[l for l in lines if abs(l["x"]-cand["x"])<4.5]
            species=cand["text"].lstrip("*").strip()
            stationparts=[l["text"] for l in same if (400<=l["y0"] and l["y1"]<=507
                           and RECORD_RE.fullmatch(l["text"]))]
            rawstations=" ".join(stationparts).strip()
            if rawstations and not RECORD_RE.fullmatch(rawstations):
                rawstations=rawstations.replace(", ",",").replace(" ,",",").strip(",")
            vals=[int(v) for v in re.findall(r"\d{1,2}",rawstations)] if rawstations else []
            life=[l["text"] for l in same if 306<=l["y0"]<=333 and l["text"] in LF]
            family=[l["text"] for l in same if 344<=l["y0"]<=412 and re.fullmatch(r"[A-Za-z]+aceae",l["text"])]
            frequency=[int(l["text"]) for l in same if 100<=l["y0"]<=185 and l["text"].isdigit()
                       and int(l["text"]) in FREQUENCIES]
            out.append({"taxon_original":species,"page":pi+1,
                        "source_station_numbers":sorted(set(vals)),
                        "source_station_list_text":rawstations,
                        "life_form":life[0] if life else "",
                        "family":family[0] if family else "",
                        "frequency_percent_source":frequency[0] if frequency else None,
                        "source_row_x":round(cand["x"],2)})
    metadata={"original_DOI":"10.20302/NC.2017.26.21","source_pdf_url":PDF_URL,
       "publisher_original_pdf_sha256":hashlib.sha256(raw).hexdigest(),
       "publisher_original_pdf_bytes":len(raw),"source_n_pdf_pages":len(doc),
       "n_record_taxon_candidates":len(out),
       "n_unique_candidate_taxa":len(set(o["taxon_original"] for o in out)),
       "n_station_lists_nonempty":sum(bool(o["source_station_numbers"]) for o in out),
       "n_life_form_nonempty":sum(bool(o["life_form"]) for o in out),
       "n_family_nonempty":sum(bool(o["family"]) for o in out),
       "n_frequency_nonempty":sum(o["frequency_percent_source"] is not None for o in out),
       "first_source_row_geometry":[{"text":l["text"],"x":round(l["x"],1),"y0":round(l["y0"],1),"y1":round(l["y1"],1)} for l in lines if abs(l["x"]-taxon_anchors[0]["x"])<5][:18] if taxon_anchors else [],
       "n_frequency_disagreements":sum(o["frequency_percent_source"] is not None
             and round(len(o["source_station_numbers"])*100/15)!=o["frequency_percent_source"]
             for o in out),
       "malformed_station_ids":[{"taxon":o["taxon_original"],"records":o["source_station_numbers"]}
              for o in out if any(v<1 or v>15 for v in o["source_station_numbers"])][:15],
       "taxa_missing_station_preview":[{"taxon":o["taxon_original"],"page":o["page"],"x":o["source_row_x"]}
              for o in out if not o["source_station_numbers"]][:30],
       "page_reports":page_reports,
       "publication_expected_taxa":359,"publication_expected_stations":15,
       "source_count_gate":len(out)==359 and len(set(o["taxon_original"] for o in out))==359,
       "full_matrix_verified":False}
    metadata["full_matrix_verified"]=(metadata["source_count_gate"]
         and metadata["n_station_lists_nonempty"]==359
         and metadata["n_life_form_nonempty"]>=350
         and metadata["n_frequency_disagreements"]==0
         and not metadata["malformed_station_ids"])
    return out,metadata

def save(pdf,outdir):
    rows,a=rows_from_original(pdf)
    p=Path(outdir);p.mkdir(parents=True,exist_ok=True)
    (p/"Croatia_Bosnia_2017_native_rotated_table_source_QA.json").write_text(
         json.dumps(a,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({k:v for k,v in a.items() if k!="page_reports"},ensure_ascii=False,indent=2))
    print("PAGE REPORTS",json.dumps(a["page_reports"],ensure_ascii=False))
    if not a["full_matrix_verified"]:
         raise ValueError("Original vector PDF parsing did not pass 359-taxon station-list/trait frequency source gates")
    with (p/"Croatia_Bosnia_2017_15station_359taxa_source_traits.csv").open("w",newline="",encoding="utf8") as f:
         w=csv.writer(f);w.writerow(["taxon_original","station_codes","life_form","family","page"])
         for r in rows:w.writerow([r["taxon_original"],";".join(map(str,r["source_station_numbers"])),
                                    r["life_form"],r["family"],r["page"]])
    with (p/"Croatia_Bosnia_2017_15station_species_incidence.csv").open("w",newline="",encoding="utf8") as f:
         w=csv.writer(f);w.writerow(["study","rail_station_number","taxon_original","presence",
                                      "life_form","family","pdf_page"])
         for r in rows:
             for station in range(1,16):
                 w.writerow(["Jasprica2017",station,r["taxon_original"],
                             int(station in r["source_station_numbers"]),
                             r["life_form"],r["family"],r["page"]])

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original-pdf",required=True);p.add_argument("--output-dir",required=True)
    z=p.parse_args();save(z.original_pdf,z.output_dir)
