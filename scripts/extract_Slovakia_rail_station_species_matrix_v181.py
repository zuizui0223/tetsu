#!/usr/bin/env python3
"""Reconstruct a SOURCE-TRACEABLE Slovak railway station incidence matrix.

Májeková, Jehlík & Zaliberová (2016), Thaiszia 26:173-188, Tab. 2.
Original survey 2012–2013 and INCOMPLETE historical records 1968–1974.
Two detailed stations only. Parser uses original glyph coordinates, not OCR;
PDF is NOT committed or rehosted. Do not call station incidence ballast-only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

PDF_URL = "https://www.upjs.sk/app/uploads/sites/15/2022/10/173-188_Majekova_et_al-upr.pdf"
FALLBACK_URL = "https://cbrb.sav.sk/usr/Janka/docs/Majekova_et_al_Thaiszia_2016.pdf"
DOI = "Thaiszia 26(2):173-188 (2016)"
STATIONS = ("Kosice", "Spisska_Nova_Ves")
PERIODS = ("historical_1968_1974_incomplete", "recent_2012_2013")
EXPECTED = {"Kosice_historical_1968_1974_incomplete":72,
            "Kosice_recent_2012_2013":192,
            "Spisska_Nova_Ves_historical_1968_1974_incomplete":56,
            "Spisska_Nova_Ves_recent_2012_2013":191}
UNION_EXPECTED = {"Kosice":207, "Spisska_Nova_Ves":204}
ABUNDANCE = re.compile(r"^(?:x|[abcde](?:[-–][abcde])?)$", re.IGNORECASE)
GENUS = re.compile(r"^[A-Z][a-zA-ZÀ-ž-]{2,}$")
SPECIES = re.compile(r"^(?:[a-z][a-zA-ZÀ-ž-]+|cf\.|sect\.|sp\.)$")
PAGE_INDICES = range(4, 11)  # printed 177–183, source 173 starts at 0


def merge_words_to_rows(words, tolerance=3.0):
    """Return visual text lines; retain each word's physical x centre."""
    groups = []
    for w in sorted(words, key=lambda w: (round((w[1]+w[3])/2, 1), w[0])):
        x0, y0, x1, y1, val = w[:5]
        ym = (y0 + y1) / 2
        if not groups or abs(groups[-1][0]-ym) > tolerance:
            groups.append([ym, [w]])
        else:
            groups[-1][1].append(w)
            groups[-1][0] = (groups[-1][0] + ym) / 2
    return [(y, sorted(ws, key=lambda w: w[0])) for y, ws in groups]


def detect_header_centres(lines):
    """Find four x positions from publisher's repeated historic/recent headers."""
    best = None
    for y, line in lines:
        hdr = [((w[0]+w[2])/2, w[4].strip().lower()) for w in line
               if w[4].strip().lower().strip(".,") in ("historic", "recent")]
        hdr.sort()
        if len(hdr) == 4:
            centres = [x for x, _ in hdr]
            if [z for _, z in hdr] == ["historic", "recent", "historic", "recent"]:
                best = centres
                break
    return best


def parse_table_page(page, centres):
    lines = merge_words_to_rows(page.get_text("words", sort=False))
    out, diagnostics = [], []
    if not centres:
        return out, [{"error":"missing original source header centres"}]
    first_center = centres[0]
    name_end = first_center - 20  # header x, name-column right margin
    for ym, words in lines:
        left = [w for w in words if (w[0]+w[2])/2 < name_end]
        vals = [w for w in words if (w[0]+w[2])/2 >= name_end]
        if not left:
            continue
        parts = [str(w[4]).strip() for w in left]
        if len(parts)<2 or not GENUS.fullmatch(parts[0]) or not SPECIES.fullmatch(parts[1]):
            continue
        name = " ".join(parts).strip()
        if len(name)>190 or name.lower().startswith(("tab.", "railway", "thaiszia")):
            continue
        slots=["","","",""]
        for w in vals:
            token=str(w[4]).strip().replace("–","-")
            if not ABUNDANCE.fullmatch(token):
                continue
            xc=(w[0]+w[2])/2
            k=min(range(4),key=lambda i:abs(centres[i]-xc))
            if abs(xc-centres[k])>27:
                diagnostics.append({"page":page.number+1,"taxon":name,
                                    "ambiguous_token":token,"x":round(xc,1)})
                continue
            if slots[k] and slots[k]!=token:
                diagnostics.append({"page":page.number+1,"taxon":name,
                                    "conflicting_abundance":[slots[k],token]})
            slots[k]=token
        if any(slots):
            out.append({"taxon_original":name,"page_pdf_index":page.number,
                        "values":slots})
        else:
            diagnostics.append({"page":page.number+1,
                                "candidate_without_presence":name[:90]})
    return out, diagnostics


def extract(pdf_path):
    import fitz
    raw=Path(pdf_path).read_bytes()
    if not raw.startswith(b"%PDF-") or len(raw)<40000:
        raise ValueError("Missing or implausibly short publisher PDF")
    f=fitz.open(stream=raw,filetype="pdf")
    if len(f)!=16:
        raise ValueError(f"Publisher PDF pages changed, expected 16 got {len(f)}")
    bypage=[merge_words_to_rows(f[i].get_text("words",sort=False)) for i in PAGE_INDICES]
    centers=[detect_header_centres(lines) for lines in bypage]
    header0=next((c for c in centers if c), None)
    if not header0:
        raise ValueError("Publisher table header historic/recent positions not detected")
    extracted,diag=[],[]
    for i in PAGE_INDICES:
        entries, errors=parse_table_page(f[i],centers[i-4] or header0)
        extracted.extend(entries);diag.extend(errors)
    totals=[0,0,0,0]
    unique=set()
    pairs=set()
    for row in extracted:
        name=row["taxon_original"]
        if name in unique:
            diag.append({"duplicate_taxon":name})
        unique.add(name)
        for j,v in enumerate(row["values"]):
            if v:
                totals[j]+=1
                pairs.add((name,j))
    unions={"Kosice":len(set(r["taxon_original"] for r in extracted if r["values"][0] or r["values"][1])),
            "Spisska_Nova_Ves":len(set(r["taxon_original"] for r in extracted if r["values"][2] or r["values"][3]))}
    actual={
        "Kosice_historical_1968_1974_incomplete":totals[0],
        "Kosice_recent_2012_2013":totals[1],
        "Spisska_Nova_Ves_historical_1968_1974_incomplete":totals[2],
        "Spisska_Nova_Ves_recent_2012_2013":totals[3]}
    return extracted,{
        "source":DOI, "source_pdf_url":PDF_URL,
        "source_pdf_sha256":hashlib.sha256(raw).hexdigest(),
        "source_pdf_bytes":len(raw),"source_pdf_page_count":len(f),
        "pdf_table_printed_pages":"177-183",
        "detected_header_x_centres":centers, "n_taxa_all":len(unique),
        "sample_species":["Echium vulgare","Amaranthus retroflexus","Bromus tectorum",
                          "Arabidopsis thaliana","Poa compressa"],
        "sample_entries":{x: next((r["values"] for r in extracted if r["taxon_original"]==x),None)
                          for x in ["Echium vulgare","Amaranthus retroflexus","Bromus tectorum",
                                    "Arabidopsis thaliana","Poa compressa"]},
        "observed_period_counts":actual,"published_expected_period_counts":EXPECTED,
        "observed_union_counts":unions,"published_union_counts":UNION_EXPECTED,
        "n_candidate_rejects":len(diag),"candidate_rejects_preview":diag[:45],
        "exact_source_count_gate":actual==EXPECTED and unions==UNION_EXPECTED,
        "historical_comparison_valid_for_change":False,
        "pure_ballast_habitat_known":False,
        "original_absence_assumed_complete_only_for_contemporary_surveys":True
    }


def write_outputs(records, audit, outdir):
    root=Path(outdir);root.mkdir(parents=True,exist_ok=True)
    (root/"Slovakia_2016_original_table2_provenance_audit.json").write_text(
        json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf8")
    if not audit["exact_source_count_gate"]:
        raise ValueError("Source counts fail published modern/historical/union totals. "
                         "Only provenance/diagnostics saved; no claimed verified incidence.")
    with (root/"Slovakia_2016_original_table2_4station_period_taxon_incidence.csv").open("w",newline="",encoding="utf8") as h:
        writer=csv.writer(h)
        writer.writerow(["country","original_study","station","period","taxon_original",
                         "source_abundance_code","presence","source_pdf_page"])
        for r in records:
            for j,v in enumerate(r["values"]):
                writer.writerow(["SK","Majekova_et_al_2016",STATIONS[j//2],
                                 PERIODS[j%2],r["taxon_original"],v,int(bool(v)),
                                 r["page_pdf_index"]+1])
    with (root/"Slovakia_2016_original_table2_contemporary_only_presence.csv").open("w",newline="",encoding="utf8") as h:
        writer=csv.writer(h)
        writer.writerow(["station","taxon_original","presence","abundance_code"])
        for r in records:
            for j,st in [(1,"Kosice"),(3,"Spisska_Nova_Ves")]:
                writer.writerow([st,r["taxon_original"],int(bool(r["values"][j])),r["values"][j]])


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-pdf",required=True)
    p.add_argument("--output-dir",required=True)
    args=p.parse_args()
    rows,audit=extract(args.source_pdf)
    write_outputs(rows,audit,args.output_dir)
    print(json.dumps({k:v for k,v in audit.items() if k not in ("candidate_rejects_preview","sample_species")},
                     ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    main()
