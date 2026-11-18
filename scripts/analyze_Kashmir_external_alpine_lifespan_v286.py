#!/usr/bin/env python3
"""v28.6: Crosswalk genuinely independent Kashmir flora life spans to Dryad rail plots.

Only same taxon names (exact text or unambiguous species binomial) are
included. Generic/family inference, fabricated trait values, and any
2021/2017 individual plot data redistribution are prohibited.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import random
import re
import statistics
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_Kashmir_annual_microsites_v285 import (
    original_taxon_headers, annual_stat, summarize, life_span_crosswalk,
)
from analyze_Kashmir_same_sites_reference_year_v254 import (
    DRYAD_SHA256, NS, PKG_NS, REL_NS, parse_xlsx_matrix,
)

ALPINE_DOI = "https://doi.org/10.5061/dryad.vmcvdnd4c"
RAIL_DOI = "https://doi.org/10.5061/dryad.0p2ngf20h"
EXPECTED_ALPINE_FILE = "Alpine_Check_List.xlsx"
PERMITTED = {"annual":"annual","perennial":"perennial","biennial":"biennial"}


def normal(name: str) -> str:
    s = re.sub(r"\s+", " ", (name or "").strip()).casefold()
    return re.sub(r"[.]$", "", s)


def binomial(name: str) -> str | None:
    """Conservative species-level two-word match, never rank/unknown aggregates."""
    chunks = re.findall(r"[A-Za-z][A-Za-z-]*",name or "")
    if len(chunks) < 2:
        return None
    genus, epi = chunks[:2]
    if (not genus[:1].isupper() or not epi[:1].islower()
            or epi.casefold() in {"sp","spp","cf","aff","hybrid"}):
        return None
    return genus.casefold()+" "+epi.casefold()


def excel_col(addr: str) -> int:
    match = re.fullmatch(r"([A-Z]+)[0-9]+", addr or "")
    if not match:
        raise ValueError("Invalid Excel ref")
    n = 0
    for ch in match.group(1):
        n = 26*n + ord(ch)-64
    return n-1


def xlsx_rows(path: Path) -> list[list[object | None]]:
    """Original XLSX rows parsed using standard library; no external workbook rewrite."""
    data=path.read_bytes()
    if not data.startswith(b"PK\x03\x04"):
        raise ValueError("External life-span source is not XLSX")
    with zipfile.ZipFile(path) as z:
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            strings=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(t.text or "" for t in si.iter("{%s}t"%NS))
                    for si in strings.findall("{%s}si"%NS)]
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        mapping={el.get("Id"):el.get("Target") for el in rel.findall("{%s}Relationship"%PKG_NS)
                 if "worksheet" in (el.get("Type") or "")}
        out=[]
        for sheet in wb.findall("{%s}sheets/{%s}sheet"%(NS,NS)):
            target=mapping.get(sheet.get("{%s}id"%REL_NS))
            if not target:
                continue
            target="xl/"+target.lstrip("/").removeprefix("xl/")
            if target not in z.namelist():
                raise ValueError("Missing workbook worksheet")
            root=ET.fromstring(z.read(target))
            for row in root.findall("{%s}sheetData/{%s}row"%(NS,NS)):
                entries={}
                for cell in row.findall("{%s}c"%NS):
                    idx=excel_col(cell.get("r"))
                    typ=cell.get("t","")
                    v=cell.find("{%s}v"%NS)
                    value=None if v is None else v.text
                    if typ=="s":
                        value=shared[int(value)]
                    elif typ=="inlineStr":
                        value="".join(t.text or "" for t in cell.iter("{%s}t"%NS))
                    entries[idx]=value
                if entries:
                    out.append([entries.get(i) for i in range(max(entries)+1)])
    return out


def detect_life_span_table(rows: list[list[object | None]]) -> tuple[list[dict],dict]:
    """Use exact source headers, avoid assuming the first row is the header."""
    seen_headers=[]
    for i, r in enumerate(rows):
        cells=[normal(str(x or "")).replace("_"," ") for x in r]
        scientific=next((j for j,x in enumerate(cells) if x in
                         {"scientific name","scientific names","species name","taxon","plant species"}),None)
        life=next((j for j,x in enumerate(cells) if x in
                   {"life span","lifespan","life-span","life cycle","life history"}),None)
        syn=next((j for j,x in enumerate(cells) if x in {"common synonym","synonym","synonyms"}),None)
        if scientific is not None and life is not None:
            seen_headers.append((i,scientific,life,syn))
    if len(seen_headers)!=1:
        raise ValueError("Expected exactly one external sheet header containing scientific name/life span")
    i, scientific, life, syn=seen_headers[0]
    taxa=[]
    for row in rows[i+1:]:
        if len(row)<=max(scientific,life) or not row[scientific]:
            continue
        name=str(row[scientific]).strip()
        span=normal(str(row[life] or ""))
        if span not in PERMITTED:
            if span:
                raise ValueError("Non-standard external source span: "+span)
            continue
        synonym=str(row[syn] or "").strip() if syn is not None and len(row)>syn else ""
        taxa.append({"name":name,"binomial":binomial(name),"synonym":synonym,
                     "span":PERMITTED[span]})
    if not 100 <= len(taxa) <= 250:
        raise ValueError("Unexpected external flora row count, inspect source")
    return taxa,{"source_original_rows_with_known_lifespan":len(taxa),
                 "external_lifespan_counts":dict(collections.Counter(x["span"] for x in taxa)),
                 "original_source_header_row_number":i+1}


def crosswalk(original: list[str], external: list[dict], allow_binomial: bool) -> tuple[list[dict],dict]:
    """Avoid non-monophyletic/ambiguous genus-only or fuzzy-mapped matches."""
    text_to_row=collections.defaultdict(list)
    binom_to_rows=collections.defaultdict(list)
    for r in external:
        text_to_row[normal(r["name"])].append(r)
        if r["binomial"]:
            binom_to_rows[r["binomial"]].append(r)
    original_binomial_counts=collections.Counter(binomial(x) for x in original)
    out=[]
    rejected=[]
    match_counts=collections.Counter()
    for name in original:
        direct=text_to_row.get(normal(name), [])
        match_class="exact_source_name"
        if direct:
            candidates=direct
        elif allow_binomial:
            bn=binomial(name)
            candidates=binom_to_rows.get(bn, []) if bn and original_binomial_counts[bn]==1 else []
            match_class="unambiguous_species_binomial"
        else:
            candidates=[]
        if candidates and len({x["span"] for x in candidates})==1 and (
            len(candidates)==1 or len({normal(x["name"]) for x in candidates})==1):
            donor=candidates[0]
            out.append({"original_taxon":name,"life_span":donor["span"],
                        "trait_source":ALPINE_DOI+" | "+match_class+" | "+donor["name"]})
            match_counts[match_class]+=1
        else:
            rejected.append(name)
    return out, {"rail_original_taxon_count":len(original),
                 "classified_original_taxa":len(out),
                 "unmatched_original_taxa":len(rejected),
                 "match_tiers":dict(match_counts),
                 "unmatched_taxon_names_not_redistributed":True,
                 "trait_identification": "original 2025 Kashmir alpine flora; no genus-level inference",
                 "binomial_matching_uses_exact_species_epithet_not_taxonomic_synonym_resolution":True}


def coverage_summary(plots, classified):
    sites=sorted({k[0] for k in plots})
    sites=[s for s in sites if all(plots.get((s,y,z)) is not None
                        for y in (2014,2017) for z in
                        ("Track","Embankment","Interior Plot 1","Interior Plot 2"))]
    cells={}
    for zone in ("Track","Embankment","Interior Plot 1","Interior Plot 2"):
        fractions=[annual_stat(plots[s,y,zone],classified)["coverage"]
                   for s in sites for y in (2014,2017)]
        cells[zone]={"eligible_plot_years":len(fractions),
                     "median_known_fraction":statistics.median(fractions),
                     "min_known_fraction":min(fractions),
                     "max_known_fraction":max(fractions),
                     "plot_years_at_least_70pct":sum(f>=.7 for f in fractions),
                     "plot_years_at_least_50pct":sum(f>=.5 for f in fractions)}
    return cells


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rail",type=Path,required=True)
    ap.add_argument("--alpine",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--crosswalk-output",type=Path,required=True)
    args=ap.parse_args()
    original=original_taxon_headers(args.rail)
    plots,raw_meta=parse_xlsx_matrix(args.rail)
    taxa,alpine_meta=detect_life_span_table(xlsx_rows(args.alpine))
    matched, match_meta=crosswalk(original,taxa,allow_binomial=True)
    exact,_=crosswalk(original,taxa,allow_binomial=False)
    args.crosswalk_output.parent.mkdir(parents=True,exist_ok=True)
    with args.crosswalk_output.open("w",newline="",encoding="utf-8") as fp:
        wr=csv.DictWriter(fp,fieldnames=["original_taxon","life_span","trait_source"])
        wr.writeheader();wr.writerows(matched)
    original_classified,check=life_span_crosswalk(args.crosswalk_output,original)
    if len(original_classified)!=len(matched):
        raise AssertionError("Generated mapping failed v28.5 fail-closed check")
    both_year=coverage_summary(plots,original_classified)
    only_exact={original.index(x["original_taxon"]):x["life_span"] for x in exact}
    diagnostics={
       "study_name":"Rashid 2021 rail x Rasray and Lone 2025 independent alpine flora",
       "original_rail_data_sha256":DRYAD_SHA256,
       "external_alpine_workbook_sha256":hashlib.sha256(args.alpine.read_bytes()).hexdigest(),
       "external_reference":ALPINE_DOI,
       "rail_reference":RAIL_DOI,
       "independent_studies":2,
       "independent_railway_effect_studies":1,
       "external_flora_is_trait_lookup_not_railway_replication":True,
       "external_raw_flora":alpine_meta,
       "crosswalk_summary":match_meta,
       "exact_name_only_count":len(exact),
       "coverage":both_year,
       "strict_70pct":summarize(plots,original_classified,min_coverage=.7,bootstrap=999),
       "strict_50pct_exploratory":summarize(plots,original_classified,min_coverage=.5,bootstrap=999),
       "exact_name_only_70pct":summarize(plots,only_exact,min_coverage=.7,bootstrap=999),
       "source_taxon_data_not_redistributed":True
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(diagnostics,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"observed_rail_sites":raw_meta["all_spatial_sites"],
                      "external_taxa_with_life_span":len(taxa),
                      "rail_taxa_exact_and_unambiguous_binomial_matches":len(matched),
                      "rail_taxa_only_exact_fullname":len(exact),
                      "rail_plot_coverage":both_year,
                      "strict_track_interior1_main":
                       diagnostics["strict_70pct"]["comparisons"]["Track__minus__Interior Plot 1"]["all"],
                      "summary":str(args.output)},ensure_ascii=False))


if __name__=="__main__":
    main()
