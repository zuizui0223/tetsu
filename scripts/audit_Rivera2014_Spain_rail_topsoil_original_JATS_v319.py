#!/usr/bin/env python3
"""First-party source JATS verification for Spain 2014 railroad embankment RCT.

No raw plot data are inferred from publication means±SE. Read DOI-locked
Europe PMC copy of original PLOS CC-BY JATS XML, output only tables needed
for explicit original author-reported treatment time comparisons.
"""
from __future__ import annotations
import argparse,hashlib,json,re,xml.etree.ElementTree as ET
from pathlib import Path

SOURCE_DOI="10.1371/journal.pone.0101413"
SOURCE_JATS_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4077791/fullTextXML"
SOURCE_ORIGINAL="https://journals.plos.org/plosone/article?id="+SOURCE_DOI

def clean(x):
    return " ".join("".join(x.itertext()).split())

def original_xml_tables(path):
    raw=path.read_bytes()
    if not (raw.lstrip().startswith(b"<?xml") or b"<article" in raw[:350]) or len(raw)<25000:
        raise ValueError("Not genuine original JATS full text")
    root=ET.fromstring(raw)
    dois={clean(x) for x in root.findall(".//article-id[@pub-id-type='doi']")}
    pmcid={clean(x) for x in root.findall(".//article-id[@pub-id-type='pmc']")}
    if SOURCE_DOI not in dois or (pmcid and "PMC4077791" not in pmcid and "4077791" not in pmcid):
        raise ValueError("Original Spain original source DOI and PMCID do not match")
    content=clean(root)
    required=["railway","topsoil","control","2010","2011","random"]
    for marker in required:
        if marker not in content.lower():
            raise ValueError("Original rail topsoil RCT missing locked source term "+marker)
    original_tables=[]
    for t in root.findall(".//table-wrap"):
        label=clean(t.find("label")) if t.find("label") is not None else t.get("id","")
        rows=[]
        for tr in t.findall(".//tr"):
            vals=[clean(c) for c in tr if c.tag.rsplit("}",1)[-1] in ("th","td")]
            if vals:rows.append(vals)
        caption=t.find("caption")
        original_tables.append({"label":label,
                               "source_table_id":t.get("id"),
                               "source_author_caption":clean(caption)[:500] if caption is not None else None,
                               "rows":rows[:55],"n_source_extracted_rows":len(rows),
                               "source_table_is_embedded_in_primary_manuscript":True})
    if not original_tables:
        raise ValueError("Original PLOS full text contains no figure/table source XML")
    context=[]
    for title in ("Study area","Experimental design","Topsoil spreading and total plant cover"):
        for sec in root.findall(".//sec"):
            node=sec.find("title")
            if node is not None and title.lower() in clean(node).lower():
                context.append({"section":title,"source_excerpt":clean(sec)[:1200]})
                break
    return {
       "version":"v319_PLOS_original_JATS_source",
       "author_paper_DOI":SOURCE_DOI,"author_original_publication":SOURCE_ORIGINAL,
       "original_open_access_full_text":SOURCE_JATS_URL,
       "true_original_article_sha256":hashlib.sha256(raw).hexdigest(),
       "original_full_text_bytes":len(raw),
       "original_source_doi_verified":True,
       "original_source_tables":original_tables,
       "original_source_context":context,
       "source_meta_review_country_error":"Mell2024 S5: Waterway and USA, but original Spain railway",
       "randomized_field_sample_units":"5 original source blocks × 3 treatments, each 6x9m 54m², repeated in 2010 and 2011",
       "treatment_control_C":"original railway embankment substrate",
       "treatment_topsoil_10":"topsoil 10cm",
       "treatment_topsoil_30":"topsoil 30cm",
       "n_new_railway_vs_nonrail_studies":0,
       "original_individual_15_plot_values_verified":False,
       "no_fake_new_DID_precision_or_confidence_interval":True,
       "status":"SOURCE_VERIFIED_TABLES_NOT_AUTOMATIC_META_OR_CAUSAL_MEDIATION"
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--xml",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    a=original_xml_tables(args.xml)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(a,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
      "original_DOI":a["author_paper_DOI"],
      "original_JATS_SHA256":a["true_original_article_sha256"],
      "source_tables":[{"label":s["label"],"caption":s["source_author_caption"],
                        "n_rows":s["n_source_extracted_rows"],"first_six_rows":s["rows"][:6]}
                      for s in a["original_source_tables"]],
      "source_stage":a["status"],
    },indent=2,ensure_ascii=False))
if __name__=="__main__":main()
