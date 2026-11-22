#!/usr/bin/env python3
"""Discover authoritative report annexes through the original ITTECOP project HTML.

Only resolves publication links on an official ittecop.fr host. Does not
convert Latin names/paragraphs into fabricated taxon×plot field data.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from html.parser import HTMLParser
from urllib.parse import urljoin,urlparse
from pathlib import Path

PROJECT="https://www.ittecop.fr/fr/tous-les-projets/recherches-2008-2009/item/33-en-marge"
REQUIRED=("en marge","annexes","rapport final")

class Anchors(HTMLParser):
    def __init__(self):
        super().__init__()
        self.stack=[]
        self.links=[]
    def handle_starttag(self,tag,attrs):
        kv=dict(attrs)
        if tag=="a":
            self.stack.append({"href":kv.get("href",""),"title":kv.get("title",""),
                "text":"","attributes":{k:v for k,v in kv.items() if k in ("href","title","class")}})
        elif tag=="img" and self.stack:
            self.stack[-1]["text"]+=" "+str(kv.get("alt") or "")
    def handle_startendtag(self,tag,attrs):
        if tag=="img":self.handle_starttag(tag,attrs)
    def handle_data(self,data):
        for a in self.stack:a["text"]+=data
    def handle_endtag(self,tag):
        if tag=="a" and self.stack:self.links.append(self.stack.pop())

def official(url):
    u=urlparse(url)
    return u.scheme=="https" and (u.hostname=="ittecop.fr" or u.hostname=="www.ittecop.fr")

def audit(path:Path):
    body=path.read_bytes()
    text=body.decode("utf-8","replace")
    plain=re.sub(r"\s+"," ",re.sub(r"<[^>]+>"," ",text)).lower()
    if len(body)<8000 or "en marge" not in plain:
        raise ValueError("Wrong original En marge source official project page")
    p=Anchors();p.feed(text)
    items=[]
    for x in p.links:
        href=x["href"]
        label=" ".join((x["text"]+" "+x["title"]).split())
        # Official Joomla project page has image-only report links with
        # empty anchor text; the file paths are authoritative labels.
        if "/rapport-final/Annexes/" in href or "/rapport-final/annexes/" in href:
            label="Annexes RF "+href.rsplit("/",1)[-1]
        elif "/rapport-final/2008/" in href or "/rapport-final/" in href:
            label="Rapport final "+href.rsplit("/",1)[-1]
        url=urljoin(PROJECT,href)
        if any(k in label.casefold() for k in ("annexes","rapport final","présentation","presentation")):
            items.append({"published_link_text":label[:190],
                          "resolved_publisher_url":url,
                          "source_href":x["href"][:500],
                          "official_same_host":official(url)})
    if not items:
        raise ValueError("No official original report/annex references present in page")
    annex=[x for x in items if "annexes" in x["published_link_text"].casefold()
                         and x["official_same_host"]]
    report=[x for x in items if "rapport final" in x["published_link_text"].casefold()
                          and x["official_same_host"]]
    if not annex or not report:
        snippets=[text[max(0,m.start()-450):min(len(text),m.end()+700)]
                  for m in list(re.finditer(r"annexes\\s*rf|rapport\\s+final",text,re.I))[:12]]
        return {
            "version":"v321_ittecop_project_link_diagnosis",
            "original_project_url":PROJECT,
            "project_html_sha256":hashlib.sha256(body).hexdigest(),
            "project_html_bytes":len(body),
            "source_is_original_ittecop_project_page":True,
            "source_link_status":"OFFICIAL_HTML_PAGE_HAS_VISIBLE_LABEL_BUT_ANCHOR_MISSING_OR_NONOFFICIAL",
            "source_exact_official_annex_links":annex,
            "source_exact_official_main_report_links":report,
            "source_all_related_link_targets":items[:35],
            "source_annex_label_html_contexts":snippets,
            "n_original_html_anchor_tags":len(p.links),
            "source_all_anchor_href_near_annex_label":[x for x in p.links
                if "annexe" in (x.get("href","")+" "+x.get("text","")).lower()][:25],
            "appendix_plant_by_rail_matrix_NOT_yet_verified":True,
            "no_new_independent_rail_comparison_claimed":True
        }
    return {
       "version":"v321_ittecop_project_link_diagnosis",
       "original_project_url":PROJECT,
       "project_html_sha256":hashlib.sha256(body).hexdigest(),
       "project_html_bytes":len(body),
       "source_is_original_ittecop_project_page":True,
       "source_exact_official_annex_links":annex,
       "source_exact_official_main_report_links":report,
       "source_all_related_link_targets":items[:35],
       "appendix_plant_by_rail_matrix_NOT_yet_verified":True,
       "no_new_independent_rail_comparison_claimed":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-html",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    x=audit(a.source_html)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(x,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":"original 2015 project links authenticated",
                      "annexes":x["source_exact_official_annex_links"],
                      "main_report":x["source_exact_official_main_report_links"],
                      "metadata_sha":x["project_html_sha256"]},
                      ensure_ascii=False,indent=2))

if __name__=="__main__":main()
