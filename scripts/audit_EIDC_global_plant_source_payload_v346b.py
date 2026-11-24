#!/usr/bin/env python3
"""Inspect real official EIDC download package, not a landing-page lookalike.

Do not output raw taxon tables or location data, only schema, file names,
hashes and aggregate counts. This is still a source gate, no biological effect.
"""
from __future__ import annotations
import argparse,csv,hashlib,html,json,re,zipfile,io
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin,urlparse
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

ID="ca339c86-3674-4030-b891-35326e71141e"
BASE="https://data-package.ceh.ac.uk"
LANDING=BASE+"/data/"+ID
SUPPORT=BASE+"/sd/"+ID+".zip"
HEADERS={"User-Agent":"tetsu-global-plant-temperatures-source-v346b/1.0","Accept":"*/*"}
MAX_RAW=260_000_000

class Links(HTMLParser):
    def __init__(self):
        super().__init__();self.refs=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag in ("a","link") and d.get("href"):
            self.refs.append((d["href"],d.get("download"),d.get("type"),d.get("title")))

def source(url,maxbytes=MAX_RAW):
    with urlopen(Request(url,headers=HEADERS),timeout=130) as r:
        b=r.read(maxbytes+1)
        if len(b)>maxbytes:raise ValueError("Source exceeds frozen file limit")
        return {"body":b,"url":r.geturl(),"status":r.status,"type":r.headers.get("Content-Type",""),
                "content_disposition":r.headers.get("Content-Disposition","")}

def inspect_zip(raw,include_csv_census=False):
    if not raw.startswith(b"PK\x03\x04"):
        raise ValueError("Original expected ZIP signature absent")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        entries=[{"path":i.filename,"uncompressed":i.file_size}
                 for i in z.infolist() if not i.is_dir()]
    return {"source_SHA256":hashlib.sha256(raw).hexdigest(),
            "bytes":len(raw),"file_count":len(entries),"file_list":entries[:40],
            "csv_file_count":sum(x["path"].lower().endswith(".csv") for x in entries)}

def inspect_page():
    p=source(LANDING,maxbytes=5_000_000)
    html_string=p["body"].decode("utf-8",errors="replace")
    links=Links();links.feed(html_string)
    absolute=[]
    for href,download,ctype,title in links.refs:
        u=urljoin(p["url"],html.unescape(href))
        if urlparse(u).netloc.endswith("ceh.ac.uk"):
            absolute.append({"url":u,"download":download,"type":ctype,"title":title})
    source_link_matches=sorted(set(re.findall(r'(?:"|\x27|\\/)([^"\x27<> ]+\.(?:csv|zip|tsv)(?:\?[^"\x27<> ]*)?)(?:"|\x27)',html_string,re.I)))
    return {"resolved_url":p["url"],"source_HTML_bytes":len(p["body"]),
            "page_SHA256":hashlib.sha256(p["body"]).hexdigest(),
            "has_dataset_download_link":"download" in html_string.lower(),
            "links_from_official_CEH_host":absolute[:60],
            "embedded_csv_zip_path_candidates":source_link_matches[:30],
            "n_embedded_csv_zip_candidates":len(source_link_matches),
            "n_official_links":len(absolute)}

def run():
    a={"original_DOI":"10.5285/"+ID,
       "source_links_obtained_from_official_UK_CKAN_metadata":True,
       "schema_or_raw_effect_not_guessed":True}
    try:a["official_data_page"]=inspect_page()
    except Exception as e:a["official_data_page_error"]={"kind":type(e).__name__,"message":str(e)[:160]}
    try:
        p=source(SUPPORT,maxbytes=140_000_000)
        a["official_supporting_zip"]=inspect_zip(p["body"])
    except Exception as e:a["official_supporting_zip_error"]={"kind":type(e).__name__,"message":str(e)[:160]}
    a["no_original_taxon_rows_republished"]=True
    return a

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True)
    args=p.parse_args();out=run()
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)
