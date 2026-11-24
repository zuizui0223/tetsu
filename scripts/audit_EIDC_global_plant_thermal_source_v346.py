#!/usr/bin/env python3
"""Check authoritative EIDC species global temperature envelope data access.

This does not infer any plant climate metric or railway ecology.
Do not accept publisher HTML landing pages as original CSV. No source raw
records exported. Metadata and source SHA only in aggregate report.
"""
from __future__ import annotations
import argparse,hashlib,json,os,re,time
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlparse
from urllib.request import Request,urlopen

DOI="10.5285/ca339c86-3674-4030-b891-35326e71141e"
IDENT="ca339c86-3674-4030-b891-35326e71141e"
REGISTRY="https://ckan.publishing.service.gov.uk/api/3/action/package_show?id=bioclimatic-envelopes-for-global-plant-species"
EIDC="https://data-package.ceh.ac.uk/data/"+IDENT
CATALOGUE="https://catalogue.ceh.ac.uk/documents/"+IDENT
MAX_HEADER=32768
AGENT="tetsu-plant-climate-original-source-feasibility/1.0 academic-study"

def probe(url,timeout=32):
    req=Request(url,headers={"User-Agent":AGENT,"Range":"bytes=0-32767","Accept":"application/json,text/csv,application/zip,*/*"})
    try:
        with urlopen(req,timeout=timeout) as response:
            b=response.read(MAX_HEADER)
            return {"url":url,"http":response.status,"redirected_url":response.geturl(),
                    "content_type":response.headers.get("Content-Type"),
                    "content_length_header":response.headers.get("Content-Length"),
                    "content_disposition":response.headers.get("Content-Disposition"),
                    "first_bytes_hex":b[:16].hex(),
                    "header":b,
                    "is_zip":b.startswith(b"PK\\x03\\x04") or b.startswith(b"PK\x03\x04"),
                    "is_html":b[:512].lstrip().lower().startswith((b"<!doctype html",b"<html")),
                   }
    except (HTTPError,URLError,TimeoutError,ValueError,OSError) as e:
        return {"url":url,"error_type":type(e).__name__,"http":getattr(e,"code",None),
                "short_error":str(e)[:130]}

def discover_resources():
    response=probe(REGISTRY)
    res={"ckan_probe":{k:v for k,v in response.items() if k!="header"},
         "official_candidate_resource_links":[]}
    if "header" in response:
        try:
            payload=json.loads(response["header"])
            if payload.get("success") and isinstance(payload.get("result"),dict):
                pkg=payload["result"]
                res["dataset_title"]=pkg.get("title")
                res["resource_list_size"]=len(pkg.get("resources",[]))
                res["official_candidate_resource_links"]=[
                    {"id":r.get("id"),"name":r.get("name"),"format":r.get("format"),
                     "size":r.get("size"),"url":r.get("url")}
                    for r in pkg.get("resources",[]) if r.get("url")]
        except (ValueError,KeyError,TypeError):
            res["ckan_truncated_or_not_json"]=True
    urls=[r["url"] for r in res["official_candidate_resource_links"]]
    urls.extend([EIDC,CATALOGUE+"?type=json"])
    verified=[]
    for url in dict.fromkeys(urls):
        p=probe(url)
        verified.append({k:v for k,v in p.items() if k!="header"})
    res["candidate_HTTP_probe"]=verified
    res["source_file_verified"]=any(x.get("is_zip") for x in verified)
    res["any_accessible_candidate"]=any(x.get("http") in (200,206) for x in verified)
    res["data_not_yet_downloaded_and_headers_not_validated"]=True
    return res

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--output",required=True)
    args=a.parse_args();r={"version":"v346-EIDC-global-plant-temperature-source-access",
        "source_DOI":DOI,"source_URL_from_metadata_only":True,
        "source_raw_plants_not_republished":True,"result":discover_resources()}
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(r,ensure_ascii=False),flush=True)
