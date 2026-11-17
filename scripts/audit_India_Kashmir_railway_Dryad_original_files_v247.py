#!/usr/bin/env python3
"""v24.7 authentic Dryad dataset discovery for non-European, replicated
railway vs adjacent vegetation controls — Kashmir Himalaya Rashid et al.
2021 J Applied Ecology DOI 10.1111/1365-2664.13961.

Study originals: 31 T-transects (16 real stations + 15 between-station),
4 genuine 2×50m zones per spatial site (track edge, embankment,
interior-1, interior-2), 2014 and 2017 resurveys, 123 vascular taxa,
~31×4×2 zone×year opportunities, original abundance recorded in
12 quadrats per site. Dryad data deposit DOI has two SPELLINGS seen in
publication/mirrors: 10.5061/dryad.0p2ngf20h vs
10.5061/dryad.0pngf20h — only an actual Dryad first-party API record
with matching dataset TITLE and contents can settle which is correct.

Discover ONLY real version/file metadata, URL, byte sizes, MD5 and title
from Dryad open REST. No file contents or fabricated occurrence matrix
if unavailable. No data from unverified PDF conversion or guessed link.
"""
from __future__ import annotations
import argparse,json,hashlib,urllib.parse,urllib.request,urllib.error,time
from pathlib import Path

DOI_VARIANTS=("10.5061/dryad.0p2ngf20h","10.5061/dryad.0pngf20h")
ROOT="https://datadryad.org/api/v2"
USER_AGENT="tetsu-original-railway-source-provenance-audit/1.0 (public ecological research)"
MAX_META=3_000_000

def getjson(url):
    req=urllib.request.Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
    last=None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=35) as r:
                data=r.read(MAX_META+1)
                if len(data)>MAX_META:raise ValueError("Source Dryad metadata unusually large")
                if r.status!=200:raise ValueError("Dryad metadata not successful HTTP")
            return json.loads(data)
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError) as e:
            last=str(e)
            if isinstance(e,urllib.error.HTTPError) and e.code in (400,404,410):break
            if attempt<2:time.sleep(2+2*attempt)
    raise RuntimeError(last)

def discover():
    checks=[];chosen=None
    for doi in DOI_VARIANTS:
        url=ROOT+"/datasets/"+urllib.parse.quote("doi:"+doi,safe="")
        try:
            d=getjson(url)
            title=str(d.get("title") or d.get("datasetTitle") or "")
            check={"doi_tested":doi,"api_url":url,"status":"ORIGINAL_DRYAD_RECORD_FOUND",
                   "returned_id":d.get("id"),"title":title,
                   "dataset_version_links":d.get("_links"),"declared_doi":d.get("identifier"),
                   "versionNumber":d.get("versionNumber"),"payload_keys":sorted(d)}
            if not all(k in title.lower() for k in ("railway","kashmir")):
                check["status"]="INVALID_TITLE_OR_DATASET_CANNOT_USE"
                checks.append(check);continue
            chosen=(doi,d);checks.append(check);break
        except Exception as e:
            checks.append({"doi_tested":doi,"api_url":url,
                           "status":"SOURCE_NOT_FETCHED","error":str(e)[:320]})
    if chosen is None:
        return {"version":"v24.7","status":"NO_FIRST_PARTY_DRYAD_DATASET_CONFIRMED","doi_variants":checks,
                "not_an_eligible_plant_site_matrix":True}
    doi,d=chosen
    dataset_id=d.get("id")
    verlist=None
    for target in (ROOT+"/datasets/"+urllib.parse.quote("doi:"+doi,safe="")+"/versions",
                   ROOT+"/datasets/"+str(dataset_id)+"/versions"):
        try:
            verlist=getjson(target);break
        except Exception as e:
            checks.append({"versions_api_attempt":target,"error":str(e)[:240]})
    # Dryad REST versions response is paginated HAL, with _embedded
    version_rows=[]
    if isinstance(verlist,dict):
        embedded=verlist.get("_embedded") or {}
        version_rows=(embedded.get("stash:versions") or embedded.get("versions") or
                      verlist.get("versions") or [])
    if not version_rows and d.get("versionNumber"):
        version_rows=[{"versionNumber":d["versionNumber"],"id":d.get("stashVersionId")}]
    files=[]
    for version in version_rows[-4:]:
        versionid=version.get("id") or version.get("stashVersionId")
        if versionid is None:continue
        url=ROOT+"/versions/"+str(versionid)+"/files"
        try:
            listing=getjson(url)
        except Exception as e:
            checks.append({"file_list_api_attempt":url,"error":str(e)[:300]});continue
        embedded=listing.get("_embedded") or {}
        listings=embedded.get("stash:files") or embedded.get("files") or listing.get("files") or []
        for file in listings:
            files.append({"version_id":versionid,
                "file_id":file.get("id"),"path":file.get("path") or file.get("filename"),
                "size":file.get("size") or file.get("length"),
                "digest":file.get("digest") or file.get("checksum"),
                "mime":file.get("mimeType") or file.get("type"),
                "links":file.get("_links")})
    return {"version":"v24.7","status":"FIRST_PARTY_DRYAD_RAILWAY_KASHMIR_DATASET_CONFIRMED",
       "verified_original_deposit_DOI":doi,"verified_title":str(d.get("title")),
       "author_first_party_original_dataset_id":dataset_id,
       "original_version_count_returned":len(version_rows),
       "original_data_files_found":files,"files_list_count":len(files),
       "source_discovery":checks,
       "cannot_identify_rail_vs_nonrail_plots_until_original_excel_rows_checked":True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",required=True)
    a=parser.parse_args();result=discover()
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print("V247_INDIA_KASHMIR_DRYAD_DISCOVERY="+json.dumps(result,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
