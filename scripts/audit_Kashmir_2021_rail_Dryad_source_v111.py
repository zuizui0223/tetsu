#!/usr/bin/env python3
"""v11.1: read-only audit of truly published 2021 Kashmir Himalayan RAIL Dryad source.

Original paper Rashid et al. (2021), J Appl Ecol DOI 10.1111/1365-2664.13961,
31 T transects 2014 and 2017. The original paper's bibliography cites a
possibly malformed data DOI ending 0pngf20h. Independent US NLM dataset
catalog records 10.5061/dryad.0p2ngf20h (with digit '2'), CC0.
Only trust live publisher API metadata after its DOI/title validate.
DO NOT invent linked data or infer true freight transport from plot records.
"""
from __future__ import annotations
import argparse,json,urllib.request,urllib.parse,urllib.error
from pathlib import Path
from datetime import datetime,timezone

PAPER_DOI="10.1111/1365-2664.13961"
NLM_CATALOG_DATA_DOI="10.5061/dryad.0p2ngf20h"
PAPER_REFERENCE_TYPO_CANDIDATE="10.5061/dryad.0pngf20h"
BASE="https://datadryad.org/api/v2/"


def url_for_dataset(doi):
    return BASE+"datasets/"+urllib.parse.quote("doi:"+doi,safe="")


def get_json(url):
    req=urllib.request.Request(url,headers={
        "Accept":"application/json",
        "User-Agent":"tetsu-academic-open-source-audit/1.0 (public metadata only)"})
    with urllib.request.urlopen(req,timeout=40) as r:
        raw=r.read(2500000)
        if len(raw)>2000000:
            raise ValueError("Publisher API metadata unexpectedly oversized")
        obj=json.loads(raw)
        if not isinstance(obj,dict):
            raise ValueError("Original Dryad dataset metadata not JSON object")
        return obj


def extract_links(obj,substring):
    links=obj.get("_links",{})
    if not isinstance(links,dict):return []
    urls=[]
    for key,value in links.items():
        if substring.lower() not in key.lower():continue
        vs=value if isinstance(value,list) else [value]
        for v in vs:
            if isinstance(v,dict) and v.get("href"):
                u=urllib.parse.urljoin("https://datadryad.org",v["href"])
                if u.startswith("https://datadryad.org/api/v2/"):
                    urls.append(u)
    return list(dict.fromkeys(urls))


def summarize_file_list(response):
    block=response.get("_embedded",{})
    arr=[]
    if isinstance(block,dict):
        for key,val in block.items():
            if "file" in key.lower() and isinstance(val,list):
                arr.extend(val)
    if isinstance(response.get("files"),list):arr.extend(response["files"])
    return [{"id":r.get("id"),
             "path":r.get("path") or r.get("filename") or r.get("name"),
             "size":r.get("size"),
             "mimeType":r.get("mimeType"),
             "digest":r.get("digest"),
             "digestType":r.get("digestType"),
             "download_urls":extract_links(r,"download")} for r in arr if isinstance(r,dict)]


def audit(source_get=get_json):
    result={
        "complete":False,"version":"v11.1",
        "original_railway_field_paper_DOI":PAPER_DOI,
        "original_NLM_catalog_Dryad_DOI":NLM_CATALOG_DATA_DOI,
        "literature_reference_malformed_candidate":PAPER_REFERENCE_TYPO_CANDIDATE,
        "catalog_source":"https://datasetcatalog.nlm.nih.gov/dataset?q=0000002490",
        "publisher_api_dataset_url":url_for_dataset(NLM_CATALOG_DATA_DOI),
        "source_data_raw_downloaded_or_ecology_fitted":False,
        "record_UTC":datetime.now(timezone.utc).isoformat(),
    }
    try:
        meta=source_get(result["publisher_api_dataset_url"])
    except Exception as ex:
        result["failure_kind"]=type(ex).__name__
        result["failure_details"]=str(ex)[:500]
        result["status"]="original_Dryad_metadata_not_verified_no_biological_result"
        return result
    title=str(meta.get("title") or "")
    doi=str(meta.get("identifier") or meta.get("doi") or "")
    titleok=("kashmir" in title.lower() and
             ("rail" in title.lower() or "species abundance" in title.lower()))
    doiok=(NLM_CATALOG_DATA_DOI.lower() in doi.lower() or
           NLM_CATALOG_DATA_DOI.lower() in
           str(meta.get("identifier","")).lower())
    result["publisher_dataset_metadata"]={
        "id":meta.get("id"),"title":title,"identifier":doi,
        "versionNumber":meta.get("versionNumber"),
        "latestVersion":meta.get("latestVersion"),
        "versionStatus":meta.get("versionStatus"),
        "doi_metadata_match":doiok,"railway_title_match":titleok,
        "available_rel_names":list(meta.get("_links",{})) if isinstance(meta.get("_links"),dict) else [],
    }
    if not titleok:
        result["status"]="Dryad_source_title_mismatch_no_biological_result"
        return result
    if not doiok:
        result["status"]="Dryad_DOI_metadata_unverified_no_biological_result"
        return result
    result["metadata_identity_confirmed"]=True
    versionurl=list(dict.fromkeys(extract_links(meta,"version")+
        [result["publisher_api_dataset_url"]+"/versions"]))
    errors=[]
    versions=[]
    for u in versionurl[:4]:
        try:
            j=source_get(u)
        except Exception as e:
            errors.append(type(e).__name__+": "+str(e)[:120])
            continue
        embed=j.get("_embedded",{})
        for key,vals in (embed.items() if isinstance(embed,dict) else []):
            if "version" in key.lower() and isinstance(vals,list):
                versions.extend(vals)
        if isinstance(j.get("versions"),list):versions.extend(j["versions"])
        if j.get("id") is not None and j.get("versionNumber") is not None:
            versions.append(j)
    summaries=[]
    for v in versions:
        if not isinstance(v,dict):continue
        id=v.get("id")
        if id is None:continue
        summaries.append({"id":id,"versionNumber":v.get("versionNumber"),"versionStatus":v.get("versionStatus")})
    result["publisher_dataset_versions"]=summaries
    candidateurls=extract_links(meta,"files")
    for v in summaries:
        candidateurls.append(BASE+"versions/"+str(v["id"])+"/files")
    candidates=list(dict.fromkeys(candidateurls))
    files=[]
    for u in candidates:
        try:
            raw=source_get(u)
            files.extend(summarize_file_list(raw))
        except Exception as e:
            errors.append(type(e).__name__+": "+str(e)[:120])
    seen=set()
    unique=[]
    for f in files:
        key=(str(f["id"]),str(f["path"]))
        if key in seen:continue
        seen.add(key);unique.append(f)
    result["file_list_source_urls"]=candidates
    result["original_publisher_dataset_files"]=unique
    result["n_original_metadata_files"]=len(unique)
    result["publisher_API_retrieval_warnings"]=errors[:8]
    result["complete"]=len(unique)>0
    result["status"]="original_Dryad_file_inventory_verified" if unique else "original_Dryad_metadata_only_files_not_verified"
    result["original_DOI_is_not_a_railway_wind_or_mowing_causal_test"]=True
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=audit()
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)


if __name__=="__main__":main()
