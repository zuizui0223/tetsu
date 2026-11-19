#!/usr/bin/env python3
"""v28.9 ORIGINAL open Figshare audit: do non-rail corridor controls hold real plant
site-by-taxon records? This is NOT a railway analysis or a railway meta-effect.

Retrieve two publisher-cited DOI article IDs, check provider MD5 for each
public file, inspect only matrix SCHEMA and broad dimensions. Raw plot/taxon
rows and field coordinates MUST NOT be written into results/ or git.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import time
import urllib.request
import urllib.error
import zipfile
from pathlib import Path
from urllib.parse import urlparse

ARTICLES={
    13286141: "10.6084/m9.figshare.13286141",
    22178651: "10.6084/m9.figshare.22178651",
}
API="https://api.figshare.com/v2/articles/"
DOWNLOAD="https://ndownloader.figshare.com/files/"
MAX_FILE_BYTES=45_000_000
MAX_UNCOMPRESSED_BYTES=100_000_000
MAX_ZIP_MEMBERS=250
MAX_TABLE_COLUMNS=9000
REQUEST_TIMEOUT=90
USER_AGENT="tetsu-source-provenance-research/28.9 (figshare public metadata; no private locations published)"

def get(url: str, *, download=False) -> bytes:
    if download:
        assert url.startswith(DOWNLOAD)
    else:
        assert url.startswith(API)
    err=None
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":USER_AGENT}),
                                        timeout=REQUEST_TIMEOUT) as resp:
                if resp.status!=200:raise ValueError(f"Publisher HTTP {resp.status}")
                announced=resp.headers.get("Content-Length")
                if announced and int(announced)>MAX_FILE_BYTES and download:
                    raise ValueError("Original source file exceeds bounded audit size")
                out=resp.read(MAX_FILE_BYTES+1)
                if download and len(out)>MAX_FILE_BYTES:
                    raise ValueError("Original file too large for bounded QA")
                return out
        except (urllib.error.URLError,TimeoutError) as exc:
            err=exc
            if k<2:time.sleep(k+1)
    raise RuntimeError("Original Figshare API/download failed") from err

def plant_columns(header: list[str])->dict:
    """Only inspect column labels; never return raw source plant/plot rows."""
    # Scientific names as wide columns should not be accidentally detected
    # when one word is an author, location or survey code.
    latin=re.compile(r"^[A-Z][a-zA-Z-]{2,}\s+[a-z][a-z-]{2,}(?:\s|$)")
    name_hits=sum(bool(latin.match(str(v).strip())) for v in header)
    h=[re.sub(r"[^a-z0-9]+","",str(v).lower()) for v in header]
    return {
        "n_columns":len(header),
        "n_scientific_binomial_like_header_columns":name_hits,
        "taxon_long_column_present":any(x in {"species","speciesname","taxon","taxonname","scientificname","plantname"}
                                         for x in h),
        "plant_occurrence_column_present":any(x in {"abundance","cover","presence","occurrence",
                                "percentcover","speciescover","count","frequency"} for x in h),
        "site_identifier_column_present":any(x in {"plot","plotid","site","siteid","transect","transectid",
                                  "landscape","landscapeid","habitat","habitatcode","habitatid","quadrat","quadratid"}
                                           for x in h),
        "plausible_wide_site_by_taxon":name_hits>=15,
        "plausible_long_site_by_taxon":any(x in {"species","speciesname","taxon","taxonname","scientificname","plantname"}
                                           for x in h) and any(x in {"plot","plotid","site","siteid","transect",
                                                 "transectid","habitat","quadrat","quadratid"} for x in h),
    }

def summarize_table(content:bytes, suffix:str) -> dict:
    if suffix in (".csv",".tsv",".txt"):
        enc=None;text=None
        for candidate in ("utf-8-sig","utf-16","cp1252"):
            try:
                txt=content.decode(candidate)
                if txt.count("\x00") > len(txt)*.01:continue
                text=txt;enc=candidate;break
            except UnicodeError:pass
        if text is None:return {"table_status":"unreadable_encoding"}
        lines=text.splitlines()
        if not lines:return {"table_status":"empty"}
        try:
            sample="\n".join(lines[:min(30,len(lines))])
            dialect=csv.Sniffer().sniff(sample,delimiters=",;\t")
            delim=dialect.delimiter
        except csv.Error:
            delim="\t" if suffix==".tsv" else ","
        reader=csv.reader(lines,delimiter=delim)
        header=next(reader,[])
        # Whole-file count is metadata only (not original values).
        n_data_rows=sum(1 for r in reader if r)
        return {"table_status":"parsed","encoding":enc,"n_rows_excluding_first_line":n_data_rows,
                "delimiter":"TAB" if delim=="\t" else delim,
                "header_structure":plant_columns(header[:MAX_TABLE_COLUMNS]),
                "n_header_columns_truncated":max(0,len(header)-MAX_TABLE_COLUMNS)}
    if suffix==".xlsx":
        try:
            import openpyxl
            with io.BytesIO(content) as f:
                wb=openpyxl.load_workbook(f,read_only=True,data_only=True)
                out=[]
                for ws in wb.worksheets[:30]:
                    hdr=next(ws.values,[])
                    out.append({"sheet_name":ws.title[:100],
                        "source_max_row":ws.max_row,
                        "source_max_column":ws.max_column,
                        "header_structure":plant_columns([str(x) if x is not None else ""
                                                         for x in hdr[:MAX_TABLE_COLUMNS]])})
                wb.close()
                return {"table_status":"parsed","sheets":out,"n_sheets":len(wb.worksheets)}
        except Exception as exc:
            return {"table_status":"xlsx_reader_error","error_type":type(exc).__name__}
    return {"table_status":"not_tabular_ext"}

def inspect_file(content:bytes, name:str)->dict:
    ext=Path(name.lower()).suffix
    if ext==".zip":
        try:
            with zipfile.ZipFile(io.BytesIO(content)) as z:
                infos=[x for x in z.infolist() if not x.is_dir()]
                if len(infos)>MAX_ZIP_MEMBERS or sum(x.file_size for x in infos)>MAX_UNCOMPRESSED_BYTES:
                    return {"container":"zip","status":"bounded_source_exceeds_limits",
                            "member_count":len(infos)}
                return {"container":"zip","member_count":len(infos),
                        "members":[{
                            "file_suffix":Path(x.filename.lower()).suffix,
                            "member_size":x.file_size,
                            "table":summarize_table(z.read(x),Path(x.filename.lower()).suffix)
                            if Path(x.filename.lower()).suffix in (".csv",".tsv",".txt",".xlsx") and x.file_size<=MAX_FILE_BYTES
                            else {"table_status":"unscanned_non_table_or_too_large"}
                         } for x in infos[:MAX_ZIP_MEMBERS]]}
        except (zipfile.BadZipFile,RuntimeError,ValueError):
            return {"container":"zip","status":"invalid_archive"}
    return {"container":"single_file","table":summarize_table(content,ext)}

def article_audit(article_id:int)->dict:
    if article_id not in ARTICLES:raise ValueError("Study ID not predeclared")
    meta=json.loads(get(API+str(article_id)))
    if int(meta.get("id",-1))!=article_id:
        raise ValueError("Figshare article ID changed")
    # Figshare often stores a versioned DOI, e.g. .v1. Validate stable prefix.
    doi=str(meta.get("doi") or meta.get("resource_doi") or "")
    if ARTICLES[article_id] not in doi.lower():
        raise ValueError(f"Original Figshare DOI no longer matches frozen source: {doi}")
    files=meta.get("files")
    if not isinstance(files,list) or not files:
        raise ValueError("Publisher dataset has no openly enumerated source files")
    result={"publisher":"figshare","article_id":article_id,"doi":doi,
            "original_title":meta.get("title"),"original_article_version":meta.get("version"),
            "citation":str(meta.get("citation") or "")[:300],"source_file_count":len(files),
            "interpretation":"original DOI corresponds to NON-RAIL road/powerline/grassland corridors; never count as rail study",
            "source_files":[]}
    for f in files:
        fid=int(f["id"]);size=int(f.get("size") or 0)
        info={"file_id":fid,"file_name":str(f.get("name") or "")[:160],
              "original_publisher_bytes":size,
              "source_publisher_md5":str(f.get("supplied_md5") or f.get("computed_md5") or ""),
              "source_file_sha256":None,"source_checksum_verified":False}
        if size>MAX_FILE_BYTES:
            info["source_status"]="EXCLUDED_SIZE";result["source_files"].append(info);continue
        raw=get(DOWNLOAD+str(fid),download=True)
        if size!=len(raw):raise ValueError("Publisher file size mismatch")
        h=hashlib.md5(raw).hexdigest()
        if info["source_publisher_md5"] and h.lower()!=info["source_publisher_md5"].lower():
            raise ValueError("Original Figshare supplied MD5 mismatch")
        info["source_checksum_verified"]=bool(info["source_publisher_md5"])
        info["source_file_sha256"]=hashlib.sha256(raw).hexdigest()
        info["source_status"]="DOWNLOADED_AND_SCHEMA_AUDITED"
        info["schema_without_raw_rows"]=inspect_file(raw,str(f.get("name") or ""))
        result["source_files"].append(info)
    result["every_included_original_source_file_downloaded"]=all(
        f["source_status"]=="DOWNLOADED_AND_SCHEMA_AUDITED" for f in result["source_files"])
    result["n_potential_plant_site_taxon_files"]=sum(
        bool(any((sh.get("header_structure") or {}).get("plausible_wide_site_by_taxon") or
                 (sh.get("header_structure") or {}).get("plausible_long_site_by_taxon")
                 for sh in (f.get("schema_without_raw_rows",{}).get("table",{}).get("sheets") or [])))
        or bool((f.get("schema_without_raw_rows",{}).get("table",{}).get("header_structure") or {})
            .get("plausible_wide_site_by_taxon"))
        or bool((f.get("schema_without_raw_rows",{}).get("table",{}).get("header_structure") or {})
            .get("plausible_long_site_by_taxon"))
        for f in result["source_files"])
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    results=[article_audit(i) for i in ARTICLES]
    out={"version":"v28.9",
         "question":"Open independent nonrail transport corridor source botanical matrix support",
         "articles":results,"n_orig_author_public_datasets":len(results),
         "not_a_railway_field_replication":True,
         "no_field_effect_size_claimed":True,
         "no_raw_original_plant_site_taxon_or_private_coordinates_exported":True,
         "source_GitHub_public_artifact_contains_metadata_and_schema_only":True}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"n_source_articles":len(results),
         "studies":[{"article":r["article_id"],"files":r["source_file_count"],
                    "all_source_downloaded":r["every_included_original_source_file_downloaded"],
                    "candidate_plant_site_taxon_file_count":r["n_potential_plant_site_taxon_files"],
                    "file_schemas":[{"name":f["file_name"],"metadata":f["schema_without_raw_rows"]}
                                     for f in r["source_files"] if "schema_without_raw_rows" in f]}
                   for r in results]},indent=2,ensure_ascii=False))

if __name__=="__main__":main()
