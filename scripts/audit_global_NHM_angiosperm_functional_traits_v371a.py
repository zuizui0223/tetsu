#!/usr/bin/env python3
"""v371a: audit genuine independently worldwide angiosperm FUNCTIONAL traits.

Natural History Museum 2026 ADEPT/WFO accepted angiosperms trait resource,
DOI 10.5519/p3dm31kc, claimed 59,369 source species and 888,328 parsed
trait values from global botanical literature. It is NOT 59,369 species
with complete growth form and dispersal data. Inspect source resource
real downloads + exact source CSV schema/trait column completeness BEFORE
joining original iNaturalist worldwide flora.

Source NHM CKAN API metadata only plus 256k sample of combined.csv, never
individual original taxon names or raw source trait text redistributed.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,urllib.error,urllib.parse,urllib.request
from collections import Counter
from pathlib import Path
DOI="10.5519/p3dm31kc"
CKAN="https://data.nhm.ac.uk/api/3/action/package_show?id=angiosperm-functional-traits-extracted-from-botanical-descriptions"
DATASET="angiosperm-functional-traits-extracted-from-botanical-descriptions"
SOURCE_YEAR=2026
PREFIX_BYTES=262144

def http(url,limit=2_000_000,method="GET",range_bytes=None):
    hdr={"User-Agent":"tetsu-v371a-global-source-botanical-functional-macro/1.0",
         "Accept":"application/json, text/csv;q=0.8, */*;q=0.2"}
    if range_bytes is not None:hdr["Range"]=f"bytes=0-{range_bytes-1}"
    req=urllib.request.Request(url,headers=hdr,method=method)
    with urllib.request.urlopen(req,timeout=75) as resp:
        raw=resp.read(limit+1)
        if len(raw)>limit:raw=raw[:limit]
        return raw,{"status":resp.status,"content_length_header":resp.headers.get("Content-Length"),
            "etag":resp.headers.get("ETag"),"last_modified":resp.headers.get("Last-Modified"),
            "content_type":resp.headers.get("Content-Type"),"resolved_URL":resp.url}

def ckan():
    raw,meta=http(CKAN,limit=6_000_000)
    obj=json.loads(raw)
    if obj.get("success") is not True or not isinstance(obj.get("result"),dict):
        raise ValueError("Not an original NHM Data Portal CKAN package source")
    result=obj["result"]
    if result.get("name")!=DATASET:
        raise ValueError("Original official NHM global trait dataset ID changed")
    resources=result.get("resources") or []
    if not resources:raise ValueError("Official NHM dataset resources missing")
    summaries=[]
    for r in resources:
        summaries.append({"name":r.get("name"),"id":r.get("id"),"format":r.get("format"),
          "size":r.get("size"),"hash":r.get("hash"),
          "last_modified":r.get("last_modified"),
          "url":r.get("url"),"datastore_active":r.get("datastore_active"),
          "created":r.get("created")})
    return result,summaries,meta

def source_sample(resource):
    name=str(resource.get("name") or "").strip().lower()
    url=resource.get("url")
    if not url or not url.startswith("https://"):
        raise ValueError("Original NHM resource URL absent or not HTTPS")
    try:
        data,meta=http(url,limit=PREFIX_BYTES,range_bytes=PREFIX_BYTES)
        lines=data.decode("utf-8-sig",errors="strict").splitlines()
        if len(lines)<3:raise ValueError("Original NHM combined.csv contains no CSV header")
        row=next(csv.reader([lines[0]]))
        data_sha=hashlib.sha256(data).hexdigest()
        source_prefix_len=len(data)
        original_source_type="Original direct NHM CSV HTTP source"
        datastore_total=None
    except urllib.error.HTTPError as err:
        if err.code not in (403,404):raise
        if not resource.get("id"):raise ValueError("Original NHM resource ID missing")
        api="https://data.nhm.ac.uk/api/3/action/datastore_search?"+urllib.parse.urlencode(
            {"resource_id":resource["id"],"limit":1})
        raw,meta=http(api,limit=3_000_000)
        obj=json.loads(raw)
        if not obj.get("success"):raise ValueError("Publisher DataStore inaccessible")
        result=obj.get("result") or {}
        row=[x["id"] for x in result.get("fields",[]) if x.get("id") not in ("_id",)]
        if not row:raise ValueError("Publisher NHM DataStore has no actual columns")
        data_sha=hashlib.sha256(raw).hexdigest()
        source_prefix_len=len(raw)
        original_source_type="Official NHM CKAN read-only datastore_search metadata, CSV direct returned 403"
        datastore_total=result.get("total")
    if len(row)<4 or len(row)!=len(set(row)):
        raise ValueError("Original NHM global trait source malformed columns")
    trait_fields=[c for c in row if any(k in c.lower() for k in
       ("growth","life","habit","dispers","seed","fruit","plant_height",
        "height","leaf","woody","annual","duration"))]
    return {"original_source_resource":name,
     "original_NHM_CSV_header_columns":row,
     "original_NHM_source_column_count":len(row),
     "original_functional_candidate_source_columns":trait_fields,
     "source_prefix_length_bytes":source_prefix_len,
     "publisher_official_datastore_total_records_reported":datastore_total,\n     "source_read_method":original_source_type,
     "source_original_prefix_SHA256":data_sha,
     "source_HTTP":{k:v for k,v in meta.items() if k!="resolved_URL"},
     "original_species_names_and_traits_SAMPLE_values_NOT_republished":True}

def run():
    pkg,res,meta=ckan()
    rec={x["name"].lower().strip():x for x in res if isinstance(x.get("name"),str)}
    primary=next((r for k,r in rec.items() if k=="combined.csv"),None)
    world=next((r for k,r in rec.items() if k=="wfo-angiosperms.csv"),None)
    if primary is None or world is None:
        raise ValueError("Original 2026 NHM global angiosperm trait and WFO source files not found")
    sample=source_sample(primary)
    if not sample["original_functional_candidate_source_columns"]:
        raise ValueError("Original NHM source combined.csv has no obvious functional trait columns; inspect schema")
    return {"version":"v371a-original-NHM-global-WFO-ADEPT-angiosperm-functional-trait-source-schema",
      "source_DOI":DOI,
      "source_original_publisher":"Natural History Museum London, Araujo Hartley Brummit 2026",
      "source_original_source_author_CC_BY_not_OpenAI_claim":pkg.get("license_title"),
      "source_official_NHM_API_resource_count":len(res),
      "source_official_WFO_angiosperms_declared_records":342625,
      "source_official_all_source_trait_species_declared":59369,
      "source_official_trait_values_declared":888328,
      "source_original_maintainer_modification_time":pkg.get("metadata_modified"),
      "original_full_global_trait_source_resources":res,
      "original_global_WFO_accepted_angiosperm_resource":world,
      "original_combined_global_traits_header_QA":sample,
      "original_taxonomic_join_WFO_species_to_iNat_species_not_verified":True,
      "original_source_trait_coverage_will_be_measured_not_all_imputed":True,
      "original_plant_species_with_trait_missing_NOT_assigned_unknown_to_other":True,
      "no_railway_microhabitat_or_CAUSAL_plants_effect_estimated":True}
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--output",required=True)
 a=p.parse_args();res=run()
 o=Path(a.output);o.parent.mkdir(parents=True,exist_ok=True)
 o.write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding="utf8")
 print(json.dumps(res,ensure_ascii=False),flush=True)
