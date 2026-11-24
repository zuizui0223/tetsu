#!/usr/bin/env python3
"""v342 original-source access gate only: genuine Australian rail comparison,
GBIF API and WorldClim temperature reference. Never emits source plot rows.
"""
from __future__ import annotations
import argparse,hashlib,io,json,sys,zipfile
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from openpyxl import load_workbook

SUPPLEMENT_URLS=[
    "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Faec.12719&file=aec12719-sup-0002-AppendixS4.xlsx",
    "https://onlinelibrary.wiley.com/doi/suppl/10.1111/aec.12719/supinfo/aec12719-sup-0002-AppendixS4.xlsx",
]
GBIF="https://api.gbif.org/v1"
CLIMATE_URL="https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip"
MAX_SOURCE_BYTES=80000000
ANCHORS=("Cynodon dactylon","Poa annua")
AGENT="tetsu-research-source-intake-v342/1.0 (source-verification, no mass harvest)"

def fetch(url,timeout=35,maxbytes=MAX_SOURCE_BYTES):
    req=Request(url,headers={"User-Agent":AGENT,"Accept":"application/octet-stream, application/json, */*"})
    with urlopen(req,timeout=timeout) as response:
        mediatype=response.headers.get("Content-Type","")
        raw=response.read(maxbytes+1)
        if len(raw)>maxbytes:
            raise ValueError("Source content exceeds source gate maxbytes; no partial file used")
        return {"content":raw,"media_type":mediatype,"status":getattr(response,"status",200),
                "final_url":response.geturl()}

def summarize_xlsx(raw):
    if not raw.startswith(b"PK\\x03\\x04") and not raw.startswith(b"PK\x03\x04"):
        raise ValueError("Not genuine XLSX ZIP signature")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        if not all(k in archive.namelist() for k in ("[Content_Types].xml","xl/workbook.xml")):
            raise ValueError("Not XLSX workbook; probable publisher error page")
    book=load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    try:
        sheets=[]
        for sheet in book.worksheets:
            h=next(sheet.values,None)
            sheets.append({"sheet_title":sheet.title,"row_count_incl_header":sheet.max_row,
                          "columns":sheet.max_column,
                          "first_row_field_names":[str(v)[:70] if v is not None else None for v in (h or ())][:24]})
        return {"sha256":hashlib.sha256(raw).hexdigest(),"size_bytes":len(raw),
                "sheets":sheets}
    finally:book.close()

def safe_fetch(url, timeout=35,maxbytes=MAX_SOURCE_BYTES):
    try:
        return fetch(url,timeout,maxbytes),None
    except (HTTPError,URLError,TimeoutError,ValueError) as e:
        code=getattr(e,"code",None)
        return None,{"type":type(e).__name__,"http_status":code,"short_message":str(e)[:180]}

def australian_gate():
    attempts=[]
    verified=None
    for url in SUPPLEMENT_URLS:
        response,error=safe_fetch(url,timeout=40,maxbytes=3000000)
        row={"requested_official_url":url,
             "publisher_HTTP_retrieval_success":response is not None,
             "http_status":response["status"] if response else error.get("http_status"),
             "error":error}
        if response:
            try:
                original=summarize_xlsx(response["content"])
                row["xlsx_source_verified"]=True
                row["original_xlsx_schema_without_source_rows"]=original
                verified=original
            except (ValueError,zipfile.BadZipFile) as exc:
                row["xlsx_source_verified"]=False
                row["schema_validation_error"]=str(exc)[:180]
        attempts.append(row)
        if verified:break
    return {"study":"Everingham Hemmings Moles 2019 doi:10.1111/aec.12719",
            "publisher_publishes_Appendix_S4":True,
            "actual_supplement_source_xlsx_retrieved_and_verified":bool(verified),
            "publisher_retrieval_attempts":attempts,
            "source_comparison_is_field_disturbance_gradient_NOT_assumed_same_site_paired":True,
            "matched_true_rail_nonrail_source_plot_matrix_authenticated":False}

def gbif_taxon_gate():
    results=[]
    for name in ANCHORS:
        url=f"{GBIF}/species/match?"+urlencode({"name":name,"strict":"true","kingdom":"Plantae"})
        response,error=safe_fetch(url,timeout=25,maxbytes=200000)
        out={"source_exact_name":name,"api_name_match_ok":False}
        if response:
            try:
                js=json.loads(response["content"])
                key=js.get("usageKey")
                out.update({"matchType":js.get("matchType"),
                            "api_status":js.get("status"),
                            "rank":js.get("rank"),
                            "scientificName":js.get("scientificName"),
                            "acceptedUsageKey":js.get("acceptedUsageKey"),
                            "name_match_usageKey":key,
                            "api_name_match_ok":bool(key) and js.get("matchType")=="EXACT"})
                if key and out["api_name_match_ok"]:
                    count_url=f"{GBIF}/occurrence/search?"+urlencode({
                        "taxonKey":key,"hasCoordinate":"true",
                        "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
                        "limit":1})
                    rc,err=safe_fetch(count_url,timeout=35,maxbytes=100000)
                    if rc:
                        docs=json.loads(rc["content"])
                        out["coordinate_quality_filtered_total_count"]=docs.get("count")
                        out["source_fields_of_first_one_record_NOT_records"]=(sorted(docs["results"][0])[:50]
                          if docs.get("results") else [])
                    else:
                        out["occurrence_api_error"]=err
            except (ValueError,KeyError,TypeError) as exc:
                out["error"]={"type":type(exc).__name__,"short_message":str(exc)[:120]}
        else:out["error"]=error
        results.append(out)
    return {"source":"GBIF taxon match and occurrence/search REST API",
            "species_anchor_requests":results,
            "api_is_only_connectivity_gate_not_global_species_niches":True,
            "gbif_occurrence_search_100k_offset_plus_limit_cap":True,
            "full_bulk_download_authentication_and_DOI_required":True}

def summarize_worldclim_zip(raw):
    if not raw.startswith(b"PK\x03\x04"):
        raise ValueError("WorldClim not ZIP")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        tif=[f for f in z.namelist() if f.lower().endswith(".tif")]
        biotemperature=any("_bio_1.tif" in x.lower() for x in tif)
        if len(tif)<19 or not biotemperature:
            raise ValueError("No official WorldClim bio1 TIFF or fewer than 19 bioclim rasters")
        return {"sha256":hashlib.sha256(raw).hexdigest(),
                "archive_zip_bytes":len(raw),"bio_tif_count":len(tif),
                "has_annual_mean_temperature_bio1_tif":True,
                "first_4_source_tif_names":sorted(tif)[:4]}

def worldclim_gate():
    response,error=safe_fetch(CLIMATE_URL,timeout=120,maxbytes=MAX_SOURCE_BYTES)
    item={"source":CLIMATE_URL,"historical_period":"1970-2000",
          "native_raster_is_not_ballast_microclimate":True,
          "download_success":response is not None}
    if error:item["error"]=error
    if response:
        try:
            item["source_tif_schema_verified"]=summarize_worldclim_zip(response["content"])
        except (ValueError,zipfile.BadZipFile) as exc:
            item["error"]={"type":type(exc).__name__,"short_message":str(exc)[:140]}
    return item

def analyze():
    return {"version":"v342-source-access-only",
            "attempts":{
                "Everingham_Australia":australian_gate(),
                "GBIF":gbif_taxon_gate(),
                "WorldClim":worldclim_gate()
            },
            "no_species_climate_means_or_Everingham_field_effect_claimed":True,
            "does_not_republish_field_species_rows_or_GPS":True,
            "independent_climate_region_rail_v_nonrail_replication_estimated":False}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=analyze()
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)
