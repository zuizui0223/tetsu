#!/usr/bin/env python3
"""v376a: authenticated ORIGINAL European Commission JRC GHSL 2020 builtup 1km source gate.

Only source HEAD, 0-4k initial ZIP signature. This does not download
full global archive nor analyze original plant/railway traits. If a
publisher disables byte range or original object cannot be retrieved,
report an honest source failure, not false urbanization control.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import urllib.request,urllib.error

SOURCE_URL=("https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/"
    "GHS_BUILT_S_GLOBE_R2023A/"
    "GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000/V1-0/"
    "GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000_V1_0.zip")
DOI="10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA"
ZIP_MAGIC=b"PK\x03\x04"
MAX_PREFIX=4096

def parse_head(code,headers):
    data={"source_original_head_status":code,
       "source_original_content_length_bytes":None,
       "source_original_etag":headers.get("ETag"),
       "source_last_modified":headers.get("Last-Modified"),
       "original_content_type":headers.get("Content-Type"),
       "original_accept_ranges":headers.get("Accept-Ranges")}
    sz=headers.get("Content-Length")
    if sz is not None:
        try:bytes_total=int(sz)
        except ValueError:raise ValueError("Unparsable publisher Content-Length")
        if bytes_total<=0:raise ValueError("Invalid publisher Content-Length")
        data["source_original_content_length_bytes"]=bytes_total
    return data

def inspect(url=SOURCE_URL):
    headers={"User-Agent":"tetsu-global-macro-independent-JRC-GHSL-source-20261124/1.0"}
    q=urllib.request.Request(url,headers=headers,method="HEAD")
    with urllib.request.urlopen(q,timeout=120) as resp:
        origin=parse_head(resp.status,resp.headers)
        if resp.status!=200:
            raise ValueError("GHSL original publisher HEAD is not HTTP 200")
    # A genuine ZIP is magic-checked on its original path, not inferred
    # from browser HTML query snippets or filename.
    get=urllib.request.Request(url,headers={**headers,"Range":f"bytes=0-{MAX_PREFIX-1}"})
    with urllib.request.urlopen(get,timeout=120) as resp:
        status=resp.status
        raw=resp.read(MAX_PREFIX)
        meta=parse_head(resp.status,resp.headers)
    if status not in (200,206):
        raise ValueError("Publisher original GHSL ZIP GET or Range is neither 200 nor 206")
    if not raw.startswith(ZIP_MAGIC):
        raise ValueError("Original JRC GHSL source did not return ZIP bytes; likely download portal or blocked proxy")
    result={"version":"v376a","original_source":"European Commission JRC Global Human Settlement Layer R2023A",
        "global_independent_builtup_confounding_dataset_DOI":DOI,
        "source_original_global_2020_1km_Mollweide_ZIP_url":url,
        "publisher_head":origin,
        "source_range_first_4096_status":status,
        "source_original_zip_prefix_bytes_checked":len(raw),
        "source_zip_first_bytes_verified":True,
        "source_prefix_SHA256":hashlib.sha256(raw).hexdigest(),
        "original_zip_file_full_SHA256_not_yet_computed":True,
        "source_range_acceptance":status==206,
        "scientific_caveats":[
          "2020 GHS built-up surface 1km values denote built m2 within cell, NOT population, road length, full human footprint or train use",
          "A source ZIP magic/HEAD check is not a complete verified raster acquisition nor a global plant ecological effect",
          "v376b will not claim global urbanization-adjusted railway effects until the full source TIFF and exact 100/200km area sum are both verified",
          "Fine railway vegetation and original plant observations' private coordinates never used"]}
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=inspect()
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
