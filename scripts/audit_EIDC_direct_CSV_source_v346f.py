#!/usr/bin/env python3
"""Strictly test EIDC official direct original CSV routes discovered from ZIP.

Only first 1KiB sampled; don't mirror raw taxon rows or treat an HTML page as
a CSV. Full ZIP is non-resumable dynamic chunked and dropped after 3min.
"""
import argparse,json
from pathlib import Path
from urllib.request import Request,urlopen
UUID="ca339c86-3674-4030-b891-35326e71141e"
FILE="global_plant_climate_envelopes.csv"
URLS=[
 f"https://catalogue.ceh.ac.uk/datastore/eidchub/{UUID}/{FILE}",
 f"https://catalogue.ceh.ac.uk/datastore/eidchub/{UUID}/data/{FILE}",
 f"https://data-package.ceh.ac.uk/data/{UUID}/data/{FILE}",
 f"https://data-package.ceh.ac.uk/data/{UUID}/{FILE}",
 f"https://catalogue.ceh.ac.uk/datastore/eidchub/{UUID}/",
]
def probe(url):
  try:
    with urlopen(Request(url,headers={"User-Agent":"tetsu-v346f-CEH-original-direct-file-test/1.0",
         "Range":"bytes=0-4095"}),timeout=28) as response:
      blob=response.read(4096)
      first=blob[:512].decode("utf-8","replace")
      names=first.split("\n",1)[0][:512]
      return {"url":url,"status":response.status,
        "content_type":response.headers.get("Content-Type"),
        "content_length":response.headers.get("Content-Length"),
        "content_range":response.headers.get("Content-Range"),
        "accept_ranges":response.headers.get("Accept-Ranges"),
        "csv_header_candidate":names if ("species" in names.lower()
          or "spp" in names.lower() or "temp" in names.lower()) else None,
        "body_prefix_is_HTML":blob[:64].lstrip().lower().startswith((b"<html",b"<!doctype")),
        "sampled_bytes":len(blob)}
  except Exception as e:
    return {"url":url,"error_type":type(e).__name__,
      "http":getattr(e,"code",None),"short_error":str(e)[:150]}
def run():
  p=[probe(url) for url in URLS]
  return {"official_source_DOI":"10.5285/"+UUID,"original_CSV_filename_verified_from_EIDC_ZIP_manifest":FILE,
    "candidate_direct_CEH_file_routes":p,
    "original_species_rows_not_republished":True,
    "global_railway_warm_affinity_effect_NOT_estimated":True}
if __name__=="__main__":
 a=argparse.ArgumentParser();a.add_argument("--output",required=True);x=a.parse_args()
 out=run();path=Path(x.output);path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(out,ensure_ascii=False,indent=2))
 print(json.dumps(out,ensure_ascii=False),flush=True)
