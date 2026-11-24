#!/usr/bin/env python3
"""Test whether EIDC official source supports resumable HTTP ranges.

Do not download whole archive. Only response headers and first ~4 KiB.
"""
import json,argparse
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
URL="https://data-package.ceh.ac.uk/data/ca339c86-3674-4030-b891-35326e71141e.zip"
REQUESTS=(("HEAD",{}),("GET",{"Range":"bytes=0-4095"}),
          ("GET",{"Range":"bytes=1048576-1052671"}))
def probe(method,hdr):
  out={"method":method,"requested_range":hdr.get("Range")}
  try:
    req=Request(URL,method=method,headers={"User-Agent":"tetsu-global-thermal-eidc-http-range-diagnostic/1.0",**hdr})
    with urlopen(req,timeout=55) as r:
      data=r.read(4096) if method=="GET" else b""
      out.update({"status":r.status,"content_length":r.headers.get("Content-Length"),
        "content_range":r.headers.get("Content-Range"),"accept_ranges":r.headers.get("Accept-Ranges"),
        "transfer_encoding":r.headers.get("Transfer-Encoding"),
        "content_type":r.headers.get("Content-Type"),"read_bytes":len(data),
        "prefix_hex":data[:12].hex()})
  except Exception as e:
    out.update({"error":type(e).__name__,"http_status":getattr(e,"code",None),
                "short_message":str(e)[:180]})
  return out
if __name__=="__main__":
  a=argparse.ArgumentParser();a.add_argument("--output",required=True);x=a.parse_args()
  r={"original_EIDC_URL":URL,"probes":[probe(*k) for k in REQUESTS],
      "headers_only_not_empirical_species_effect":True}
  p=Path(x.output);p.parent.mkdir(parents=True,exist_ok=True)
  p.write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
