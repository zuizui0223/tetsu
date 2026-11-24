#!/usr/bin/env python3
"""Inspect official EIDC manifest metadata for original plant CSV filename.

No full source attempt; bounded 1 MiB source beginning, no raw species data
reported. May enable direct file retrieval if server cannot stream ZIP.
"""
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
from collections import Counter
import argparse,io,json,re,struct,zipfile
ID="ca339c86-3674-4030-b891-35326e71141e"
URLS={
"support":"https://data-package.ceh.ac.uk/sd/"+ID+".zip",
"catalogue":"https://catalogue.ceh.ac.uk/documents/"+ID+"?type=json",
"original_zip_first":"https://data-package.ceh.ac.uk/data/"+ID+".zip",
}
def get(url,bound=2**20):
  req=Request(url,headers={"User-Agent":"tetsu-v346e-CEH-file-name-audit/1.0"})
  with urlopen(req,timeout=100) as r:
    raw=r.read(bound)
    return {"raw":raw,"status":r.status,"content_type":r.headers.get("Content-Type"),"url":r.geturl()}
def audit():
  o={}
  for key,url in URLS.items():
    try:
      r=get(url,25000 if key=="original_zip_first" else 2**20)
      row={"url":url,"http":r["status"],"content_type":r["content_type"],
           "bytes_probed":len(r["raw"]),"file_paths":[]}
      if key=="original_zip_first":
        b=r["raw"]
        for m in list(re.finditer(b"PK\x03\x04",b))[:10]:
          i=m.start()
          if i+30>len(b):continue
          n,extra=struct.unpack_from("<HH",b,i+26)
          name=b[i+30:i+30+n].decode("utf-8","replace")
          if n<256 and ".." not in name:
            row["file_paths"].append(name)
      elif key=="support":
        with zipfile.ZipFile(io.BytesIO(r["raw"])) as z:
          for name in z.namelist():
            if not name.lower().endswith((".json",".md",".html")):continue
            txt=z.read(name).decode("utf-8","replace")
            row["file_paths"]+=list(set(re.findall(r"[a-zA-Z0-9_./-]+\.(?:csv|tsv|zip|parquet)",txt,re.I)))[:60]
            row["source_document_name_count"]=len(z.namelist())
      elif key=="catalogue":
        txt=r["raw"].decode("utf-8","replace")
        row["file_paths"]+=list(set(re.findall(r"[a-zA-Z0-9_./-]+\.(?:csv|tsv|zip|parquet)",txt,re.I)))[:50]
        try:
          obj=json.loads(txt)
          row["first_level_source_metadata_keys"]=sorted(obj)[:100]
          row["original_eidc_distribution"]=str(obj.get("distribution",""))[:1800]
        except Exception:row["metadata_not_parseable"]=True
      o[key]=row
    except Exception as e:o[key]={"url":url,"error":type(e).__name__,"details":str(e)[:220]}
  return {"source":"https://doi.org/10.5285/"+ID,"probes":o,
          "no_original_species_rows_exported":True}
if __name__=="__main__":
  p=argparse.ArgumentParser();p.add_argument("--output",required=True);a=p.parse_args()
  out=audit();target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
  target.write_text(json.dumps(out,ensure_ascii=False,indent=2))
  print(json.dumps(out,ensure_ascii=False),flush=True)
