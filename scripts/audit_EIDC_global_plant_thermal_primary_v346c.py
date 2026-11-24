#!/usr/bin/env python3
"""Authenticate true EIDC Harris et al global plant climatic envelope original ZIP.

Actual original data at official EIDC /data/<DOI-UUID>.zip endpoint discovered
from primary CEH page. Report exact source ZIP SHA/bytes, CSV column names and
row census, possible T in K and original record-count field names.
No original species-level data rows exported or railway effects calculated.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,os,re,zipfile
from collections import Counter
from pathlib import Path
from urllib.request import Request,urlopen

ID="ca339c86-3674-4030-b891-35326e71141e"
SOURCE_URL="https://data-package.ceh.ac.uk/data/"+ID+".zip"
MAX_ARCHIVE=300_000_000

def download(path):
    h=hashlib.sha256();size=0
    req=Request(SOURCE_URL,headers={"User-Agent":"tetsu-source-verified-climatepref-v346c/1.0"})
    with urlopen(req,timeout=200) as resp,path.open("wb") as out:
        if resp.status!=200:raise ValueError("Original CEH dataset response != 200")
        for block in iter(lambda:resp.read(2**20),b""):
            size+=len(block)
            if size>MAX_ARCHIVE:raise ValueError("Source archive exceeds 300M bound")
            h.update(block);out.write(block)
    return {"bytes":size,"sha256":h.hexdigest(),"url":SOURCE_URL}

def audit(path,meta):
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise ValueError("Original data archive CRC invalid")
        entries=[{"name":x.filename,"uncompressed_bytes":x.file_size} for x in z.infolist() if not x.is_dir()]
        csv_files=[x for x in z.infolist() if not x.is_dir() and x.filename.lower().endswith((".csv",".tsv"))]
        if not csv_files:
            raise ValueError("Original CEH data ZIP did not include expected CSV")
        if len(csv_files)>6:
            raise ValueError("Unclear primary source CSV among >6 original files")
        result=[]
        for f in csv_files:
            with z.open(f) as src:
                wrapper=io.TextIOWrapper(src,encoding="utf-8-sig",errors="replace",newline="")
                delimiter="\t" if f.filename.lower().endswith(".tsv") else ","
                reader=csv.DictReader(wrapper,delimiter=delimiter)
                fields=reader.fieldnames or []
                counts=Counter()
                thermo=[name for name in fields if re.search(r"temp|temperature|tavg|ambient|Tmean",name,re.I)]
                taxa=[name for name in fields if re.search(r"scientific|taxon|species|genus|family",name,re.I)]
                records=[name for name in fields if re.search(r"occurr|records|^n$|^count|sample",name,re.I)]
                examples={}
                temps={}
                for row in reader:
                    counts["rows"]+=1
                    for key in thermo[:16]+records[:8]:
                        value=(row.get(key) or "").strip()
                        if not value:continue
                        try:v=float(value)
                        except (ValueError,OverflowError):continue
                        if not math.isfinite(v):continue
                        k=temps.setdefault(key,{"min":v,"max":v,"valid_n":0})
                        k["min"]=min(k["min"],v);k["max"]=max(k["max"],v)
                        k["valid_n"]+=1
                if counts["rows"]<100000:
                    raise ValueError("Expected 200k+ global species dataset, fewer than 100k rows")
                result.append({"filename":f.filename,"uncompressed_bytes":f.file_size,
                   "record_rows":counts["rows"],"fields":fields,
                   "candidate_temperature_columns":thermo,"candidate_taxonomy_columns":taxa,
                   "candidate_original_record_count_columns":records,
                   "numeric_range_for_candidate_columns":temps,
                   "no_original_taxon_rows_reprinted":True})
        return {"original_zip":meta,"source_file_entries":entries,
                "tables":result,"source_taxon_rows_never_exported":True}

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True)
    a=p.parse_args()
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    archive=dst.parent/"EIDC_global_species_climate_original_2022.zip"
    try:
        meta=download(archive)
        out=audit(archive,meta)
        dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
        print(json.dumps(out,ensure_ascii=False),flush=True)
    finally:
        if archive.exists():archive.unlink()

if __name__=="__main__":main()
