#!/usr/bin/env python3
"""Audit first Tokyo 2024 official road centreline mesh (N13).

No plant effect or wind exposure inference. Uses official CC BY 4.0 GIS.
N13 is line geometry, not routing-ready node/link topology; its data
registration date is NOT road opening date.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
import requests
import shapefile

DATASET_PAGE = "https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-2024.html"
N13_URL_PATTERN = "https://nlftp.mlit.go.jp/ksj/gml/data/N13/N13-24/N13-24_{mesh}_SHP.zip"


def audit(mesh="5339"):
    url=N13_URL_PATTERN.format(mesh=mesh)
    s=requests.Session()
    s.headers.update({"User-Agent":"tetsu-rail-research-audit/0.5 (https://github.com/zuizui0223/tetsu)"})
    with s.get(url, stream=True, timeout=(25,180)) as r:
        r.raise_for_status()
        raw=bytearray()
        for chunk in r.iter_content(chunk_size=1024*1024):
            if chunk:
                raw.extend(chunk)
            if len(raw)>150_000_000:
                raise ValueError("Road archive exceeds 150 MB pilot bound")
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError(f"Official road URL gave invalid zip: {url}")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        files=z.namelist()
        shp=[n for n in files if n.lower().endswith(".shp")]
        if not shp:
            raise ValueError(f"No shapefile found in road archive: {files[:5]}")
        # Prefer UTF-8 directory where available, but only inspect, do not
        # convert road features to an active graph.
        preferred=next((x for x in shp if "UTF-8" in x or "utf-8" in x.lower()),shp[0])
        stem=preferred[:-4]
        mem={k.lower():k for k in files}
        def buff(ext):
            return io.BytesIO(z.read(mem[(stem+ext).lower()]))
        reader=shapefile.Reader(shp=buff(".shp"),shx=buff(".shx"),dbf=buff(".dbf"))
        names=[f[0] for f in reader.fields[1:]]
        result={
            "source":url,
            "official_page":DATASET_PAGE,
            "mesh":mesh,
            "sha256":hashlib.sha256(raw).hexdigest(),
            "downloaded_bytes":len(raw),
            "n_road_shapes":len(reader),
            "bbox":reader.bbox,
            "shp_member":preferred,
            "attribute_fields":names,
            "all_members":files,
            "license":"CC BY 4.0 with original source attribution and survey-approval acknowledgement",
            "use_limitation":"Proximity only. Not an active road traffic or network connectivity dataset."
        }
        reader.close()
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mesh",default="5339")
    p.add_argument("--output",default="results/audit/n13_5339_source.json")
    args=p.parse_args()
    result=audit(args.mesh)
    result["audited_at_utc"]=datetime.now(timezone.utc).isoformat()
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k!="all_members"},indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
