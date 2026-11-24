#!/usr/bin/env python3
"""v360b: source-access probe for available official Geofabrik OSM PBF.

The previously frozen 2026-10-16 daily source returned HTTP404. No OSM
route names, way-memberships, biological outcomes or source records inferred.
Probe predeclared standard monthly snapshot and rolling latest from official
download hosts with HEAD (no 300MB download), plus tiny official md5 sidecar.
"""
import argparse,json,urllib.request
from pathlib import Path
URLS=[
 "https://download.geofabrik.de/europe/slovenia-261001.osm.pbf",
 "https://download.geofabrik.de/europe/slovenia-260901.osm.pbf",
 "https://download.geofabrik.de/europe/slovenia-260101.osm.pbf",
 "https://download.geofabrik.de/europe/slovenia-latest.osm.pbf",
 "https://download-ext2.geofabrik.de/europe/slovenia-latest.osm.pbf",
]
def probe(url):
    row={"official_url":url}
    try:
        req=urllib.request.Request(url,method="HEAD",headers={
            "User-Agent":"tetsu-v360b-source-access-original-rail-route/1.0"})
        with urllib.request.urlopen(req,timeout=28) as r:
            row.update({"http":r.status,"resolved_url":r.geturl(),
              "source_content_length_header":r.headers.get("Content-Length"),
              "source_last_modified_header":r.headers.get("Last-Modified"),
              "source_content_type":r.headers.get("Content-Type")})
    except Exception as e:
        row.update({"HEAD_error_type":type(e).__name__,
           "HEAD_http_status":getattr(e,"code",None),"HEAD_error_message":str(e)[:130]})
    for ext in (".md5",):
        try:
            req=urllib.request.Request(url+ext,headers={"User-Agent":"tetsu-v360b-source-QA/1.0"})
            with urllib.request.urlopen(req,timeout=20) as r:
                raw=r.read(512)
                row["official_md5_http"]=r.status
                row["official_md5_hex"]=raw.decode("utf8","replace")[:100]
        except Exception as e:
            row["official_md5_status_or_error"]=getattr(e,"code",type(e).__name__)
    return row
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);a=p.parse_args()
    out={"version":"v360b-original-Geofabrik-static-vs-latest-source-route",
         "original_Slovenia_oct16_daily_PBF_404_not_a_botanical_failure":True,
         "original_source_file_count_full_railway_route_membership_NOT_YET_AUTHENTICATED":True,
         "results":[probe(u) for u in URLS]}
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf8")
    print(json.dumps(out,ensure_ascii=False),flush=True)
