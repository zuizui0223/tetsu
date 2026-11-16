#!/usr/bin/env python3
"""Small route-label/transition provenance for 2010-25 N05 rail event candidates.

Only factual names, years, and change-code snippets are exported, NOT the
third-party source N05 geometry or station coordinates. Route-end labels
remain CANDIDATE and require independent railway-company verification.
The source data are licensed NONCOMMERCIAL.
"""
from __future__ import annotations
import argparse
from collections import defaultdict,Counter
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import zipfile

from audit_n05_historical_railway import (
    SOURCE,known_year,is_current_end,bounded_zip_download,shape_reader
)
from audit_n05_event_city_overlap import bbox_overlap
from audit_inat_rail_phenophases import AREAS


def classify_candidates(reader,areas=AREAS):
    names=[f[0] for f in reader.fields[1:]]
    required={"N05_004","N05_005b","N05_005e","N05_006","N05_002"}
    if not required.issubset(names):
        raise ValueError(f"Missing line attributes {required-set(names)}")
    groups=defaultdict(lambda:{
        "opens":set(),"ends":set(),"is_current":False,
        "areas":set(),"names":set(),"operators":set(),"revisions":Counter(),
        "n_rows":0
    })
    for pair in reader.iterShapeRecords():
        row=dict(zip(names,pair.record))
        gid=str(row.get("N05_006") or "").strip()
        if not gid:
            continue
        d=groups[gid]
        open_year=known_year(row.get("N05_004"))
        end_year=known_year(row.get("N05_005e"))
        if open_year:d["opens"].add(open_year)
        if end_year:d["ends"].add(end_year)
        if is_current_end(row.get("N05_005e")):d["is_current"]=True
        route=str(row.get("N05_002") or "").strip()
        company=str(row.get("N05_003") or "").strip()
        if route:d["names"].add(route)
        if company:d["operators"].add(company)
        change=str(row.get("N05_007") or "").strip()
        d["revisions"][change]+=1
        d["n_rows"]+=1
        if pair.shape.shapeType != 0 and hasattr(pair.shape,"bbox"):
            for city,b in areas.items():
                if bbox_overlap(pair.shape.bbox,b):
                    d["areas"].add(city)
    candidates=[]
    for group,d in groups.items():
        opened=min(d["opens"]) if d["opens"] else None
        ended=max(d["ends"]) if d["ends"] and not d["is_current"] else None
        events=[]
        if opened is not None and 2010<=opened<=2025:
            events.append({"candidate_type":"documented_open_year",
                           "year":opened})
        if ended is not None and 2010<=ended<=2024:
            events.append({"candidate_type":"terminal_record_UNVERIFIED_closure",
                           "year":ended})
        for e in events:
            candidates.append({
                "group_hash_id":hashlib.sha256(group.encode()).hexdigest()[:12],
                "name_variants":sorted(d["names"])[:4],
                "operator_variants":sorted(d["operators"])[:3],
                "metro_bbox_with_5km_halo":sorted(d["areas"]),
                "n_attribute_versions":d["n_rows"],
                "continues_as_record_2025":d["is_current"],
                "transition_code_histogram":dict(d["revisions"]),
                **e
            })
    return {
        "n_railway_line_relation_groups":len(groups),
        "n_all_candidate_events":len(candidates),
        "n_candidate_types":dict(Counter(c["candidate_type"] for c in candidates)),
        "candidate_event_facts":sorted(
            candidates,key=lambda x:(x["candidate_type"],x["year"],
                                     x["group_hash_id"])
        ),
        "interpretation":"Line names and years are provenance candidates. Terminal version can be an attribute/ID replacement and MUST NOT be called confirmed closure.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/n05_recent_route_names_for_source_verification.json")
    args=p.parse_args()
    raw=bounded_zip_download()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        reader,_=shape_reader(z,"RailroadSection2")
        report=classify_candidates(reader)
        reader.close()
    report.update({
        "source":SOURCE,
        "sha256":hashlib.sha256(raw).hexdigest(),
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "license":"NONCOMMERCIAL N05 2025; contains small derived name/year facts, no licensed coordinates, no source geometry.",
    })
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),
                   encoding="utf-8")
    print(json.dumps({
        "n_railway_line_groups":report["n_railway_line_relation_groups"],
        "candidate_types":report["n_candidate_types"],
        "events_in_pilot_metros":[
            x for x in report["candidate_event_facts"]
            if x["metro_bbox_with_5km_halo"]
        ],
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
