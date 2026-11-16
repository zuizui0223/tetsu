#!/usr/bin/env python3
"""Coarse geographic feasibility of historical N05 route events near fixed JP metros.

All railway open/end groups selected by official temporal *attributes*,
NOT by floral observations or effect direction. City bounding boxes were
frozen in tetsu v0.4; a 5-km halo avoids edge misclassification. Results
are aggregate counts only, not original noncommercial N05 line geometry.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import zipfile

from audit_n05_historical_railway import (
    SOURCE, INFO, known_year, is_current_end,
    bounded_zip_download, shape_reader
)
from audit_inat_rail_phenophases import AREAS


def bbox_overlap(shape_bbox, area, padding=0.05):
    """SHP bbox x_min, y_min, x_max, y_max vs source lat/lon rectangle."""
    lat0, lon0, lat1, lon1 = area
    ax0,ay0,ax1,ay1=shape_bbox
    return not (
        ax1 < lon0-padding or ax0 > lon1+padding
        or ay1 < lat0-padding or ay0 > lat1+padding
    )


def summarize_n05_city_event_support(reader, areas=AREAS):
    names=[x[0] for x in reader.fields[1:]]
    required={"N05_004","N05_005b","N05_005e","N05_006"}
    if not required.issubset(names):
        raise ValueError(f"N05 missing fields: {required-set(names)}")
    groups=defaultdict(lambda:{
        "opening_years":set(), "ends":set(),
        "has_current":False, "metro":set(),
        "n_records":0
    })
    bad_geometry=0
    for item in reader.iterShapeRecords():
        row=dict(zip(names,item.record))
        group=str(row.get("N05_006") or "").strip()
        if not group:
            continue
        record=groups[group]
        record["n_records"]+=1
        opening=known_year(row.get("N05_004"))
        ending=known_year(row.get("N05_005e"))
        if opening:
            record["opening_years"].add(opening)
        if ending:
            record["ends"].add(ending)
        if is_current_end(row.get("N05_005e")):
            record["has_current"]=True
        geometry=item.shape
        if geometry.shapeType in (0,):
            bad_geometry+=1
            continue
        if not hasattr(geometry,"bbox"):
            bad_geometry+=1
            continue
        for name, rect in areas.items():
            if bbox_overlap(geometry.bbox,rect):
                record["metro"].add(name)
    cohort=Counter()
    by_year=defaultdict(Counter)
    current_hist=Counter()
    for group,info in groups.items():
        recent_new=min(info["opening_years"]) if info["opening_years"] else None
        tentative_end=max(info["ends"]) if info["ends"] and not info["has_current"] else None
        for metro in info["metro"]:
            current_hist[metro]+=int(info["has_current"])
            if recent_new is not None and 2010<=recent_new<=2025:
                cohort[(metro,"open_2010_2025")]+=1
                by_year[(metro,"opened")][recent_new]+=1
            if tentative_end is not None and 2010<=tentative_end<=2024:
                cohort[(metro,"possible_terminal_2010_2024")]+=1
                by_year[(metro,"potential_end")][tentative_end]+=1
    return {
        "type":"coarse bbox intersection of passenger N05 history relation IDs, 5km halo",
        "n_relation_groups":len(groups),
        "n_rail_history_records_missing_geometry":bad_geometry,
        "n_current_relation_groups_overlapping_city_bbox_halo":{
            name:current_hist[name] for name in areas
        },
        "city_event_candidate_groups":{
            name:{
                "opening_known_2010_2025":cohort[(name,"open_2010_2025")],
                "possible_terminal_not_yet_verified_2010_2024":
                    cohort[(name,"possible_terminal_2010_2024")],
                "open_year_distribution":{
                    str(y):n for y,n in sorted(by_year[(name,"opened")].items())
                },
                "possible_terminal_year_distribution":{
                    str(y):n for y,n in sorted(by_year[(name,"potential_end")].items())
                },
            } for name in areas
        },
        "independence_caveat":"One N05 group is a GIS relation grouping, not a whole railway corridor or independent city; city boundaries may overlap a line group.",
        "event_identification_caveat":"Current endpoint=9999 means a continuing GIS record; a terminal recorded year may still reflect replacement of line/group, not a verified physical closure.",
        "botanical_effect_caveat":"City intersections only: no iNaturalist station-level incidence, no growth or seed results, no rail-wind or mowing exposure; rail opening can shift sampling and traffic."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/n05_recent_events_metropolitan_overlap.json")
    args=p.parse_args()
    archive=bounded_zip_download()
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        reader,member=shape_reader(z,"RailroadSection2")
        result=summarize_n05_city_event_support(reader)
        reader.close()
    result.update({
        "source":SOURCE,
        "source_sha256":hashlib.sha256(archive).hexdigest(),
        "source_shapefile_name":Path(member).name,
        "computed_utc":datetime.now(timezone.utc).isoformat(),
        "license":"MLIT N05 noncommercial-only; original geometry NOT published",
        "place_rectangles":{k:list(v) for k,v in AREAS.items()}
    })
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(result["city_event_candidate_groups"],ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
