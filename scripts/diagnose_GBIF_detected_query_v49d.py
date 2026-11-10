#!/usr/bin/env python3
"""v4.9d reconcile GBIF DETECTED / NOT_DETECTED source filters.

GBIF current OCCURRENCE_STATUS FACET reports DETECTED and NOT_DETECTED,
while legacy occurrenceStatus=PRESENT search returns zero. Determine
whether current search supports DETECTED/NOT_DETECTED directly, and
inspect ONLY field names (not values of personal identifiers) from
actual row search to locate the status schema representation.

No raw human observer, occurrence key, or coordinate retained.
"""
from __future__ import annotations
import json,time
from datetime import datetime,timezone
from pathlib import Path
import requests
from audit_gbif_equalarea_macroblocks_v29 import make_grid
from harvest_griis_macro_50km_facets_v31 import species_params
from audit_gbif_original_row_quality_v48 import source_page

LOCATIONS=("us_midwest","fr_burgundy")


def run():
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-GBIF-DETECTED-query-support-check-v4.9d"})
    out={"utc":datetime.now(timezone.utc).isoformat(),"regions":{}}
    for block in LOCATIONS:
        cell=next(v for v in make_grid(block) if v["cell_id"]==f"{block}_x1_y1")
        p=species_params(cell,(2018,2025),0)
        for key in ("facet","facetLimit","facetOffset",
                    "speciesKey.facetLimit","speciesKey.facetOffset","occurrenceStatus"):
            p.pop(key,None)
        cases={
            "all_status":{**p,"limit":5},
            "status_DETECTED":{**p,"occurrenceStatus":"DETECTED","limit":5},
            "status_NOT_DETECTED":{**p,"occurrenceStatus":"NOT_DETECTED","limit":5},
            "status_PRESENT_legacy":{**p,"occurrenceStatus":"PRESENT","limit":5},
            "status_ABSENT_legacy":{**p,"occurrenceStatus":"ABSENT","limit":5},
        }
        vals={}
        for label,params in cases.items():
            try:
                data=source_page(session,params)
                first=(data.get("results") or [{}])[0]
                if not isinstance(first,dict):first={}
                keypaths=[k for k in first if any(z in k.lower() for z in (
                    "occurrence","status","detect","classification","taxon"))]
                vals[label]={
                    "http_success":True,
                    "source_indexed_count":data.get("count"),
                    "n_original_indexed_records_returned":len(data.get("results") or []),
                    "nonidentifying_status_related_fields_in_first_row":sorted(keypaths),
                    "status_like_record_enum_values":{
                        k:str(first[k]) for k in keypaths
                        if "status" in k.lower() or "detect" in k.lower()
                    },
                }
            except requests.HTTPError as exc:
                vals[label]={"http_success":False,
                             "http_status":exc.response.status_code,
                             "error_type":"HTTP_ERROR"}
            except Exception as exc:
                vals[label]={"http_success":False,
                             "error_type":type(exc).__name__}
            print(json.dumps({"block":block,"query":label,
                              **vals[label]},ensure_ascii=False),flush=True)
            time.sleep(.12)
        out["regions"][block]=vals
    out["complete"]=all(
        v["all_status"]["http_success"] for v in out["regions"].values()
    )
    out["source_limitation"]="Only fixed 50km US/France source cells, explicit first five raw records; not a representative survey or seed transport signal."
    return out


def main():
    result=run()
    file=Path("results/audit/gbif_DETECTED_search_filter_support_v49d.json")
    file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":result["complete"],
                      "status_search_behaviors":{k:{x:[y.get("http_success"),y.get("source_indexed_count")] for x,y in v.items()}
                      for k,v in result["regions"].items()}},ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("GBIF source baseline did not answer properly")


if __name__=="__main__":
    main()
