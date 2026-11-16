#!/usr/bin/env python3
"""v4.9c GBIF occurrenceStatus existential data contract.

Prior 50km source taxon facet count is reproducible ONLY when
occurrenceStatus=PRESENT is omitted (v4.9b), whereas specifying PRESENT
now returns zero. Determine whether source records are actually PRESENT,
ABSENT or unknown. Contrasting source status FILTER vs FACET and
the raw first-page indexed occurrenceStatus field is required before
any event-level botanical interpretation.

No raw records, individual names/IDs or coordinates are exported.
"""
from __future__ import annotations

import json,time
from datetime import datetime,timezone
from pathlib import Path
import requests

from harvest_griis_macro_50km_facets_v31 import species_params
from audit_gbif_equalarea_macroblocks_v29 import make_grid
from audit_gbif_original_row_quality_v48 import source_page

BLOCKS=("us_midwest","fr_burgundy")
PERIOD=(2018,2025)


def run():
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-critical-GBIF-presence-status-index-reconciliation-v4.9c"})
    out={"utc":datetime.now(timezone.utc).isoformat(),
         "source":"GBIF official v1 occurrence API",
         "regions":{},"errors":{}}
    for block in BLOCKS:
        cell=next(r for r in make_grid(block) if r["cell_id"]==f"{block}_x1_y1")
        base=species_params(cell,PERIOD,0)
        base.pop("occurrenceStatus",None)
        base.pop("facet",None)
        base.pop("facetOffset",None)
        base.pop("facetLimit",None)
        base.pop("speciesKey.facetLimit",None)
        base.pop("speciesKey.facetOffset",None)
        versions={
            "no_status_faceted_OCCURRENCE_STATUS":{
                **base,"limit":0,"facet":"OCCURRENCE_STATUS",
                "facetLimit":10
            },
            "no_status_actual_first_300_rows":{**base,"limit":300},
            "filtered_status_PRESENT":{**base,"occurrenceStatus":"PRESENT","limit":0},
            "filtered_status_ABSENT":{**base,"occurrenceStatus":"ABSENT","limit":0},
            "filtered_status_NOT_DETECTED":{**base,"occurrenceStatus":"NOT_DETECTED","limit":0},
        }
        answers={}
        for name,params in versions.items():
            try:
                data=source_page(session,params)
                facet=[]
                for f in data.get("facets",[]) or []:
                    if isinstance(f,dict):
                        facet=[{"label":str(x.get("name")),"indexed_source_record_count":x.get("count")}
                                for x in f.get("counts",[])]
                by_status={}
                for x in data.get("results",[]) or []:
                    status=str(x.get("occurrenceStatus") or "MISSING")
                    by_status[status]=by_status.get(status,0)+1
                answers[name]={
                    "count":data.get("count"),
                    "n_returned_original_rows":len(data.get("results",[]) or []),
                    "source_faceted_status_counts":facet,
                    "raw_first_page_observed_status_categories":by_status,
                    "source_responded":True,
                }
                print(json.dumps({"block":block,"variant":name,
                    **answers[name]},ensure_ascii=False),flush=True)
            except Exception as exc:
                answers[name]={"source_responded":False,
                    "error":f"{type(exc).__name__}: {str(exc)[:180]}"}
                out["errors"][block+"_"+name]=answers[name]["error"]
            time.sleep(.1)
        out["regions"][block]={"source_cell":cell["cell_id"],"country":cell["country"],
                               "results":answers}
    out["complete"]=not out["errors"]
    out["limits"]=[
        "Only two pre-frozen 50km geographic source cells sampled for API semantics; not a countrywide study.",
        "A GBIF status facet is source data but not independent botanical evidence.",
        "A no-filter FIRST 300 occurrence-row page is not a random sample when total >300.",
        "No source-level railway seed dispersal, observer quality or field survey absence inferred.",
    ]
    return out


def main():
    result=run()
    path=Path("results/audit/gbif_occurrence_status_API_contract_v49c.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":result["complete"],"errors":result["errors"]},
                     ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("GBIF occurrence status source contract failed")

if __name__=="__main__":
    main()
