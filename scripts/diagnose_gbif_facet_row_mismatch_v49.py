#!/usr/bin/env python3
"""v4.9: critical GBIF count-vs-row source reconciliation.

Within 2 fixed outcome-blind 50km source cells US Midwest and France
Burgundy, compare exact same occurrence filters with:
- SOURCE FACET + limit=0 (previous GRIIS speciesKey census)
- SOURCE FACET + limit=1
- NO FACET + limit=0, 1, 300
- only country/kingdom/year + limit=1 as broad positive control.

Do not inspect species, railway effects, or original record names.
This is a source API failure/contract diagnostic; row counts and
geography must AGREE before treating 2026 GBIF facets as real incidence.
"""
from __future__ import annotations
import json,sys
from datetime import datetime,timezone
from pathlib import Path
import requests
import time

from audit_gbif_equalarea_macroblocks_v29 import make_grid
from audit_gbif_original_row_quality_v48 import source_params as row_params
from harvest_griis_macro_50km_facets_v31 import species_params

LOCATIONS=("us_midwest","fr_burgundy")
PERIOD=(2018,2025)
API="https://api.gbif.org/v1/occurrence/search"


def try_request(session,name,params):
    try:
        r=session.get(API,params=params,timeout=(20,70))
        r.raise_for_status()
        data=r.json()
        return {
            "label":name,
            "http":r.status_code,
            "count":data.get("count"),
            "n_results":len(data.get("results") or []),
            "endOfRecords":data.get("endOfRecords"),
            "n_facets":len(data.get("facets") or []),
            "api_url_param_names":sorted(params),
            "returned_keys":sorted(data.keys()),
            "source_responded":True,
        }
    except Exception as exc:
        return {"label":name,"source_responded":False,
                "error":f"{type(exc).__name__}: {str(exc)[:220]}"}


def run():
    s=requests.Session()
    s.headers.update({"User-Agent":"tetsu-critical-GBIF-facet-vs-real-occurrence-reconcile-v4.9"})
    out={"utc":datetime.now(timezone.utc).isoformat(),
         "fixes_exact_source_cell_before_examining_results":True,
         "regions":{},"source_inconsistent":False}
    for block in LOCATIONS:
        cell=next(x for x in make_grid(block) if x["cell_id"]==f"{block}_x1_y1")
        params=row_params(cell,"2018_2025",limit=1)
        from_source=species_params(cell,PERIOD,offset=0)
        from_source["limit"]=0
        versions=[
            ("source_species_facet_limit0",from_source),
            ("source_species_facet_limit1",{**from_source,"limit":1}),
            ("original_row_no_facet_limit0",{**params,"limit":0}),
            ("original_row_no_facet_limit1",params),
            ("original_row_no_facet_limit300",{**params,"limit":300}),
            ("broad_country_control",{
                "country":cell["country"],
                "kingdomKey":6,
                "year":"2018,2025",
                "hasCoordinate":"true",
                "limit":1
            }),
        ]
        answers=[]
        for label,p in versions:
            result=try_request(s,label,p)
            answers.append(result)
            print(json.dumps({"block":block,"variant":label,"count":result.get("count"),
                "n_results":result.get("n_results"),"http":result.get("http"),
                "source_responded":result["source_responded"],
                "error":result.get("error")},ensure_ascii=False),flush=True)
            time.sleep(.12)
        a={x["label"]:x for x in answers}
        required=[a[k] for k in ("source_species_facet_limit0","original_row_no_facet_limit1",
                                  "original_row_no_facet_limit300","broad_country_control")]
        inconsistent=(any(not v["source_responded"] for v in required)
                      or len({v.get("count") for v in required[:3]})>1
                      or (required[1].get("count",0)>0 and required[1].get("n_results")==0))
        out["source_inconsistent"]=out["source_inconsistent"] or inconsistent
        out["regions"][block]={"country":cell["country"],
                              "cell":cell["cell_id"],
                              "variant_responses":answers,
                              "inconsistent":inconsistent}
    out["complete"]=not out["source_inconsistent"]
    out["interpretation"]="If facet limit0 source and raw record source disagree, previous source-level taxon count cannot be interpreted as verified occurrence-row population."
    return out


def main():
    result=run()
    path=Path("results/audit/gbif_facet_vs_rows_v49.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"source_consistent":result["complete"],
                      "source_inconsistent":result["source_inconsistent"]},ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("Critical GBIF count-vs-row API mismatch; old biological source inference suspended")
if __name__=="__main__":
    main()
