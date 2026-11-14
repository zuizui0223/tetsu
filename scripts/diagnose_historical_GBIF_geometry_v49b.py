#!/usr/bin/env python3
"""v4.9b mandatory provenance-vs-current GBIF 50km query source diagnostic.

The previous v3.6 file (source #42761140912) reports indexed COUNT >0
in several pre-frozen 50km cells. Today both faceted limit=0 and raw
limit=1 sometimes give COUNT=0 for exactly those same source cells.
This is a critical contradiction; diagnose all query dimensions, not
a biological plant extinction. No species/rail results are computed.
"""
from __future__ import annotations
from datetime import datetime,timezone
import json,sys,time
from pathlib import Path

import requests
from shapely import wkt
from shapely.geometry import box
from shapely.geometry.polygon import orient

from audit_gbif_equalarea_macroblocks_v29 import make_grid
from harvest_griis_macro_50km_facets_v31 import species_params
from audit_gbif_original_row_quality_v48 import source_page

REGIONS=("us_midwest","fr_burgundy")
SOURCE_YEAR=(2018,2025)


def source_rows(cell):
    params=species_params(cell,SOURCE_YEAR,0)
    return {
        "full_current_50km_facet":params,
        "without_geographic_polygon":{k:v for k,v in params.items() if k!="geometry"},
        "without_year":{k:v for k,v in params.items() if k!="year"},
        "without_kingdomKey":{k:v for k,v in params.items() if k!="kingdomKey"},
        "without_occurrenceStatus":{k:v for k,v in params.items() if k!="occurrenceStatus"},
        "without_taxonRank":{k:v for k,v in params.items() if k!="taxonRank"},
        "without_hasGeospatialIssue":{k:v for k,v in params.items() if k!="hasGeospatialIssue"},
        "minimal_geometry_country_year":{k:params[k] for k in ("country","geometry","year","limit")},
        "minimal_geometry_country_only":{k:params[k] for k in ("country","geometry","limit")},
    }


def run(flora):
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-geospatial-GBIF-index-contradiction-debug-v4.9b"})
    report={"utc":datetime.now(timezone.utc).isoformat(),
            "previous_real_GRIIS_source_run":42761140912,
            "cell_locations_fixed_from_existing_v2_9_source_not_outcomes":True,
            "regions":{},"errors":{}}
    for block in REGIONS:
        row=next(r for r in make_grid(block) if r["cell_id"]==f"{block}_x1_y1")
        p=flora["regions"][block]["cells"][row["cell_id"]]
        previous=p["periods"]["2018_2025"]
        prev=previous.get("indexed_all_plantae_source_records")
        cases=source_rows(row)
        shp=wkt.loads(row["gbif_geometry_wgs84"])
        cases["simple_geometry_bbox"]={**cases["full_current_50km_facet"],
                                      "geometry":box(*shp.bounds).wkt}
        cases["opposite_winding_dense_polygon"]={**cases["full_current_50km_facet"],
                                             "geometry":orient(shp,sign=-1.0).wkt}
        outcomes={}
        for label,params in cases.items():
            try:
                data=source_page(session,params)
                count=data.get("count")
                outcomes[label]={
                    "gbif_indexed_count":count,
                    "source_http_success":True,
                    "n_first_page_original_records_returned":len(data.get("results") or []),
                    "has_source_facets":bool(data.get("facets")),
                }
                print(json.dumps({"block":block,"prior_2026_source_indexed_count":prev,
                    "variant":label,"current_indexed_count":count,
                    "source_quality_flag":p["source_pre_outcome_both_period_recording_eligible"]
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                outcomes[label]={"source_http_success":False,
                     "error":f"{type(exc).__name__}: {str(exc)[:180]}"}
                report["errors"][block+"_"+label]=outcomes[label]["error"]
            time.sleep(.1)
        recent=outcomes["full_current_50km_facet"].get("gbif_indexed_count")
        report["regions"][block]={
            "country":flora["regions"][block]["country"],
            "source_cell_id":row["cell_id"],
            "source_cell_polygon_valid":shp.is_valid,
            "prior_full_facet_source_indexed_count":prev,
            "prior_species_key_categories":previous.get("n_all_observed_species_key_categories"),
            "prior_country_GRIIS_species_keys":previous.get("n_observed_GRIIS_introduced_plant_species_keys"),
            "prior_two_period_source_eligible":p["source_pre_outcome_both_period_recording_eligible"],
            "current_full_source_count":recent,
            "current_indexed_count_drift":recent-prev if isinstance(recent,int) and isinstance(prev,int) else None,
            "query_controls":outcomes,
        }
    report["required_source_reconciliation_succeeded"]=not report["errors"] and all(
        v["current_full_source_count"]==v["prior_full_facet_source_indexed_count"] for v in report["regions"].values()
    )
    return report


def main():
    if len(sys.argv)!=2:
        raise SystemExit("Usage: script <original_v36_flora_source_JSON>")
    flora=json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    out=run(flora)
    file=Path("results/audit/GBIF_old_species_facet_current_query_mismatch_v49b.json")
    file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"source_consistent":out["required_source_reconciliation_succeeded"],
        "prior_vs_current":{k:[x["prior_full_facet_source_indexed_count"],x["current_full_source_count"]]
                            for k,x in out["regions"].items()},
        "errors":out["errors"]},ensure_ascii=False),flush=True)
    if not out["required_source_reconciliation_succeeded"]:
        raise SystemExit("PRIOR 50km recorded-plant source count discrepancy: do not infer railway ecological effect")


if __name__=="__main__":
    main()
