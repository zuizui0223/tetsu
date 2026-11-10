#!/usr/bin/env python3
"""v4.8b presence-aware bounded GBIF ORIGINAL occurrence-row quality audit.

GBIF indexed OCCURRENCE_STATUS FACET currently uses DETECTED/
NOT_DETECTED, while GET occurrenceStatus=DETECTED is HTTP 400 and
GET occurrenceStatus=PRESENT returns zero (live v4.9b/c/d contract tests).
Therefore this source diagnostic deliberately does NOT filter status via
GET, but audits the exact same no-status source count by STATUS FACET,
then examines first <=300 original rows for metadata completeness.

This is a matched-source **quality audit**, not a DOI-bearing download,
not an unbiased sample of full record corpus, not a biological
vegetation survey, and not a train-mediated seed transport test.
"""
from __future__ import annotations

import argparse,json,time,math
from datetime import datetime,timezone
from pathlib import Path

import requests

from audit_gbif_original_row_quality_v48 import (
    frozen_cells,source_params,source_page,audit_source_page,
    BLOCKS_SIX,PERIODS,PAGE_LIMIT
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS


def params_without_broken_legacy_status(cell,period):
    params=source_params(cell,period,limit=PAGE_LIMIT)
    removed=params.pop("occurrenceStatus",None)
    if removed!="PRESENT":
        raise ValueError("Original GBIF legacy status search expectation changed")
    return params


def parse_status_facet(payload):
    n=payload.get("count")
    if isinstance(n,bool) or not isinstance(n,int) or n<0:
        raise ValueError("No verified GBIF source status count")
    facts=payload.get("facets")
    if not isinstance(facts,list):
        raise ValueError("Missing OCCURRENCE_STATUS facet in actual API")
    candidates=[x for x in facts if isinstance(x,dict) and
                "".join(ch for ch in str(x.get("field","")).lower() if ch.isalnum())
                =="occurrencestatus"]
    if len(candidates)!=1 or not isinstance(candidates[0].get("counts"),list):
        if n==0 and not candidates:
            return {"indexed_source_count":0,"statuses":{},"unclassified_count":0}
        raise ValueError("No single valid OCCURRENCE_STATUS facet")
    counts={}
    for item in candidates[0]["counts"]:
        key=str(item.get("name") or "").strip()
        val=item.get("count")
        if not key or key in counts or isinstance(val,bool) or not isinstance(val,int) or val<0:
            raise ValueError("Malformed/duplicate occurrence status category")
        counts[key]=val
    summation=sum(counts.values())
    if summation>n:
        raise ValueError("GBIF facet statuses exceed indexed record count")
    return {"indexed_source_count":n,
            "statuses":counts,
            "unclassified_count":n-summation}


def run(flora,query=source_page,sleep=time.sleep):
    if flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("Prior GBIF/GRIIS original 50km source was incomplete")
    if not set(BLOCKS_SIX).issubset(flora.get("regions",{})):
        raise ValueError("Fixed six-country source frame omitted")
    sess=requests.Session()
    sess.headers.update({"User-Agent":"tetsu-GBIF-indexed-DETECTED-plus-bounded-raw-occurrence-quality-v4.8b"})
    output={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source":"Official GBIF occurrence/search; status facet DETECTED, raw FIRST PAGE intentionally no status filter",
        "source_status_legacy_broken_contract_evidence_runs":[42793141570,42793539814,42793849132],
        "previous_GRIIS_country_source_run":42761140912,
        "original_fixed_country_blocks":list(BLOCKS_SIX),
        "frozen_sample_cell_offsets":[[1,1],[2,2]],
        "regions":{},"errors":{},
        "no_raw_personal_names_or_precise_geocoordinates_exported":True,
    }
    for block in BLOCKS_SIX:
        orig=flora["regions"][block]
        if orig.get("country")!=BLOCKS[block][0]:
            raise ValueError("GBIF source region country drift")
        cells={}
        for c in frozen_cells(block):
            cid=c["cell_id"];prior=orig["cells"].get(cid)
            periods={}
            for period in PERIODS:
                try:
                    no_status=params_without_broken_legacy_status(c,period)
                    facets={**no_status,"limit":0,"facet":"OCCURRENCE_STATUS",
                            "facetLimit":8}
                    status=parse_status_facet(query(sess,facets))
                    page=query(sess,no_status)
                    audit=audit_source_page(c,period,page,prior_flora=prior)
                    if audit["indexed_source_search_total"]!=status["indexed_source_count"]:
                        raise ValueError("Indexed count changed between status facet and actual occurrence page")
                    detected=status["statuses"].get("DETECTED",0)
                    notdetected=status["statuses"].get("NOT_DETECTED",0)
                    detected_share=detected/status["indexed_source_count"] if status["indexed_source_count"] else None
                    audit["full_source_status_categories"]=status["statuses"]
                    audit["full_source_unclassified_occurrence_status_count"]=status["unclassified_count"]
                    audit["DETECTED_indexed_positive_record_count"]=detected
                    audit["NOT_DETECTED_indexed_nonoccurrence_record_count"]=notdetected
                    audit["source_DETECTED_record_share"]=round(detected_share,7) if detected_share is not None else None
                    audit["original_row_first_page_cannot_individually_filter_DETECTED_via_GET"]=True
                    baseline=audit.get("prior_indexed_source_record_count")
                    audit["prior_DETECTED_indexed_source_reconciled_to_5pct"]=(
                        baseline is None or
                        abs(detected-baseline)<=max(5,baseline*.05)
                    )
                    if baseline is not None and baseline>0 and
                       (detected==0 or not audit["prior_DETECTED_indexed_source_reconciled_to_5pct"]):
                        raise ValueError(f"Prior indexed DETECTED source denominator drift: prior {baseline}, current {detected}")
                    periods[period]=audit
                except Exception as exc:
                    output["errors"][cid+"_"+period]=f"{type(exc).__name__}: {str(exc)[:240]}"
                    periods[period]={"source_error":True}
                sleep(.15)
            cells[cid]=periods
        ok=[p for rows in cells.values() for p in rows.values()
            if p.get("actual_original_record_rows_examined") is not None]
        output["regions"][block]={
            "country":BLOCKS[block][0],
            "n_fixed_source_50km_cells":len(cells),
            "n_valid_cell_period_original_status_source_bins":len(ok),
            "n_original_first_page_rows_examined":sum(p["actual_original_record_rows_examined"] for p in ok),
            "n_full_source_DETECTED_indexed_records_at_fixed_cells":sum(p["DETECTED_indexed_positive_record_count"] for p in ok),
            "n_full_source_NOT_DETECTED_indexed_records_at_fixed_cells":sum(p["NOT_DETECTED_indexed_nonoccurrence_record_count"] for p in ok),
            "n_original_rows_missing_coordinate_uncertainty":sum(p["n_missing_coordinateUncertaintyInMeters_first_page"] for p in ok),
            "n_original_rows_missing_recordedBy":sum(p["n_missing_recordedBy_first_page"] for p in ok),
            "n_first_pages_nonrandom_due_to_truncation":sum(p["first_page_is_nonrandom_if_more_records"] for p in ok),
            "cells":cells
        }
        result=output["regions"][block]
        print(json.dumps({
            "block":block,
            "country":result["country"],
            "source_bins":result["n_valid_cell_period_original_status_source_bins"],
            "bounded_original_record_rows":result["n_original_first_page_rows_examined"],
            "source_DETECTED_records":result["n_full_source_DETECTED_indexed_records_at_fixed_cells"],
            "source_NOT_DETECTED_records":result["n_full_source_NOT_DETECTED_indexed_records_at_fixed_cells"],
            "source_missing_uncertainty":result["n_original_rows_missing_coordinate_uncertainty"],
            "source_missing_recorder":result["n_original_rows_missing_recordedBy"],
            "nonrandom_first_pages":result["n_first_pages_nonrandom_due_to_truncation"],
            "errors":len(output["errors"])
        },ensure_ascii=False),flush=True)
    output["complete"]=not output["errors"] and len(output["regions"])==len(BLOCKS_SIX)
    output["summary"]={
        "n_independent_country_blocks":len(output["regions"]),
        "n_frozen_source_cell_period_slots":len(BLOCKS_SIX)*2*len(PERIODS),
        "n_facet_plus_original_row_source_bins_complete":sum(
            x["n_valid_cell_period_original_status_source_bins"] for x in output["regions"].values()),
        "n_firstpage_actual_source_occurrences_inspected":sum(
            x["n_original_first_page_rows_examined"] for x in output["regions"].values()),
        "n_facet_DETECTED_indexed_records_full_24_selected_slots":sum(
            x["n_full_source_DETECTED_indexed_records_at_fixed_cells"] for x in output["regions"].values()),
        "n_facet_NOT_DETECTED_indexed_records_full_24_selected_slots":sum(
            x["n_full_source_NOT_DETECTED_indexed_records_at_fixed_cells"] for x in output["regions"].values()),
        "n_firstpage_missing_recordedBy":sum(x["n_original_rows_missing_recordedBy"] for x in output["regions"].values()),
        "n_firstpage_missing_coordinate_uncertainty":sum(
            x["n_original_rows_missing_coordinate_uncertainty"] for x in output["regions"].values()),
        "n_firstpage_nonrandom_source_bins":sum(
            x["n_first_pages_nonrandom_due_to_truncation"] for x in output["regions"].values()),
    }
    output["critical_interpretation_limits"]=[
        "OCCURRENCE_STATUS faceted DETECTED and NOT_DETECTED counts are exhaustive for FIXED inspected source bins, but the raw first-page original records are NOT individually validated as DETECTED because GET status filters are incompatible with the current indexed enum.",
        "Only 2 prespecified 50km cells per country and two years-wide periods are inspected; this cannot establish source metadata quality for all 105 region-eligible cells or globally.",
        "First <=300 records of GET search are nonrandom where count>300; their missingness/concentration rates must NOT be generalized as unbiased population properties.",
        "Raw GBIF source recorder names, individual occurrence identifiers and precise geolocation are never saved or uploaded; only safe source aggregates are retained.",
        "Even DETECTED indexed records are not an independent plot survey or plant seed arrival at a railway.",
        "For a citable full event-level occurrence baseline, authenticated bulk GBIF download with DOI is still necessary.",
    ]
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--flora",required=True)
    p.add_argument("--output",default="results/audit/GBIF_detected_status_original_rows_6country_v48b.json")
    a=p.parse_args()
    flora=json.loads(Path(a.flora).read_text(encoding="utf-8"))
    res=run(flora)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":res["complete"],"summary":res["summary"],
                      "errors":res["errors"]},ensure_ascii=False),flush=True)
    if not res["complete"]:
        raise SystemExit("GBIF indexed detected source and bounded raw metadata did not reconcile")
if __name__=="__main__":
    main()
