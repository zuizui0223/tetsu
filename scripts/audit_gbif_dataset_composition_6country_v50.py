#!/usr/bin/env python3
"""v5.0: GBIF dataset-specific source recording concentration for six country blocks.

This is a REAL indexed source-composition control, not downloaded GBIF
original occurrence records. For every country-GRIIS-eligible 50km cell
in one preselected country-block and both historic periods, fetch full
datasetKey facet counts from official GBIF occurrence/search WITHOUT
broken legacy occurrenceStatus=PRESENT filter; fetch the matching
OCCURRENCE_STATUS facet on the same request. The same GBIF query is
used for source time, country, EPSG6933-reprojected cell geometry,
kingdomKey=6, taxonRank=SPECIES, coordinates and geospatial issue gates.

Dataset counts are all Plantae; GRIIS introduced taxa are NOT dataset-
specific in this source. The resulting controls are imperfect proxies
for source recorder effort. Country-area and positive-introduced flora
gates are reused from the ORIGINAL v3.6 independent source study, no
reranking of cities by the ecological response.

Output stores dataset UUID/counts only; no personal observer ID, raw
coordinates or original individual GBIF occurrence rows.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import time
import requests

from audit_gbif_equalarea_macroblocks_v29 import BLOCKS,make_grid
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX,PERIODS
from harvest_griis_macro_50km_facets_v31 import species_params
from harvest_gbif_cell_species_facets_v20 import source_json

GBIF_FACET_PAGE=5000
MAX_DATASET_FACET_PAGES=4
MAX_STATUS_UNCLASSIFIED_SHARE=0.02
MAX_STATUS_NOT_DETECTED_SHARE=0.02
MAX_SOURCE_INDEXED_DRIFT_SHARE=0.10
MIN_DATASET_COVERAGE=0.98
SOURCE_SLEEP=0.1
# Frozen before US Southeast RAIL/ROAD and GRIIS ecological source response.
ADDITIONAL_OUT_OF_ANALYSIS_REGION=("us_southeast",)
ALLOWED_SOURCE_BLOCKS=BLOCKS_SIX+ADDITIONAL_OUT_OF_ANALYSIS_REGION


def source_query(cell,period,offset=0):
    if period not in PERIODS or offset<0 or offset%GBIF_FACET_PAGE:
        raise ValueError("Nonfrozen country/time or invalid dataset source offset")
    year=tuple(int(x) for x in period.split("_"))
    source=species_params(cell,year,0)
    for field in ("occurrenceStatus","facet","facetOffset","facetLimit",
                  "speciesKey.facetLimit","speciesKey.facetOffset"):
        source.pop(field,None)
    source.update({
        "limit":0,
        "facet":["datasetKey","occurrenceStatus"],
        "facetLimit":GBIF_FACET_PAGE,
        "datasetKey.facetLimit":GBIF_FACET_PAGE,
        "datasetKey.facetOffset":offset,
        "occurrenceStatus.facetLimit":10,
    })
    return source


def _facet_rows(data,field):
    if not isinstance(data,dict):
        raise ValueError("GBIF occurrence source response is not a JSON object")
    fs=data.get("facets")
    if not isinstance(fs,list):
        raise ValueError("GBIF faceted source response has no facets")
    target="".join(x for x in field.lower() if x.isalnum())
    hits=[f for f in fs if isinstance(f,dict) and
          "".join(x for x in str(f.get("field","")).lower() if x.isalnum())==target]
    if len(hits)!=1 or not isinstance(hits[0].get("counts"),list):
        raise ValueError(f"GBIF facet field {field} missing/duplicated")
    rows=hits[0]["counts"]
    parsed=[]
    for record in rows:
        if not isinstance(record,dict):
            raise ValueError("Malformed GBIF source category")
        key=str(record.get("name") or "").strip()
        n=record.get("count")
        if not key or isinstance(n,bool) or not isinstance(n,int) or n<0:
            raise ValueError("Malformed category name/count in GBIF source")
        parsed.append((key,n))
    if len({key for key,_ in parsed})!=len(parsed):
        raise ValueError("Duplicate facet categories")
    return parsed


def full_dataset_profile(client,cell,period,prior,request=source_json,sleep=time.sleep):
    if not isinstance(prior,int) or prior<=0:
        raise ValueError("Prior original source GBIF plant-record denominator unavailable")
    datasets={}
    total=None
    status={}
    pages=0
    for page in range(MAX_DATASET_FACET_PAGES):
        params=source_query(cell,period,page*GBIF_FACET_PAGE)
        payload=request(client,params)
        count=payload.get("count")
        if isinstance(count,bool) or not isinstance(count,int) or count<0:
            raise ValueError("Invalid actual GBIF occurrence count")
        if total is None:
            total=count
            status=dict(_facet_rows(payload,"occurrenceStatus"))
        elif total!=count:
            raise ValueError("GBIF indexed total changed during dataset facet pagination")
        entries=_facet_rows(payload,"datasetKey")
        if len(entries)>GBIF_FACET_PAGE:
            raise ValueError("GBIF dataset facet page larger than source cap")
        for key,n in entries:
            if key in datasets:
                raise ValueError("Duplicate datasetKey across GBIF facet source pages")
            datasets[key]=n
        pages+=1
        if len(entries)<GBIF_FACET_PAGE:
            break
        sleep(SOURCE_SLEEP)
    else:
        raise ValueError("GBIF dataset facet page cap reached; source incomplete")
    counted=sum(datasets.values())
    if counted>total or sum(status.values())>total:
        raise ValueError("Faceted source category totals exceed indexed records")
    if total==0:
        raise ValueError("Previous GRIIS-eligible botanical source cell has zero current Plantae records")
    drift=abs(total-prior)/prior
    dataset_coverage=counted/total
    detected=status.get("DETECTED",0)
    not_detected=status.get("NOT_DETECTED",0)
    other_status=total-sum(status.values())
    bad_detection=(not_detected+other_status)/total
    # Reject unknown status rather than treating unclassified as plant
    # occurrence records. Do not invent dataset×DETECTED join counts:
    # datasetKey and occurrenceStatus are separate marginal facets.
    enough=(dataset_coverage>=MIN_DATASET_COVERAGE and
            bad_detection<=MAX_STATUS_NOT_DETECTED_SHARE and
            drift<=MAX_SOURCE_INDEXED_DRIFT_SHARE and
            len(datasets)>=3 and detected>0)
    probs={k:n/total for k,n in datasets.items()}
    hhi=sum(p*p for p in probs.values())
    other_count=total-counted
    # Missing dataset IDs are not shared between cells in comparisons.
    hhi_including_unknown=hhi+(other_count/total)**2
    sorted_rows=sorted(datasets.items(),key=lambda x:(-x[1],x[0]))
    return {
        "source_indexed_all_Plantae_records":total,
        "source_reference_v36_indexed_all_Plantae_records":prior,
        "relative_indexed_source_count_drift":round(drift,6),
        "dataset_key_record_counts":datasets,
        "n_distinct_source_dataset_keys":len(datasets),
        "n_facet_source_pages":pages,
        "source_dataset_facet_fully_paged":True,
        "sum_records_with_known_dataset_keys":counted,
        "missing_datasetKey_indexed_records":other_count,
        "dataset_source_coverage_fraction":round(dataset_coverage,7),
        "top_datasetKey_source_share":round(sorted_rows[0][1]/total,7) if sorted_rows else None,
        "top_datasetKey_source_id":sorted_rows[0][0] if sorted_rows else None,
        "top_three_datasets_record_share":round(sum(n for _,n in sorted_rows[:3])/total,7),
        "dataset_source_HHI_including_unattributed_bucket":round(hhi_including_unknown,7),
        "effective_datasets_by_indexed_record_inverse_HHI":round(1/hhi_including_unknown,4) if hhi_including_unknown>0 else None,
        "status_source_counts":status,
        "status_DETECTED_share":round(detected/total,7),
        "status_NOT_DETECTED_share":round(not_detected/total,7),
        "status_unclassified_source_count":other_status,
        "source_profile_valid_for_recording_bias_covariates":enough,
        "quality_gates":{"minimum_dataset_source_coverage":MIN_DATASET_COVERAGE,
                         "maximum_non_detected_or_unclassified_fraction":MAX_STATUS_NOT_DETECTED_SHARE,
                         "maximum_old_to_current_indexed_record_count_drift":MAX_SOURCE_INDEXED_DRIFT_SHARE},
        "raw_GBIF_occurrence_ID_recordedBy_coordinates_NOT_stored":True,
        "NOT_GRIIS_aliens_only_source_dataset_distribution":True,
    }


def run_one_block(block,flora,request=source_json,sleep=time.sleep):
    if flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("Original country-GRIIS source incomplete")
    if block not in ALLOWED_SOURCE_BLOCKS:
        raise ValueError("Nonpredeclared independent or frozen holdout country source block")
    row=flora["regions"][block]
    if row.get("country")!=BLOCKS[block][0]:
        raise ValueError("Source country identity mismatch")
    grid={r["cell_id"]:r for r in make_grid(block)}
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-macro-GBIF-dataset-facets-recording-bias-v5.0"})
    source={
        "utc":datetime.now(timezone.utc).isoformat(),
        "block":block,
        "country":BLOCKS[block][0],
        "source":"Official GBIF occurrence/search fully paged DATASET_KEY and OCCURRENCE_STATUS facets",
        "status":"GBIF facet DETECTED != legacy GET occurrenceStatus=PRESENT; the latter is NEVER used",
        "source_eligible_cells_pre_frozen_from_v36":True,
        "cells":{},"errors":{},
    }
    for cid in sorted(grid):
        original=row["cells"][cid]
        if original.get("cell_id")!=cid or original.get("country")!=BLOCKS[block][0]:
            raise ValueError("Original GRIIS source cell geographic frame changed")
        both_eligible=(original.get("source_pre_outcome_both_period_recording_eligible") is True
            and all(original.get("periods",{}).get(p,{}).get(
                "alien_observed_minimum_5_species_and_20_records") is True for p in PERIODS))
        item={"country":BLOCKS[block][0],
              "original_GRIIS_two_period_source_eligible":both_eligible,
              "periods":{}}
        if both_eligible:
            for period in PERIODS:
                old=original["periods"][period]
                try:
                    p=full_dataset_profile(session,grid[cid],period,
                        old["indexed_all_plantae_source_records"],
                        request=request,sleep=sleep)
                    item["periods"][period]=p
                except Exception as exc:
                    item["periods"][period]={"source_error":True}
                    source["errors"][cid+"_"+period]=f"{type(exc).__name__}: {str(exc)[:240]}"
                sleep(SOURCE_SLEEP)
        source["cells"][cid]=item
    good=[v for v in source["cells"].values() if v["original_GRIIS_two_period_source_eligible"]]
    source["n_original_country_GRIIS_source_eligible_cells"]=len(good)
    source["n_source_dataset_facet_complete_both_period_cells"]=sum(
        all(v["periods"].get(p,{}).get("source_dataset_facet_fully_paged") for p in PERIODS)
        for v in good)
    source["n_source_covariate_valid_both_period_cells"]=sum(
        all(v["periods"].get(p,{}).get("source_profile_valid_for_recording_bias_covariates") for p in PERIODS)
        for v in good)
    source["complete"]=not source["errors"] and source["n_source_dataset_facet_complete_both_period_cells"]==len(good)
    source["limits"]=[
        "Faceted datasetKey counts are all-Plantae marginal source data, not GRIIS-introduced-only nor independent observer/survey effort; dataset×species joint counts unavailable.",
        "OCCURRENCE_STATUS facet is a separate indexed marginal; dataset-key counts may include a small NOT_DETECTED fraction and cannot be corrected dataset-by-dataset without raw downloads.",
        "Record-based concentration is an imperfect sampling source proxy, NOT a causal adjustment for monitoring bias.",
        "No individual recordedBy, site coordinates, GBIF occurrence IDs, species-by-dataset or survey protocols stored.",
        "Data selection inherited original 50km source and GRIIS record opportunities, not railway association outcomes.",
    ]
    top=[v["periods"].get("2018_2025",{}).get("top_datasetKey_source_share")
         for v in good if v["periods"].get("2018_2025",{}).get("top_datasetKey_source_share") is not None]
    from statistics import median
    print(json.dumps({
        "block":block,"country":source["country"],
        "original_GRIIS_eligible_cells":len(good),
        "complete_dataset_facet_cells":source["n_source_dataset_facet_complete_both_period_cells"],
        "dataset_covariate_valid_cells":source["n_source_covariate_valid_both_period_cells"],
        "median_recent_top1_dataset_record_share":round(median(top),6) if top else None,
        "source_errors":len(source["errors"])
    },ensure_ascii=False),flush=True)
    return source


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--block",choices=ALLOWED_SOURCE_BLOCKS,required=True)
    p.add_argument("--flora",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    flora=json.loads(Path(a.flora).read_text(encoding="utf-8"))
    result=run_one_block(a.block,flora)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":result["complete"],"block":a.block,
                      "source_errors":result["errors"]},ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("GBIF source dataset facet audit incomplete; country biological inference prohibited")


if __name__=="__main__":
    main()
