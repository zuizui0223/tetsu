#!/usr/bin/env python3
"""v2.0 spatial GBIF COMPLETE species-key facet matrix (source-side census).

Obtain *every* server-reported speciesKey facet, paging until exhaustion,
across fixed 3x3 geographical cells in ten preselected national windows
and two fixed periods. This returns indexed taxon counts per cell, not
a DOI-citable original GBIF occurrence-row download or ecological
colonization data. Missing species remain UNOBSERVED, not biological
absences. Source conservation: sequential requests and retry/backoff.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import time
import requests

from audit_overture_gbif_overlap_v19 import REGIONS,GBIF

PERIODS=((2010,2017),(2018,2025))
SPLIT_N=3
FACET_PAGE=200
MAX_FACET_PAGES=80
MIN_SOURCE_RECORDS_PER_CELL_PERIOD=100
MIN_DATASET_CATEGORIES=3


def grid_cells(region_name,divisions=SPLIT_N):
    if region_name not in REGIONS or divisions!=3:
        raise ValueError("Non-frozen source cell frame")
    country,w,s,e,n=REGIONS[region_name]
    dx=(e-w)/divisions
    dy=(n-s)/divisions
    items=[]
    for yi in range(divisions):
        for xi in range(divisions):
            x0=w+xi*dx
            x1=w+(xi+1)*dx
            y0=s+yi*dy
            y1=s+(yi+1)*dy
            if not(x1>x0 and y1>y0):
                raise ValueError("Malformed cell polygon")
            items.append({
                "cell_id":f"{region_name}_x{xi}_y{yi}",
                "country":country,"bbox":(round(x0,7),round(y0,7),round(x1,7),round(y1,7)),
            })
    return items


def polygon(bbox):
    w,s,e,n=bbox
    if not all(math.isfinite(x) for x in bbox) or not(e>w and n>s):
        raise ValueError("Bad spatial source region")
    return f"POLYGON(({w} {s},{e} {s},{e} {n},{w} {n},{w} {s}))"


def query_params(cell,period,facet_name,offset=0,limit=FACET_PAGE):
    if period not in PERIODS or facet_name not in ("speciesKey","datasetKey"):
        raise ValueError("Non-frozen source query")
    if not isinstance(offset,int) or offset<0 or limit<1 or limit>FACET_PAGE:
        raise ValueError("Unbounded facet request")
    return {
        "country":cell["country"],"geometry":polygon(cell["bbox"]),
        "kingdomKey":6,"taxonRank":"SPECIES",
        "year":f"{period[0]},{period[1]}",
        "occurrenceStatus":"PRESENT",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "limit":0,"facet":facet_name,
        "facetLimit":limit,"facetOffset":offset
    }


def parse_page(payload,facet_name):
    if not isinstance(payload,dict):
        raise ValueError("Malformed GBIF source response")
    n=payload.get("count")
    if isinstance(n,bool) or not isinstance(n,int) or n<0:
        raise ValueError("Invalid GBIF source record count")
    fs=payload.get("facets")
    if not isinstance(fs,list):
        raise ValueError("GBIF source facets missing")
    target=facet_name.lower()
    out=[f for f in fs if isinstance(f,dict) and
         "".join(ch for ch in str(f.get("field","")).lower() if ch.isalnum())==target]
    if n==0 and not out:
        return 0,[]
    if len(out)!=1 or not isinstance(out[0].get("counts"),list):
        raise ValueError(f"Expected exactly one source {facet_name} facet")
    rows=out[0]["counts"]
    if len(rows)>FACET_PAGE:
        raise ValueError("GBIF violated facet page size")
    data=[]
    for row in rows:
        if not isinstance(row,dict):
            raise ValueError("Non-object GBIF facet category")
        key=str(row.get("name") or "").strip()
        count=row.get("count")
        if (not key or isinstance(count,bool) or
            not isinstance(count,int) or count<0):
            raise ValueError("Invalid GBIF facet category or record count")
        data.append((key,count))
    if len(set(k for k,_ in data))!=len(data):
        raise ValueError("Duplicate categories in one facet page")
    return n,data


def source_json(client,params):
    for retry in range(5):
        try:
            resp=client.get(GBIF,params=params,timeout=(20,80))
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException,ValueError):
            if retry==4: raise
            time.sleep(min(16,2**retry))
    raise RuntimeError("GBIF source retry exhausted")


def full_species_facet(client,cell,period,request=source_json,sleep=time.sleep):
    seen={}
    total=None
    page_count=0
    for page in range(MAX_FACET_PAGES):
        params=query_params(cell,period,"speciesKey",offset=page*FACET_PAGE)
        payload=request(client,params)
        n,items=parse_page(payload,"speciesKey")
        if total is None:
            total=n
        elif n!=total:
            raise ValueError("Source count changed within facet pagination")
        for key,c in items:
            if key in seen:
                raise ValueError("Repeated facet speciesKey across pages: offset not honored")
            seen[key]=c
        page_count+=1
        if len(items)<FACET_PAGE:
            break
        sleep(.12)
    else:
        raise ValueError("GBIF facet maximum page budget exceeded; cannot claim complete species list")
    if sum(seen.values())>total:
        raise ValueError("Summed species facet counts exceed source occurrence count")
    return {
        "indexed_occurrence_count":total,
        "n_species_key_categories":len(seen),
        "species_key_occurrence_counts":seen,
        "facet_pages":page_count,
        "full_facet_pagination_completed":True,
        "not_original_occurrence_rows":True,
    }


def dataset_coverage(client,cell,period,request=source_json):
    params=query_params(cell,period,"datasetKey",offset=0,limit=4)
    count,facets=parse_page(request(client,params),"datasetKey")
    return {"source_occurrence_count":count,
            "n_source_datasets_at_least":len(facets),
            "n_dataset_categories_capped_at":4}


def compute_cell_period(client,cell,period,request=source_json,sleep=time.sleep):
    species=full_species_facet(client,cell,period,request=request,sleep=sleep)
    datasets=dataset_coverage(client,cell,period,request=request)
    if datasets["source_occurrence_count"]!=species["indexed_occurrence_count"]:
        raise ValueError("GBIF different source counts between species and dataset facets")
    good=(species["indexed_occurrence_count"]>=MIN_SOURCE_RECORDS_PER_CELL_PERIOD
          and datasets["n_source_datasets_at_least"]>=MIN_DATASET_CATEGORIES)
    return {
        **species,**datasets,
        "source_eligible_for_recording_effort_screen":good,
        "species_keys_are_observed_only_not_biological_absence":True,
    }


def comparison(rows_by_period):
    a,b=PERIODS
    sa=rows_by_period[a];sb=rows_by_period[b]
    if not(sa["source_eligible_for_recording_effort_screen"]
           and sb["source_eligible_for_recording_effort_screen"]):
        return {"eligible":False,"indexed_taxon_jaccard":None,
                "recording_coverage_not_biological_turnover":True}
    x=set(sa["species_key_occurrence_counts"])
    y=set(sb["species_key_occurrence_counts"])
    union=x|y
    jacc=len(x&y)/len(union) if union else None
    return {"eligible":True,"indexed_taxon_jaccard":(round(jacc,6) if jacc is not None else None),
            "n_shared_observed_species_keys":len(x&y),
            "n_union_observed_species_keys":len(union),
            "recording_coverage_not_biological_turnover":True}


def run(names):
    if not names or len(names)!=len(set(names)) or any(n not in REGIONS for n in names):
        raise ValueError("Unexpected regional study frame")
    out={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "source":"GBIF /v1/occurrence/search faceted geographic counts only",
        "source_type":"completely paginated speciesKey x cell x period indexed RECORD aggregates",
        "full_ten_country_frozen_frame":list(REGIONS),
        "requested_region_names":list(names),
        "grid":"3x3 non-equal-area WGS84 boxes per city window; ONLY within-city comparisons",
        "periods":[list(x) for x in PERIODS],
        "eligibility":{"min_records":MIN_SOURCE_RECORDS_PER_CELL_PERIOD,
                       "min_source_datasets":MIN_DATASET_CATEGORIES},
        "regions":{},"errors":{},
    }
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-gbif-macro-cell-taxon-coverage-v2.0; source facets only"})
    for name in names:
        results={}
        for cell in grid_cells(name):
            row={"country":cell["country"],"cell_id":cell["cell_id"],
                 "bbox":cell["bbox"],"by_period":{},"observation_jaccard":None}
            cache={}
            for period in PERIODS:
                key=f"{period[0]}_{period[1]}"
                try:
                    record=compute_cell_period(client,cell,period)
                    row["by_period"][key]=record
                    cache[period]=record
                except Exception as exc:
                    out["errors"][f"{cell['cell_id']}_{key}"]=f"{type(exc).__name__}: {str(exc)[:200]}"
                    row["by_period"][key]={"source_error":True}
                time.sleep(.15)
            if len(cache)==len(PERIODS):
                row["observation_jaccard"]=comparison(cache)
            results[cell["cell_id"]]=row
            print(json.dumps({
                "cell":cell["cell_id"],"country":cell["country"],
                "period_record_counts":{
                    k:v.get("indexed_occurrence_count") for k,v in row["by_period"].items()},
                "period_species_key_counts":{
                    k:v.get("n_species_key_categories") for k,v in row["by_period"].items()},
                "eligible_both":bool(row["observation_jaccard"] and row["observation_jaccard"]["eligible"]),
                "n_errors":len(out["errors"])
            },ensure_ascii=False),flush=True)
        out["regions"][name]={
            "country":REGIONS[name][0],
            "n_geo_cells":len(results),"cells":results,
            "n_two_period_source_eligible_cells":sum(
                bool(x["observation_jaccard"] and x["observation_jaccard"]["eligible"])
                for x in results.values()
            )
        }
    out["source_complete"]=(not out["errors"] and len(out["regions"])==len(names))
    out["n_cell_period_source_queries_succeeded"]=sum(
        sum(1 for v in cell["by_period"].values() if "indexed_occurrence_count" in v)
        for region in out["regions"].values() for cell in region["cells"].values()
    )
    out["do_not_infer"]=[
        "Observed species keys are indexed source records and not verified biological occupation or absences.",
        "Within-region year-period record differences can reflect changing citizen-science efforts and duplicated observations.",
        "This source data has no alien-native country classification, vascular-only screening, operator route traffic, freight or railway network exposure.",
        "A geographic source facet census is not a DOI-bearing full occurrence-row download and lacks record coordinate uncertainty and observer de-duplication.",
        "Grid is WGS84 3x3 per pilot city, not an equal-area national grid; do not pool as macro population sample.",
        "No plant–rail coefficient, seed dispersal, mowing, flowering or genetic selection effect inferred."
    ]
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--regions",nargs="+",choices=list(REGIONS),
                    default=["jp_sendai","de_leipzig"])
    ap.add_argument("--output",default="results/audit/gbif_species_3x3_jp_de_v20.json")
    a=ap.parse_args()
    report=run(a.regions)
    path=Path(a.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":report["source_complete"],
                      "n_regions":len(report["regions"]),
                      "n_cell_period_source_queries_succeeded":report["n_cell_period_source_queries_succeeded"],
                      "errors":report["errors"]},ensure_ascii=False),flush=True)
    if not report["source_complete"]:
        raise SystemExit("Incomplete GBIF speciesKey × spatial cell source audit; no ecological inference")


if __name__=="__main__":
    main()
