#!/usr/bin/env python3
"""v7.7 outcome-blind US Pacific third FRA crossing macro block GBIF GRIIS
source opportunity, using pre-frozen v2.9 16 equal-area 50km cells.

Current 2026 GBIF search indexed DETECTED status facets (legacy
occurrenceStatus=PRESENT filter is inconsistent with live source);
speciesKey counts are recorded taxon categories NOT true presence/absence.
Independent Natural Earth country polygon >=80% and 2010–17/2018–25
source opportunity gate before any crossing train-richness correlation.
This is not a railway ecological effect or original GBIF event download.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from datetime import datetime,timezone
from hashlib import sha256
from statistics import median

import requests
from audit_gbif_equalarea_macroblocks_v29 import make_grid, BLOCKS,PERIODS
from audit_naturalearth_macro_country_v35 import (
    ZIP_URL,MAX_ZIP_BYTES,source_countries,_cell_geometry,COUNTRY_FRAC_MIN)
from audit_gbif_griis_dominant_dataset_leaveout_v56 import whole_source_species
from analyze_griis_alien_rail_source_v25 import valid_griis
from harvest_gbif_cell_species_facets_v20 import parse_page
from audit_gbif_detected_original_rows_v48b import parse_status_facet
from audit_gbif_macro_country_coverage_v18 import get_json

BLOCK="us_pacific"
COUNTRY="US"
MIN_ALL_PLANT_RECORDS=100
MIN_PUBLISHED_DATASETS=3
MIN_COUNTRY_ALIEN_SPECIES=5
MIN_COUNTRY_ALIEN_INDEXED_RECORDS=20
MAX_STATUS_UNKNOWN_RATE=.02
MAX_BETWEEN_QUERY_COUNT_DRIFT=.02
GBIF_URL="https://api.gbif.org/v1/occurrence/search"
MAX_DETAILED_SOURCE_QUERIES=32

def country_polygon_fraction(raw):
    poly=source_countries(raw)[COUNTRY]
    cells=_cell_geometry(BLOCK)
    if len(cells)!=16:
        raise ValueError("Changed predeclared US Pacific source 50km grid")
    return {k:round(min(1.,max(0.,p.intersection(poly).area/p.area)),6)
            for k,p in cells.items()}

def dataset_source_params(cell,period):
    if cell["country"]!=COUNTRY or period not in PERIODS:
        raise ValueError("Frozen Pacific macro country/year source mismatch")
    return {
       "country":COUNTRY, "geometry":cell["gbif_geometry_wgs84"],
       "kingdomKey":6,"taxonRank":"SPECIES",
       "year":str(period[0])+","+str(period[1]),
       "hasCoordinate":"true","hasGeospatialIssue":"false",
       "limit":0,"facet":["datasetKey","occurrenceStatus"],
       "facetLimit":20,"datasetKey.facetLimit":20,
       "occurrenceStatus.facetLimit":10
    }

def dataset_facets(client,cell,period):
    payload=get_json(client,GBIF_URL,dataset_source_params(cell,period))
    count,items=parse_page(payload,"datasetKey")
    status=parse_status_facet(payload)
    if count!=status["indexed_source_count"]:
        raise ValueError("GBIF dataset/status source count mismatch")
    if count==0:
        return {"indexed":0,"n_dataset_source_categories_at_least":0,
                "top_dataset_fraction":None,"status_not_detected_or_unknown":None}
    non=status["statuses"].get("NOT_DETECTED",0)+status["unclassified_count"]
    if non>count:raise ValueError("GBIF DETECTED status marginal exceeds total")
    return {
      "indexed":count,
      "n_dataset_source_categories_at_least":len(items),
      "top_dataset_fraction":round(max((n for key,n in items),default=0)/count,6),
      "status_not_detected_or_unknown":round(non/count,7)
    }

def source_cell(client,cell,period,keys):
    category=dataset_facets(client,cell,period)
    total=category["indexed"]
    if total==0:
        return {**category,"complete":True,"source_eligible":False,
                "observed_GRIIS_country_specieskey_records":{},
                "exclude_reason":"no_indexed_source_records_not_plant_absence"}
    # Fetch ALL species-key facet pages, never sample only the first 5000.
    full=whole_source_species(client,cell,
                             str(period[0])+"_"+str(period[1]),sleep=lambda _:None)
    n=full["source_indexed_all_Plantae_records"]
    if abs(n-total)/total>MAX_BETWEEN_QUERY_COUNT_DRIFT:
        raise ValueError("GBIF indexed count drift during independent facet pulls >2%")
    non=full["source_non_DETECTED_or_unknown_fraction"]
    if non is None or non>MAX_STATUS_UNKNOWN_RATE:
        reason="unknown_or_not_detected_source_gt2pct"
    elif category["n_dataset_source_categories_at_least"]<MIN_PUBLISHED_DATASETS:
        reason="less_than_3_dataset_UUID_sources"
    elif n<MIN_ALL_PLANT_RECORDS:
        reason="less_than_100_Plantae_indexed_records"
    else:
        reason=None
    source={key:num for key,num in full["speciesKey_count_by_all_taxa"].items()
            if key in keys and num>0}
    if reason is None and (len(source)<MIN_COUNTRY_ALIEN_SPECIES or
                           sum(source.values())<MIN_COUNTRY_ALIEN_INDEXED_RECORDS):
        reason="less_than_5_GRIIS_taxa_or_20_indexed_records"
    return {
       **category,
       "complete":True,"source_eligible":reason is None,
       "n_all_Plantae_indexed_source_records":n,
       "n_GRIIS_introduced_source_taxon_categories":len(source),
       "n_GRIIS_introduced_indexed_source_records":sum(source.values()),
       "observed_GRIIS_country_specieskey_records":source,
       "non_DETECTED_or_unknown_source_fraction":non,
       "source_full_speciesKey_facet_pagination":full["source_speciesKey_full_pagination"],
       "source_species_facet_pages":full["n_source_speciesKey_facet_pages"],
       "exclude_reason":reason,
       "NOT_independent_observer_survey_or_confirmed_plant_presence":True
    }

def run(session,griis,download=requests.get):
    if BLOCK not in BLOCKS or BLOCKS[BLOCK][0]!=COUNTRY:
        raise ValueError("Original frozen Pacific country geography missing")
    accepted=valid_griis(griis).get(COUNTRY)
    if not isinstance(accepted,set) or len(accepted)<100:
        raise ValueError("Original official country GRIIS accepted taxonomy source missing")
    r=download(ZIP_URL,timeout=(20,120))
    r.raise_for_status()
    raw=r.content
    if not 500000<len(raw)<MAX_ZIP_BYTES:
        raise ValueError("Natural Earth country polygon source unexpected size")
    fractions=country_polygon_fraction(raw)
    geom=make_grid(BLOCK)
    if set(fractions)!={x["cell_id"] for x in geom}:
        raise ValueError("Pacific fixed grid polygons differ from source cell identity")
    out={
       "complete":True,"date_utc":datetime.now(timezone.utc).isoformat(),
       "study":"Genuine 2026 GBIF current DETECTED speciesKey facet source for third frozen US Pacific FRA train-count macro region",
       "predefined_original_source_block":BLOCK,"country":COUNTRY,
       "n_original_50km_cells":16,
       "GBIF_status_source_is_DETECTED_facet_NOT_legacy_present_filter":True,
       "Natural_Earth_zip_sha256":sha256(raw).hexdigest(),
       "country_fraction_minimum":COUNTRY_FRAC_MIN,
       "original_GRIIS_country_taxa_source":len(accepted),
       "cells":{},"errors":{}}
    for cell in geom:
        cid=cell["cell_id"]
        valid=fractions[cid]>=COUNTRY_FRAC_MIN
        row={"country":COUNTRY,"country_fraction":fractions[cid],
             "country_80percent_quality":valid,
             "periods":{},
             "n_possible_original_periods":len(PERIODS)}
        if not valid:
            row["pre_outcome_source_both_period_eligible"]=False
            row["excluded_without_GBIF_query_reason"]="country_area_less_than_80pct_not_species_absence"
            out["cells"][cid]=row
            continue
        for period in PERIODS:
            name=str(period[0])+"_"+str(period[1])
            try:
                row["periods"][name]=source_cell(session,cell,period,accepted)
            except Exception as exc:
                out["errors"][cid+"_"+name]=type(exc).__name__+": "+str(exc)[:150]
                row["periods"][name]={"source_error":True}
        row["pre_outcome_source_both_period_eligible"]=all(
            row["periods"].get(str(a)+"_"+str(b),{}).get("source_eligible") is True
            for a,b in PERIODS)
        out["cells"][cid]=row
        print(json.dumps({"cell":cid,"country_fraction":fractions[cid],
                           "both_period_qualified":row["pre_outcome_source_both_period_eligible"],
                           "n_source_errors":len(out["errors"])},ensure_ascii=False),flush=True)
    out["complete"]=not bool(out["errors"]) and len(out["cells"])==16
    elig=[k for k,v in out["cells"].items() if v.get("pre_outcome_source_both_period_eligible")]
    out["summary"]={
       "n_50km_cells_country_fraction_80percent":sum(v["country_80percent_quality"] for v in out["cells"].values()),
       "n_50km_cells_with_two_period_indexed_flora_and_country_GRIIS_source":len(elig),
       "n_source_periods_queried":sum(len(v["periods"]) for v in out["cells"].values()),
       "n_source_errors":len(out["errors"]),
       "n_indexed_plant_taxon_categories_2018_2025_qualified":len(
           {t for k in elig for t in out["cells"][k]["periods"]["2018_2025"]["observed_GRIIS_country_specieskey_records"]})
    }
    out["limits"]=[
       "Current GBIF DETECTED status facets cannot assign NOT_DETECTED to each speciesKey; up to 2% unknown source status is tolerated, so this is NOT true species presence, occurrence-level validation or biological absence.",
       "Source US Pacific was selected BEFORE looking at the train-vs-plant coefficient, but is a prechosen geographical anchor not random ecoregion replication.",
       "A source dataset UUID is not an independent observer, site, standardized plot or plant individual.",
       "Two GBIF periods with CURRENT 2026 source API are not genuine detection-corrected historic time series.",
       "Current FRA 2026 crossing traffic is not actual freight tonne-km, historical freight route utilization, train-carried seeds or mechanical gust measurements.",
       "Country fraction is not ground-level vegetation. Source crossing intersections reflect road density.",
       "This completes only source-eligibility acquisition; ecological rail × flora effect remains unestimated until independent source and trait coverage gates are verified."
    ]
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--griis",required=True);ap.add_argument("--output",required=True)
    a=ap.parse_args()
    g=json.loads(Path(a.griis).read_text(encoding="utf-8"))
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-pacific-frozen-GRIIS-source-v77"})
    out=run(client,g)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],"errors":out["errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Pacific GBIF source hard error; no ecological inference")

if __name__=="__main__":main()
