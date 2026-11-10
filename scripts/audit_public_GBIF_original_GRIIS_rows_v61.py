#!/usr/bin/env python3
"""v6.1 REAL occurrence-event row QA in one GBIF anonymous Parquet file.

The ONE SOURCE FILE was fixed in v6.0 metadata/result before inspecting
taxon outcomes:
 s3://gbif-public-data/development/2026-08-01-DWCA/occurrence/
 p_taxon_a5/taxonpartition=Plantae/
 a5_r2=1621295865853378560/data_1.parquet

This is NOT the entire GBIF snapshot, NOT all the original 16 grid
cells or an independent presence/absence survey. Nevertheless, this
reads real individual interpreted event rows (not GBIF search facets),
assigns geographic points to exact original frozen 2500km2 EPSG6933
cells, audits row-grain presence, source duplicates, coordinate
uncertainty and country GRIIS accepted species keys, and outputs only
aggregated privacy-safe SOURCE counts.

The source snapshot is a GBIF development preview derived from the
published August 2026 DWCA DOI 10.15468/dl.8yrbe7. No new DOI is
created. No train dispersal, active freight or causal flora effect.
"""
from __future__ import annotations

import argparse,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import re

from pyproj import Transformer
from shapely import wkt

from audit_gbif_equalarea_macroblocks_v29 import GRID_CRS, BLOCKS, make_grid
from ingest_GBIF_DWCA_original_events_v58 import _cell_of_lonlat,source_cell_frame
from analyze_griis_alien_rail_source_v25 import valid_griis

BLOCK="us_midwest"
SOURCE_FILE=(
    "s3://gbif-public-data/development/2026-08-01-DWCA/"
    "occurrence/p_taxon_a5/taxonpartition=Plantae/"
    "a5_r2=1621295865853378560/data_1.parquet"
)
UNDERLYING_SOURCE_DOI="10.15468/dl.8yrbe7"
MAX_ROWS=500000
BATCH_SIZE=25000
SOURCE_YEARS=(2010,2025)
STRICT_UNCERTAINTY_M=5000.


def _to_float(v):
    try:
        if v is None:return None
        f=float(v)
        return f if math.isfinite(f) else None
    except (ValueError,TypeError):return None


def _present(v):
    return str(v or "").strip().upper() in ("PRESENT","DETECTED")


def source_query(db,block=BLOCK):
    if block!=BLOCK:
        raise ValueError("One prespecified source Parquet pilot only")
    frame=source_cell_frame(block)
    polys=[wkt.loads(row["gbif_geometry_wgs84"]) for row in make_grid(block)]
    south=min(p.bounds[1] for p in polys);north=max(p.bounds[3] for p in polys)
    west=min(p.bounds[0] for p in polys);east=max(p.bounds[2] for p in polys)
    uri=SOURCE_FILE
    if "'" in uri or ".." in uri:
        raise ValueError("Untrusted GBIF source parquet file")
    f=f"read_parquet('{uri}',hive_partitioning=true,hive_types={{'a5_r2':UBIGINT}})"
    schema=db.execute(f"DESCRIBE SELECT * FROM {f}").fetchall()
    names={str(row[0]).lower():str(row[0]) for row in schema}
    mandatory=("gbifid","datasetkey","countrycode","year","specieskey",
               "occurrencestatus","decimallongitude","decimallatitude")
    if not all(x in names for x in mandatory):
        raise ValueError("Actual public GBIF interpreted Parquet source lacks mandatory event fields")
    def fld(name):
        if name not in names:return "NULL"
        return '"'+names[name].replace('"','""')+'"'
    fields={
        "gbifid":"VARCHAR","datasetkey":"VARCHAR","countrycode":"VARCHAR",
        "year":"INTEGER","specieskey":"VARCHAR","occurrencestatus":"VARCHAR",
        "decimallongitude":"DOUBLE","decimallatitude":"DOUBLE",
        "coordinateuncertaintyinmeters":"DOUBLE","occurrenceid":"VARCHAR",
        "eventid":"VARCHAR","eventdate":"VARCHAR","taxonrank":"VARCHAR",
    }
    selected=", ".join(f"TRY_CAST({fld(name)} AS {dtype}) AS {name}"
                       for name,dtype in fields.items())
    query=f"""SELECT {selected} FROM {f}
     WHERE {fld('countrycode')}='US'
       AND TRY_CAST({fld('year')} AS INTEGER) BETWEEN 2010 AND 2025
       AND TRY_CAST({fld('decimallongitude')} AS DOUBLE) BETWEEN {west:.9f} AND {east:.9f}
       AND TRY_CAST({fld('decimallatitude')} AS DOUBLE) BETWEEN {south:.9f} AND {north:.9f}"""
    return query,list(fields),frame,{"bbox_west":west,"bbox_east":east,
                                  "bbox_south":south,"bbox_north":north,
                                  "available_source_interpreted_fields":sorted(names)}


def audit_rows(rows,columns,accepted,block=BLOCK):
    if block!=BLOCK:raise ValueError("Unexpected source region")
    if not isinstance(accepted,set) or not accepted:
        raise ValueError("GRIIS accepted species keys required; nonmembers unknown")
    frame=source_cell_frame(block)
    to_proj=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    out={
        p:{cid:{
            "n_original_source_rows_GRIIS_PRESENCE":0,
            "n_unique_raw_GBIF_source_IDs":0,
            "n_distinct_source_dataset_UUIDs":0,
            "n_removed_duplicate_dataset_occurrenceID":0,
            "n_removed_duplicate_dataset_event_species":0,
            "n_missing_coordinate_uncertainty":0,
            "n_known_uncertainty_gt5000m":0,
            "n_source_GRIIS_keys":0,
            "country_GRIIS_specieskeys_original_row_counts":{},
            "country_GRIIS_specieskeys_precise_5000m_counts":{}
        } for cid in frame["wkt"]}
        for p in ("2010_2017","2018_2025")
    }
    counters=Counter()
    seen_ids=set();seen_occurrence=set();seen_event=set()
    seen_dataset=defaultdict(set)
    for record in rows:
        row=dict(zip(columns,record))
        counters["raw_source_rows_in_original_bounding_rectangle"]+=1
        if counters["raw_source_rows_in_original_bounding_rectangle"]>MAX_ROWS:
            raise ValueError("Bounded nonrepresentative original-source sample exceeds 500k records")
        lon=_to_float(row["decimallongitude"]);lat=_to_float(row["decimallatitude"])
        if lon is None or lat is None:continue
        cid=_cell_of_lonlat(block,lon,lat,to_proj,frame)
        if cid is None:
            counters["source_rows_outside_exact_original_50km_cells"]+=1
            continue
        counters["source_rows_inside_exact_original_50km_footprints"]+=1
        year=row["year"]
        period="2010_2017" if year<=2017 else "2018_2025"
        status=str(row["occurrencestatus"] or "").upper()
        if not _present(status):
            counters["non_PRESENCE_or_unknown_status_rows_not_used_as_plant_presence"]+=1
            continue
        if row["taxonrank"] is not None and str(row["taxonrank"]).upper()!="SPECIES":
            counters["not_specifically_species_rank_rows"]+=1
            continue
        sp=str(row["specieskey"] or "")
        counters["source_numeric_speciesKey_records"]+=bool(re.fullmatch(r"\d+",sp))
        counters["source_NONNUMERIC_speciesKey_records"]+=not bool(re.fullmatch(r"\d+",sp))
        if sp not in accepted:
            counters["speciesKey_not_in_country_GRIIS_status_UNKNOWN"]+=1
            continue
        dataset=str(row["datasetkey"] or "")
        if not dataset:
            counters["missing_dataset_UUID"]+=1
            continue
        gbifid=str(row["gbifid"] or "")
        if not gbifid or gbifid in seen_ids:
            counters["missing_or_duplicate_GBIF_unique_source_key"]+=1
            continue
        seen_ids.add(gbifid)
        summary=out[period][cid]
        occ=str(row["occurrenceid"] or "")
        if occ:
            k=(dataset,occ)
            if k in seen_occurrence:
                summary["n_removed_duplicate_dataset_occurrenceID"]+=1
                continue
            seen_occurrence.add(k)
        evt=str(row["eventid"] or "")
        if evt:
            # Real source event and species repeated within one location
            # is only one event-species category; NO fake observer key.
            k=(dataset,evt,str(row["eventdate"] or ""),sp,cid,period)
            if k in seen_event:
                summary["n_removed_duplicate_dataset_event_species"]+=1
                continue
            seen_event.add(k)
        summary["n_original_source_rows_GRIIS_PRESENCE"]+=1
        summary["n_unique_raw_GBIF_source_IDs"]+=1
        seen_dataset[(period,cid)].add(dataset)
        counts=summary["country_GRIIS_specieskeys_original_row_counts"]
        counts[sp]=counts.get(sp,0)+1
        unc=_to_float(row["coordinateuncertaintyinmeters"])
        if unc is None or unc<=0:
            summary["n_missing_coordinate_uncertainty"]+=1
        else:
            if unc>STRICT_UNCERTAINTY_M:
                summary["n_known_uncertainty_gt5000m"]+=1
            else:
                strict=summary["country_GRIIS_specieskeys_precise_5000m_counts"]
                strict[sp]=strict.get(sp,0)+1
    for (p,cid),setds in seen_dataset.items():
        out[p][cid]["n_distinct_source_dataset_UUIDs"]=len(setds)
    for p,data in out.items():
        for cid,c in data.items():
            c["n_source_GRIIS_keys"]=len(c["country_GRIIS_specieskeys_original_row_counts"])
    return out,dict(counters)


def run(griis_source):
    import duckdb
    if not griis_source.get("complete"):
        raise ValueError("Original country accepted GRIIS source list incomplete")
    accepted=valid_griis(griis_source).get("US")
    if not accepted:raise ValueError("US country GRIIS accepted Plantae key source unavailable")
    db=duckdb.connect()
    db.execute("INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1';")
    db.execute("SET threads=2; SET memory_limit='3GB';")
    try:
        query,names,frame,meta=source_query(db)
        cur=db.execute(query)
        def stream():
            while True:
                items=cur.fetchmany(BATCH_SIZE)
                if not items:break
                for record in items:yield record
        out,counters=audit_rows(stream(),names,accepted)
    finally:
        db.close()
    byperiod={}
    for period,data in out.items():
        byperiod[period]={
            "n_exact_50km_source_cells_with_introduced_events":sum(
                x["n_original_source_rows_GRIIS_PRESENCE"]>0 for x in data.values()),
            "n_country_GRIIS_intro_positive_event_rows_in_this_ONE_source_Parquet_object":
                sum(x["n_original_source_rows_GRIIS_PRESENCE"] for x in data.values()),
            "n_country_GRIIS_intro_strict_5000m_precision_event_rows":
                sum(sum(x["country_GRIIS_specieskeys_precise_5000m_counts"].values())
                    for x in data.values()),
            "n_source_cells_with_K20_original_row_events_for_recording_feasibility":sum(
                x["n_original_source_rows_GRIIS_PRESENCE"]>=20 and x["n_source_GRIIS_keys"]>=5
                for x in data.values()),
            "cells":data,
        }
    report={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source":"ONE prespecified GBIF August2026 public A5-r2 Plantae interpreted Parquet OBJECT",
        "fixed_original_file":SOURCE_FILE,
        "underlying_GBIF_source_DOI":UNDERLYING_SOURCE_DOI,
        "source_object_selected_before_species_GRIIS_outcome_inspection":True,
        "country":"US","regional_block":BLOCK,
        "exact_original_50km_equal_area_cells":16,
        "real_original_occurrence_rows_streamed_from_public_AWS":True,
        "source_schema":meta,
        "n_country_US_accepted_GRIIS_plant_specieskeys":len(accepted),
        "summary":counters,
        "periods":byperiod,
        "not_a_complete_geographic_plant_census":True,
        "not_a_NEW_authenticated_GBIF_download_or_DOI":True,
        "not_a_railway_seed_transport_result":True,
        "limits":[
            "This single source A5-r2 file/one Plantae taxonpartition is NOT all GBIF records or a random sample of all global/country indexed records.",
            "Exact 50km cell membership is verified geometrically but A5 centroids and source file selection do not guarantee completeness across every needed geographic partition.",
            "Species keys from August 2026 development snapshot may involve mixed GBIF/Catalogue of Life taxonomy; compare numeric join support before any country-GRIIS ecological analysis.",
            "Non-GRIIS speciesKey is UNKNOWN status, not native.",
            "GBIF record ID, eventID and dataset/occurrenceID dedup is source provenance QA, not proof of independently surveyed biological plants.",
            "Missing >5000m coordinate uncertainty and 2026 source map geometry prevent historical location or active railway freight identification.",
            "No train propagules, railway verge adaptation or independent country population effect inferred."
        ]
    }
    print(json.dumps({
        "complete":True,"region":BLOCK,
        "raw_sample_rows":counters.get("raw_source_rows_in_original_bounding_rectangle",0),
        "raw_inside_original_50km_grid":counters.get("source_rows_inside_exact_original_50km_footprints",0),
        "numeric_source_species_keys":counters.get("source_numeric_speciesKey_records",0),
        "nonnumeric_source_species_keys":counters.get("source_NONNUMERIC_speciesKey_records",0),
        "GRIIS_source_events_both_periods":sum(
            v["n_country_GRIIS_intro_positive_event_rows_in_this_ONE_source_Parquet_object"]
            for v in byperiod.values()),
        "periods":{p:{
            "n_GRIIS_source_event_rows":v["n_country_GRIIS_intro_positive_event_rows_in_this_ONE_source_Parquet_object"],
            "n_strict_precision_rows":v["n_country_GRIIS_intro_strict_5000m_precision_event_rows"],
            "cells_K20":v["n_source_cells_with_K20_original_row_events_for_recording_feasibility"]
        } for p,v in byperiod.items()},
        "ONE_FILE_NOT_FULL":True
    },ensure_ascii=False),flush=True)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--griis-keys",required=True)
    p.add_argument("--output",default="results/audit/GBIF_public_original_A5_GRIIS_US_row_QA_v61.json")
    args=p.parse_args()
    keys=json.loads(Path(args.griis_keys).read_text(encoding="utf-8"))
    result=run(keys)
    dst=Path(args.output)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")


if __name__=="__main__":
    main()
