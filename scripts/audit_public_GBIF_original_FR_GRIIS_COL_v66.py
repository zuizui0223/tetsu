#!/usr/bin/env python3
"""v6.6: actual public GBIF August2026 interpreted Plantae rows, France.

Source files are the EXACT 48 publicly tested in v6.5 (fixed BEFORE
country GRIIS flora membership or RAIL effect inspection). The original
200km France Burgundy 50km equal-area source frame and 2010-2025
periods remain unchanged.

Link 2023 GBIF GRIIS checklist accepted plant keys to 2026 Catalogue
of Life XR current source species IDs, via independently verified
official v6.3 strict accepted map; synonyms are separately reported.
No old numeric keys silently matched to COL XR IDs.

Original interpreted row-grain PRESENT/DETECTED, source GBIF ID,
dataset+occurrenceID and (where available) eventID+species
de-duplication is order-independent at coordinate precision bands:
when duplicate event rows have different stated uncertainty, retain
the best stated precision for the event rather than using row order.
Only privacy-safe summaries leave this script; temporary hashed IDs
are destroyed with SQLite scratch.

This public Parquet development preview is a bounded non-exhaustive
source subset of GBIF DWCA DOI 10.15468/dl.8yrbe7, NOT a newly issued
DOI nor a verified whole-country GBIF archive. No train dispersal,
actual freight, or ecological presence/absence inference.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
from hashlib import blake2b,sha256
import json,math,re,secrets,sqlite3,tempfile
from pathlib import Path

from pyproj import Transformer
from shapely import wkt
from shapely.ops import unary_union

from audit_gbif_equalarea_macroblocks_v29 import BLOCKS,GRID_CRS,make_grid
from ingest_GBIF_DWCA_original_events_v58 import (
    _cell_of_lonlat,source_cell_frame,resolve_country_species_key_set
)

BLOCK="fr_burgundy"
PREVIEW_ROOT="development/2026-08-01-DWCA/occurrence/p_taxon_a5/"
PUBLIC_S3_BUCKET="gbif-public-data"
ORIGINAL_UNDERLYING_DOI="10.15468/dl.8yrbe7"
MAX_SOURCE_PARQUET_OBJECTS=48
MAX_BOUNDING_ROWS=2_000_000
BATCH_SIZE=25000
MAX_SOURCE_COUNTRY_GRIIS_EVENTS=1_500_000
STRICT_UNCERTAINTY_M=5000
SECONDARY_UNCERTAINTY_M=25000
STRICT_SOURCE_MIN_TAXA=5
STRICT_SOURCE_MIN_EVENTS=20

FIELDS={
    "gbifid":"VARCHAR","datasetkey":"VARCHAR","countrycode":"VARCHAR",
    "year":"INTEGER","specieskey":"VARCHAR","occurrencestatus":"VARCHAR",
    "decimallongitude":"DOUBLE","decimallatitude":"DOUBLE",
    "coordinateuncertaintyinmeters":"DOUBLE","occurrenceid":"VARCHAR",
    "eventid":"VARCHAR","eventdate":"VARCHAR","taxonrank":"VARCHAR",
}


def _numeric(value):
    try:
        if value is None or str(value).strip()=="":return None
        x=float(value)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):
        return None


def source_manifest(manifest):
    if manifest.get("complete") is not True:
        raise ValueError("Actual original-GBIF v6.5 metadata source incomplete")
    if set(manifest.get("regions",{}))!={BLOCK}:
        raise ValueError("Source must be exactly original fixed Burgundy block")
    row=manifest["regions"][BLOCK]
    if row.get("genuine_GBIF_original_interpreted_rows_verified") is not True or row.get("source_errors"):
        raise ValueError("Previous Parquet source object audit invalid")
    if row.get("n_source_objects_actually_scanned")!=MAX_SOURCE_PARQUET_OBJECTS:
        raise ValueError("Original pre-outcome 48 object source frame changed")
    selected=row.get("source_object_results")
    if not isinstance(selected,list) or len(selected)!=MAX_SOURCE_PARQUET_OBJECTS:
        raise ValueError("Original source object list missing")
    if sum(r.get("original_interpreted_source_rows_in_200km_bounding_envelope",0)
           for r in selected)!=row.get("n_actual_interpreted_event_source_rows_in_selected_files_and_bbox"):
        raise ValueError("Previous independent source row sum mismatches object list")
    keys=[]
    for item in selected:
        key=item.get("source_key")
        if not isinstance(key,str) or not key.startswith(PREVIEW_ROOT) or (
                not re.fullmatch(r"[A-Za-z0-9_./=-]+",key)) or (
                not key.endswith(".parquet")):
            raise ValueError("Unapproved external Parquet source object key")
        keys.append(key)
    if len(set(keys))!=len(keys):
        raise ValueError("Duplicate source Parquet objects in pre-outcome manifest")
    return keys,row


def sql_source(db,keys):
    if len(keys)!=MAX_SOURCE_PARQUET_OBJECTS:
        raise ValueError("Only the original exactly 48 tested source objects can be queried")
    files="["+",".join("'s3://"+PUBLIC_S3_BUCKET+"/"+key+"'" for key in keys)+"]"
    source=f"read_parquet({files}, hive_partitioning=true, hive_types={{'a5_r2':UBIGINT}}, union_by_name=true)"
    schema=db.execute(f"DESCRIBE SELECT * FROM {source}").fetchall()
    names={str(row[0]).lower():str(row[0]) for row in schema}
    required=("gbifid","datasetkey","countrycode","year","specieskey",
              "occurrencestatus","decimallongitude","decimallatitude")
    if not all(k in names for k in required):
        raise ValueError("Original public GBIF Parquet source lacks required interpreted fields")
    def col(k):
        n=names.get(k)
        return '"'+n.replace('"','""')+'"' if n else "NULL"
    selected=", ".join(f"TRY_CAST({col(k)} AS {t}) AS {k}" for k,t in FIELDS.items())
    cells=[wkt.loads(r["gbif_geometry_wgs84"]) for r in make_grid(BLOCK)]
    footprint=unary_union(cells)
    west,south,east,north=footprint.bounds
    predicate=(f"{col('countrycode')}='FR' AND TRY_CAST({col('year')} AS INTEGER) BETWEEN 2010 AND 2025 "
               f"AND TRY_CAST({col('decimallongitude')} AS DOUBLE) BETWEEN {west:.9f} AND {east:.9f} "
               f"AND TRY_CAST({col('decimallatitude')} AS DOUBLE) BETWEEN {south:.9f} AND {north:.9f}")
    sql=f"SELECT {selected} FROM {source} WHERE {predicate}"
    return sql,list(FIELDS),{"source_original_bbox":[west,south,east,north],
                             "n_true_Parquet_files":len(keys),
                             "original_interpreted_fields_in_source_schema":sorted(names)}


def audit_rows(rows,columns,strict_set,strict_plus_synonyms,
               expected_late_bbox_count=None,max_rows=MAX_BOUNDING_ROWS):
    if not isinstance(strict_set,set) or not strict_set or not strict_set.issubset(strict_plus_synonyms):
        raise ValueError("GRIIS accepted COL key source groups invalid")
    if max_rows<1 or max_rows>MAX_BOUNDING_ROWS:
        raise ValueError("Out-of-contract original source row bound")
    frame=source_cell_frame(BLOCK)
    proj=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    privacy_salt=secrets.token_bytes(16)
    def hashed(*parts):
        return blake2b("\x1e".join(str(x) for x in parts).encode("utf-8","replace"),
                       key=privacy_salt,digest_size=16).hexdigest()
    quality=Counter()
    datasets=defaultdict(set)
    rejected=Counter()
    with tempfile.TemporaryDirectory(prefix="tetsu_fr_gbif_original_") as tmp:
        db=sqlite3.connect(str(Path(tmp)/"events.sqlite"))
        try:
            db.execute("CREATE TABLE unique_gbif(h TEXT PRIMARY KEY)")
            db.execute("CREATE TABLE source_occ(h TEXT PRIMARY KEY, event_hash TEXT NOT NULL)")
            db.execute("""CREATE TABLE concepts (
                h TEXT PRIMARY KEY, period TEXT, cell TEXT, species TEXT,
                strict INTEGER, strict5 INTEGER, strict25 INTEGER,
                extended INTEGER, extended5 INTEGER, extended25 INTEGER
            )""")
            for n,values in enumerate(rows,1):
                if n>max_rows:
                    raise ValueError("Source original public 200km bbox row budget exceeded, abort inference")
                row=dict(zip(columns,values))
                quality["n_original_public_bounding_envelope_rows_2010_to_2025"]+=1
                year=row.get("year")
                if not isinstance(year,int) or isinstance(year,bool) or year<2010 or year>2025:
                    quality["n_invalid_source_year"]+=1
                    continue
                period="2010_2017" if year<=2017 else "2018_2025"
                if period=="2018_2025":
                    quality["n_original_late_period_bounding_envelope_rows"]+=1
                longitude=_numeric(row.get("decimallongitude"))
                latitude=_numeric(row.get("decimallatitude"))
                if longitude is None or latitude is None:
                    quality["n_source_missing_coordinates"]+=1
                    continue
                cell=_cell_of_lonlat(BLOCK,longitude,latitude,proj,frame)
                if cell is None:
                    quality["n_source_rows_outside_exact_frozen_50km_cell_polygons"]+=1
                    continue
                quality["n_source_rows_inside_exact_50km_cells"]+=1
                status=str(row.get("occurrencestatus") or "").strip().upper()
                if status not in ("PRESENT","DETECTED"):
                    quality["n_source_non_PRESENT_or_unknown_status"]+=1
                    continue
                if str(row.get("taxonrank") or "").upper().strip()!="SPECIES":
                    quality["n_source_records_not_interpreted_species_rank"]+=1
                    continue
                taxon=str(row.get("specieskey") or "").strip()
                if taxon not in strict_plus_synonyms:
                    quality["n_source_speciesKey_not_in_country_GRIIS_accepted_unknown"]+=1
                    continue
                is_strict=taxon in strict_set
                quality["n_source_GRIIS_matched_interpreted_record_rows_before_dedup"]+=1
                if is_strict:quality["n_source_strict_country_GRIIS_matched_rows_before_dedup"]+=1
                original=str(row.get("gbifid") or "").strip()
                dataset=str(row.get("datasetkey") or "").strip()
                if not original or not dataset:
                    quality["n_source_missing_record_ID_or_dataset"]+=1
                    continue
                if db.execute("INSERT OR IGNORE INTO unique_gbif VALUES (?)",
                              (hashed("GBIF",original),)).rowcount==0:
                    quality["n_source_duplicate_GBIF_unique_ID"]+=1
                    continue
                occ=str(row.get("occurrenceid") or "").strip()
                evt=str(row.get("eventid") or "").strip()
                date=str(row.get("eventdate") or "").strip()
                # The repeated dataset+occurrence ID must refer to the
                # SAME biological source event/taxon/cell; otherwise no
                # ordering convention can safely choose a canonical row.
                # When eventID is absent, use a stable dataset/occurrence
                # key, not a read-order-selected GBIF row ID.
                unit=("EVENT",dataset,evt,date,taxon,cell,period) if evt else (
                    "DATASET_OCC",dataset,occ,taxon,cell,period) if occ else (
                    "GBIF_OCC",original,taxon,cell,period)
                identifier=hashed(*unit)
                if occ:
                    hashed_occ=hashed("DATASET_OCCURRENCE",dataset,occ)
                    prior_event=db.execute(
                        "SELECT event_hash FROM source_occ WHERE h=?",(hashed_occ,)
                    ).fetchone()
                    if prior_event is not None:
                        if prior_event[0]!=identifier:
                            raise ValueError("Conflicting event/taxon/cell under one dataset+occurrenceID, cannot resolve source by input row order")
                        quality["n_source_duplicate_same_dataset_occurrenceID"]+=1
                    else:
                        db.execute("INSERT INTO source_occ(h,event_hash) VALUES (?,?)",
                                   (hashed_occ,identifier))
                        quality["n_country_GRIIS_presence_records_after_occurrenceID_dedup"]+=1
                else:
                    quality["n_country_GRIIS_presence_records_after_occurrenceID_dedup"]+=1
                datasets[(period,cell)].add(dataset)
                uncertainty=_numeric(row.get("coordinateuncertaintyinmeters"))
                if uncertainty is None or uncertainty<=0:
                    quality["n_source_missing_or_nonpositive_coordinate_uncertainty"]+=1
                    has5=False;has25=False
                else:
                    has5=uncertainty<=STRICT_UNCERTAINTY_M
                    has25=uncertainty<=SECONDARY_UNCERTAINTY_M
                    if not has5:quality["n_source_known_uncertainty_gt5000m"]+=1
                    if not has25:quality["n_source_known_uncertainty_gt25000m"]+=1
                flags=(int(is_strict),int(is_strict and has5),
                       int(is_strict and has25),1,int(has5),int(has25))
                current=db.execute("SELECT strict,strict5,strict25,extended,extended5,extended25 FROM concepts WHERE h=?",
                                   (identifier,)).fetchone()
                if current is None:
                    db.execute("INSERT INTO concepts VALUES (?,?,?,?,?,?,?,?,?,?)",
                               (identifier,period,cell,taxon,*flags))
                else:
                    quality["n_source_duplicate_same_event_species_group"]+=1
                    improved=tuple(max(x,y) for x,y in zip(current,flags))
                    if improved!=current:
                        db.execute("""UPDATE concepts SET strict=?,strict5=?,strict25=?,
                            extended=?,extended5=?,extended25=? WHERE h=?""",
                            (*improved,identifier))
                if n%25000==0:db.commit()
            db.commit()
            if (expected_late_bbox_count is not None and
                quality["n_original_late_period_bounding_envelope_rows"]!=expected_late_bbox_count):
                raise ValueError(f"Actual original Parquet 2018-25 source row count disagrees with v6.5 pre-outcome file audit: {quality['n_original_late_period_bounding_envelope_rows']} vs {expected_late_bbox_count}")
            per={
                p:{f"{BLOCK}_x{i}_y{j}":{
                    "n_GRIIS_strict_COL_event_species":0,
                    "n_GRIIS_strict_COL_event_species_known_uncertainty_le5000m":0,
                    "n_GRIIS_strict_COL_event_species_known_uncertainty_le25000m":0,
                    "n_GRIIS_strict_plus_synonym_COL_event_species":0,
                    "n_GRIIS_strict_plus_synonym_event_species_known_uncertainty_le5000m":0,
                    "n_GRIIS_strict_plus_synonym_event_species_known_uncertainty_le25000m":0,
                    "n_distinct_source_dataset_UUIDs_with_any_GRIIS_event":len(datasets[(p,f"{BLOCK}_x{i}_y{j}")]),
                    "COL_strict_country_GRIIS_speciesKey_event_counts":{},
                    "COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts":{},
                    "COL_strict_plus_synonym_country_GRIIS_speciesKey_event_counts":{}
                } for j in range(4) for i in range(4)}
                for p in ("2010_2017","2018_2025")
            }
            cursor=db.execute("""SELECT period,cell,species,COUNT(*),
                SUM(strict),SUM(strict5),SUM(strict25),
                SUM(extended),SUM(extended5),SUM(extended25)
                FROM concepts GROUP BY period,cell,species""")
            for period,cell,taxon,n,strict,s5,s25,ext,e5,e25 in cursor:
                item=per[period][cell]
                item["n_GRIIS_strict_COL_event_species"]+=strict
                item["n_GRIIS_strict_COL_event_species_known_uncertainty_le5000m"]+=s5
                item["n_GRIIS_strict_COL_event_species_known_uncertainty_le25000m"]+=s25
                item["n_GRIIS_strict_plus_synonym_COL_event_species"]+=ext
                item["n_GRIIS_strict_plus_synonym_event_species_known_uncertainty_le5000m"]+=e5
                item["n_GRIIS_strict_plus_synonym_event_species_known_uncertainty_le25000m"]+=e25
                if strict:
                    item["COL_strict_country_GRIIS_speciesKey_event_counts"][taxon]=strict
                if s5:
                    item["COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts"][taxon]=s5
                if ext:
                    item["COL_strict_plus_synonym_country_GRIIS_speciesKey_event_counts"][taxon]=ext
            summaries={}
            for p,cells in per.items():
                summaries[p]={
                    "n_exact_source_cells_with_strict_GRIIS_event_species":sum(
                        bool(x["n_GRIIS_strict_COL_event_species"]) for x in cells.values()),
                    "n_strict_country_GRIIS_distinct_event_species_counts":sum(
                        x["n_GRIIS_strict_COL_event_species"] for x in cells.values()),
                    "n_strict_country_GRIIS_event_species_known_uncertainty_le5000m":sum(
                        x["n_GRIIS_strict_COL_event_species_known_uncertainty_le5000m"] for x in cells.values()),
                    "n_country_GRIIS_strict_plus_synonym_event_species":sum(
                        x["n_GRIIS_strict_plus_synonym_COL_event_species"] for x in cells.values()),
                    "n_source_cells_strict_GRIIS_at_least5_species_and20_event_species":sum(
                        len(x["COL_strict_country_GRIIS_speciesKey_event_counts"])>=STRICT_SOURCE_MIN_TAXA and
                        x["n_GRIIS_strict_COL_event_species"]>=STRICT_SOURCE_MIN_EVENTS for x in cells.values()),
                    "n_source_cells_strict_GRIIS_at_least5_species_and20_5km_precise_event_species":sum(
                        len(x["COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts"])>=STRICT_SOURCE_MIN_TAXA and
                        x["n_GRIIS_strict_COL_event_species_known_uncertainty_le5000m"]>=STRICT_SOURCE_MIN_EVENTS for x in cells.values()),
                    "cells":cells
                }
            return summaries,dict(quality)
        finally:
            db.close()


def run(preflight,bridge,checklists):
    source_keys,source_info=source_manifest(preflight)
    strict,taxmeta=resolve_country_species_key_set(checklists,BLOCK,bridge,False)
    with_syn,taxmeta2=resolve_country_species_key_set(checklists,BLOCK,bridge,True)
    if not strict.issubset(with_syn):
        raise ValueError("Official source COL accepted taxon sets inconsistent")
    import duckdb
    db=duckdb.connect()
    db.execute("SET threads=2; SET memory_limit='3GB'; SET enable_progress_bar=false")
    db.execute("INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1'")
    try:
        query,columns,info=sql_source(db,source_keys)
        cur=db.execute(query)
        def stream():
            while True:
                rows=cur.fetchmany(BATCH_SIZE)
                if not rows:break
                for r in rows:yield r
        periods,quality=audit_rows(stream(),columns,strict,with_syn,
                        expected_late_bbox_count=source_info["n_actual_interpreted_event_source_rows_in_selected_files_and_bbox"])
    finally:
        db.close()
    result={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source":"ACTUAL original interpreted GBIF August2026 A5 development preview Parquet event rows",
        "underlying_official_original_GBIF_DWCA_DOI":ORIGINAL_UNDERLYING_DOI,
        "this_is_NOT_a_new_GBIF_download_or_issued_DOI":True,
        "geographic_block":BLOCK,
        "source_country":"FR",
        "original_exact_16_equal_area_50km_source_cells":True,
        "n_pre_outcome_verified_source_Parquet_objects":len(source_keys),
        "n_unscanned_candidate_source_Parquet_objects":source_info["source_candidate_objects_in_centroid_A5_bound"]-len(source_keys),
        "n_positive_2018to2025_original_source_Parquet_objects":source_info["n_source_objects_with_nonzero_actual_source_occurrences"],
        "source_Parquet_object_keys_sha256":sha256("\n".join(source_keys).encode()).hexdigest(),
        "source_selection_verified_before_any_country_GRIIS_taxon_outcome":True,
        "source_original_2018to2025_bbox_rows_verified_against_v65":True,
        "GRIIS_col_taxonomy":taxmeta,
        "GRIIS_col_synonym_sensitivity_taxonomy":taxmeta2,
        "source_sql_metadata":info,
        "source_quality":quality,
        "periods":periods,
        "no_source_complete_census_or_actual_train_seed_transport_inferred":True,
        "limits":[
            "Only 48 of 71 source small Parquet files were scanned in ORIGINAL A5 centroid partitions, with more source objects and >=400 MB files excluded; missing published original plant records may change species/source sample.",
            "Source A5 r2 geographic selection used original cell centroids, not all overlapping tile polygons, so a 200km-region complete plant census is NOT established.",
            "Original subset may include multiple events from the same biological observer or survey even after record and event/species ID dedup; absence remains unknown.",
            "Official source GRIIS old GBIF numeric Backbone2023 -> May2026 COL XR accepted/synonym crosswalk is not guaranteed identical to the August2026 preview taxonomy snapshot; unmapped keys stay unknown.",
            "Source mapped 2026 railway/road physical graph is not an active freight exposure or a cause of plant presence; no dispersal rate, colonization or adaptation demonstrated.",
            "No raw GBIF occurrence IDs, human recorder names, eventIDs, precise event coordinates or user GBIF credentials stored in outputs.",
        ]
    }
    for period,vals in periods.items():
        print(json.dumps({
            "block":BLOCK,"period":period,
            "n_true_strict_COL_country_GRIIS_events":vals["n_strict_country_GRIIS_distinct_event_species_counts"],
            "n_known_coordinate_uncertainty_le5km_events":vals["n_strict_country_GRIIS_event_species_known_uncertainty_le5000m"],
            "n_including_synonym_map_events":vals["n_country_GRIIS_strict_plus_synonym_event_species"],
            "n_source_cells_GRIIS_5species_20events":vals["n_source_cells_strict_GRIIS_at_least5_species_and20_event_species"],
            "n_source_cells_precise5km_GRIIS_5species_20events":vals["n_source_cells_strict_GRIIS_at_least5_species_and20_5km_precise_event_species"],
            "source_nonexhaustive":True
        },ensure_ascii=False),flush=True)
    print(json.dumps({"complete":True,
        "FR_country_source_pre_outcome_selected_objects":len(source_keys),
        "all_plants_late_bbox_record_rows":quality["n_original_late_period_bounding_envelope_rows"],
        "validated_2018to2025_source_denominator":source_info["n_actual_interpreted_event_source_rows_in_selected_files_and_bbox"],
        "raw_true_source_50km_rows":quality["n_source_rows_inside_exact_50km_cells"],
        "COL_strict_GRIIS_source_not_train_dispersal":True},ensure_ascii=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pre-outcome-original-source-scan",required=True)
    p.add_argument("--sixcountry-COL-taxonomy-bridge",required=True)
    p.add_argument("--griis-original-country-keys",required=True)
    p.add_argument("--output",default="results/audit/public_GBIF_original_FR_GRIIS_COL_50km_event_quality_v66.json")
    args=p.parse_args()
    original=json.loads(Path(args.pre_outcome_original_source_scan).read_text(encoding="utf-8"))
    bridge=json.loads(Path(args.sixcountry_COL_taxonomy_bridge).read_text(encoding="utf-8"))
    keys=json.loads(Path(args.griis_original_country_keys).read_text(encoding="utf-8"))
    result=run(original,bridge,keys)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")


if __name__=="__main__":
    main()
