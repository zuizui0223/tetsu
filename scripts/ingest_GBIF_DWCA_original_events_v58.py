#!/usr/bin/env python3
"""v5.8 streaming GBIF DOI-bearing DWCA original-event QA (PREPARED, NOT RUN).

Input: a REAL completed GBIF occurrence DOWNLOAD ZIP with DWCA
occurrence.txt, plus the corresponding public GBIF download metadata
JSON containing status=SUCCEEDED, DOI and download key. The archive
must have been downloaded by the user; this script does not request
auth, create a DOI, download a zip or claim a real observed result.

Reads GBIF INTERPRETED occurrence.txt in UTF-8 TSV streaming mode.
Separately audits row-level PRESENT/ABSENT (or new DETECTED/
NOT_DETECTED) status, original precise coordinate uncertainty,
country/taxon, and SOURCE dataset/occurrence/event duplicates.
GBIF source IDs are never output: only salted/in-session hash digests
are used in a temporary SQLite QA database and discarded at exit.
Only country-GRIIS-listed accepted plant species keys count as
introduced; nonmembership stays UNKNOWN, never native.

The source 2026 mapped railway network is NOT a treatment exposure
here: this is a prerequisite independent GBIF original-record quality
ledger for pre-frozen EPSG6933 50km equal-area cells.

This source parser is intentionally conservative and may exclude
actual records without required metadata. Its quality estimates are
valid only for the downloaded (and publicly DOI-identified) GBIF
archive, not all global true plant occurrences.
"""
from __future__ import annotations

import argparse
from collections import defaultdict,Counter
import csv
from datetime import datetime,timezone
from hashlib import sha256,blake2b
import io,json,math
from pathlib import Path
import re,sqlite3,tempfile
import zipfile
from urllib.parse import urlparse

from pyproj import Transformer
from shapely import wkt
from shapely.geometry import Point

from audit_gbif_equalarea_macroblocks_v29 import BLOCKS,PILOT,make_grid,GRID_CRS
from analyze_griis_alien_rail_source_v25 import valid_griis

DOWNLOAD_KEY_PATTERN=re.compile(r"^[0-9]{7,}-[0-9]{8,}$")
DOI_PATTERN=re.compile(r"^10\.15468/dl\.[A-Za-z0-9_-]+$")
MAX_UNCERTAINTY_PRECISE_M=5000.
MAX_UNCERTAINTY_SENSITIVITY_M=25000.
PROGRESS_N=100000


def metadata_contract(metadata):
    """Validate supplied complete GBIF download identity. Not network-verified."""
    if not isinstance(metadata,dict):
        raise ValueError("GBIF download metadata JSON required")
    status=str(metadata.get("status",""))
    doi=str(metadata.get("doi") or "")
    key=str(metadata.get("key") or "")
    if status!="SUCCEEDED" or not DOI_PATTERN.fullmatch(doi):
        raise ValueError("Official GBIF downloaded status SUCCEEDED and download DOI required")
    if not DOWNLOAD_KEY_PATTERN.fullmatch(key):
        raise ValueError("Original GBIF download request key invalid")
    link=metadata.get("downloadLink") or metadata.get("downloadlink")
    host=urlparse(link).hostname if isinstance(link,str) else None
    if not isinstance(link,str) or urlparse(link).scheme!="https" or (
            host!="gbif.org" and (host is None or not host.endswith(".gbif.org"))):
        raise ValueError("GBIF download metadata must contain genuine HTTPS gbif.org hostname")
    return {"source_download_key":key,"source_download_doi":doi,
            "original_download_status":status,"user_supplied_metadata_not_verified_against_live_GBIF_API":True}


def _digest(*items):
    data="\x1e".join(str(x) for x in items).encode("utf-8",errors="replace")
    return blake2b(data,digest_size=16).hexdigest()


def _int(x):
    try:
        if x is None or not str(x).strip():return None
        return int(str(x).strip())
    except (ValueError,TypeError):
        return None


def _num(x):
    try:
        if x is None or not str(x).strip():return None
        v=float(str(x).strip())
        return v if math.isfinite(v) else None
    except (ValueError,TypeError):
        return None


def _cell_of_lonlat(block,longitude,latitude,transformed,cells):
    if not (-180<=longitude<=180 and -90<=latitude<=90):
        return None
    x,y=transformed.transform(longitude,latitude)
    from math import floor
    xx=floor((x-cells["xmin"])/50000)
    yy=floor((y-cells["ymin"])/50000)
    if not 0<=xx<4 or not 0<=yy<4:
        return None
    cid=f"{block}_x{xx}_y{yy}"
    poly=cells["wkt"][cid]
    if not poly.covers(Point(longitude,latitude)):
        return None
    return cid


def source_cell_frame(block):
    if block not in PILOT:
        raise ValueError("Only ten original geographical frozen study blocks allowed")
    country,lon,lat=BLOCKS[block]
    to_eq=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    x,y=to_eq.transform(lon,lat)
    polygons={c["cell_id"]:wkt.loads(c["gbif_geometry_wgs84"]) for c in make_grid(block)}
    if len(polygons)!=16:
        raise ValueError("Original frozen 16-cell geographic frame changed")
    return {"xmin":x-100000.,"ymin":y-100000.,
            "wkt":polygons,"country":country}


def _one_event_bin():
    return {
        "n_GBIF_interpreted_raw_rows":0,
        "n_source_raw_PRESENCE_rows_geographically_valid":0,
        "n_source_distinct_GBIF_occurrence_ids":0,
        "n_records_rejected_duplicate_GBIF_ids":0,
        "n_records_rejected_duplicate_dataset_occurrenceID":0,
        "n_occurrenceID_missing":0,
        "n_eventID_missing":0,
        "n_recorder_name_missing":0,
        "n_coord_uncertainty_missing":0,
        "n_coord_uncertainty_gt5000m":0,
        "n_coord_uncertainty_gt25000m":0,
        "n_distinct_dataset_keys":0,
        "country_GRIIS_introduced_event_speciesKey_count_all_valid_geometry":{},
        "country_GRIIS_introduced_event_speciesKey_count_precise_5000m":{},
        "country_GRIIS_introduced_event_speciesKey_count_sensitivity_25000m":{},
        "n_event_species_raw_composites_deduplicated":0,
        "n_event_species_precise_composites_deduplicated":0,
        "n_event_species_sensitivity_composites_deduplicated":0,
    }


def resolve_country_species_key_set(checklists,block,col_taxonomy_bridge=None,
                                    include_synonym_resolution=False):
    """Explicitly declare taxonomy key namespace before source-row inference.

    An old GRIIS country list has numeric GBIF Backbone accepted IDs.
    The August 2026 public original DWCA often uses alpha-numeric COL
    XR species IDs. Never treat old IDs as valid new keys by string
    comparison, or derive zero species to mean botanical absence.
    """
    if block not in PILOT:
        raise ValueError("Unregistered 50km country source block")
    country=BLOCKS[block][0]
    original=valid_griis(checklists).get(country)
    if not original:
        raise ValueError("Missing original GRIIS accepted country plant species keys")
    if col_taxonomy_bridge is None:
        if include_synonym_resolution:
            raise ValueError("Synonym resolution requires verified COL crosswalk")
        return original,{
            "species_key_namespace":"legacy_GBIF_Backbone_numeric",
            "n_original_country_GRIIS_specieskeys":len(original),
            "COL_XR_mapping_used":False,
        }
    if not isinstance(col_taxonomy_bridge,dict) or col_taxonomy_bridge.get(
            "official_mapping_url")!="https://download.checklistbank.org/col/gbif/gbif-col265.tsv.gz":
        raise ValueError("COL XR source bridge provenance missing or unverified")
    row=col_taxonomy_bridge.get("countries",{}).get(country)
    if not isinstance(row,dict) or row.get("country")!=country:
        raise ValueError("Country-specific GRIIS -> COL source bridge absent")
    direct=row.get("old_GRIIS_to_COL_strict_accepted_map")
    syn=row.get("old_GRIIS_to_COL_synonym_acceptedID_map")
    if not isinstance(direct,dict) or not isinstance(syn,dict):
        raise ValueError("COL source taxon mapping lacks accepted/synonym ID groups")
    if not set(direct).issubset(original) or not set(syn).issubset(original):
        raise ValueError("COL source taxonomic mapping does not originate in this country's GRIIS checklist")
    if set(direct)&set(syn):
        raise ValueError("Old country GRIIS taxon cannot be both accepted and synonym")
    if any(not isinstance(x,str) or not x for x in (*direct.values(),*syn.values())):
        raise ValueError("Empty/invalid current COL XR accepted species ID")
    values=set(direct.values())
    if include_synonym_resolution:
        values.update(syn.values())
    if not values:
        raise ValueError("No current COL XR accepted plant keys for this country")
    return values,{
        "species_key_namespace":"COL26.5_XR_accepted_taxon_ID",
        "source_version":col_taxonomy_bridge.get("source_taxonomy_version"),
        "source_official_mapping_url":col_taxonomy_bridge["official_mapping_url"],
        "include_synonym_resolution":bool(include_synonym_resolution),
        "n_original_country_GRIIS_specieskeys":len(original),
        "n_resolved_country_COL_species_IDs":len(values),
        "n_original_GRIIS_keys_unmapped_or_rejected":row.get("n_old_keys_unmapped_or_no_usable_COL"),
        "COL_XR_mapping_used":True,
        "COL_source_bridge_is_May2026_not_exact_original_snapshot_version":True,
    }


def analyze_archive(path,block,checklists,metadata,source_metadata_verified_by_API=False,
                    col_taxonomy_bridge=None,include_synonym_resolution=False):
    if block not in PILOT:
        raise ValueError("Original source country-block frame invalid")
    identity=metadata_contract(metadata)
    keys,taxonomy=resolve_country_species_key_set(
        checklists,block,col_taxonomy_bridge,include_synonym_resolution
    )
    frame=source_cell_frame(block)
    to_eq=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    f=Path(path)
    if not f.is_file() or not zipfile.is_zipfile(f):
        raise ValueError("Actual user-provided GBIF DWCA ZIP not found/valid")
    result={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "country":frame["country"],
        "frozen_region":block,
        "country_GRIIS_specieskey_taxonomy_namespace":taxonomy,
        "GBIF_original_download_identity":identity,
        "actual_original_archive_file_read":True,
        "source_metadata_matches_live_GBIF_API":bool(source_metadata_verified_by_API),
        "source_archive_sha256":None,
        "periods":{},
        "source_quality_counters":Counter(),
        "limits":[
            "GBIF original DWCA archive was user-provided; supplied DOI/status metadata is not independently authenticated unless API crosscheck succeeded.",
            "Source PRESENT / DETECTED parsed at ROW grain and unknown statuses do not enter recorded plant presence.",
            "Missing/greater than 5km coordinate uncertainty is excluded from the strict precise 50km source cohort, not treated as location absence.",
            "EventID grouping is a conservative source event-species key proxy; no uniform botanical sampling design, repeat observers or true field plot absences established.",
            "The source records' country GRIIS acceptance does not identify railway carriage or active historic freight.",
            "The importer streams and avoids exporting personal recorder names/occurrenceIDs; hashed QA source identifiers are deleted with temporary SQLite scratch.",
        ]
    }
    for period in ("2010_2017","2018_2025"):
        result["periods"][period]={
            cid:_one_event_bin() for cid in sorted(frame["wkt"])
        }
    hash_bytes=sha256()
    with f.open("rb") as inp:
        while True:
            part=inp.read(2*1024*1024)
            if not part:break
            hash_bytes.update(part)
    result["source_archive_sha256"]=hash_bytes.hexdigest()
    with zipfile.ZipFile(f) as archive:
        if "occurrence.txt" not in archive.namelist():
            raise ValueError("GBIF DWCA ZIP lacks interpreted occurrence.txt (not SIMPLE_CSV)")
        if "meta.xml" not in archive.namelist():
            raise ValueError("GBIF DWCA lacks meta.xml source schema descriptor")
        # The published GBIF DWCA archive occurrence.txt has a UTF-8
        # header row and tab-delimited interpreted values. Do NOT parse
        # verbatim.txt for accepted speciesKey/coordinates.
        with tempfile.TemporaryDirectory(prefix="tetsu_dwca_qa_") as temp:
            conn=sqlite3.connect(str(Path(temp)/"dedup.sqlite"))
            try:
                conn.execute("CREATE TABLE ids (hash TEXT PRIMARY KEY)")
                conn.execute("CREATE TABLE dataset_occurrence (hash TEXT PRIMARY KEY)")
                conn.execute("CREATE TABLE event_species (hash TEXT PRIMARY KEY)")
                seen_dataset=defaultdict(set)
                with archive.open("occurrence.txt","r") as input_file:
                    wrapped=io.TextIOWrapper(input_file,encoding="utf-8-sig",errors="strict",newline="")
                    reader=csv.DictReader(wrapped,delimiter="\t",quotechar='"')
                    required={"gbifID","speciesKey","decimalLongitude","decimalLatitude",
                              "datasetKey","occurrenceStatus","year","countryCode"}
                    if not reader.fieldnames or not required.issubset(reader.fieldnames):
                        raise ValueError("GBIF interpreted DWCA schema missing mandatory source columns")
                    for ix,row in enumerate(reader,1):
                        result["source_quality_counters"]["n_archive_rows_processed"]+=1
                        if not isinstance(row,dict) or None in row:
                            raise ValueError("Malformed GBIF DWCA source TSV row")
                        if str(row.get("countryCode") or "").strip()!=frame["country"]:
                            result["source_quality_counters"]["n_nonmatching_source_country"]+=1
                            continue
                        year=_int(row.get("year"))
                        if year is None or year<2010 or year>2025:
                            result["source_quality_counters"]["n_invalid_or_outside_source_year"]+=1
                            continue
                        period="2010_2017" if year<=2017 else "2018_2025"
                        lon=_num(row.get("decimalLongitude"))
                        lat=_num(row.get("decimalLatitude"))
                        if lon is None or lat is None:
                            result["source_quality_counters"]["n_missing_source_coordinates"]+=1
                            continue
                        cid=_cell_of_lonlat(block,lon,lat,to_eq,frame)
                        if cid is None:
                            result["source_quality_counters"]["n_outside_original_50km_WKT"]+=1
                            continue
                        entry=result["periods"][period][cid]
                        entry["n_GBIF_interpreted_raw_rows"]+=1
                        gbif_id=str(row.get("gbifID") or "").strip()
                        if not gbif_id:
                            result["source_quality_counters"]["n_missing_GBIF_record_primary_id"]+=1
                            continue
                        h=_digest("GBIF",gbif_id)
                        if conn.execute("INSERT OR IGNORE INTO ids (hash) VALUES (?)",(h,)).rowcount==0:
                            entry["n_records_rejected_duplicate_GBIF_ids"]+=1
                            continue
                        dataset=str(row.get("datasetKey") or "").strip()
                        if not dataset:
                            result["source_quality_counters"]["n_missing_dataset_uuid"]+=1
                            continue
                        occurrence_id=str(row.get("occurrenceID") or "").strip()
                        if occurrence_id:
                            hc=_digest(dataset,occurrence_id)
                            if conn.execute("INSERT OR IGNORE INTO dataset_occurrence (hash) VALUES (?)",(hc,)).rowcount==0:
                                entry["n_records_rejected_duplicate_dataset_occurrenceID"]+=1
                                continue
                        else:
                            entry["n_occurrenceID_missing"]+=1
                        status=str(row.get("occurrenceStatus") or "").strip().upper()
                        if status not in ("PRESENT","DETECTED"):
                            result["source_quality_counters"]["n_nonpresence_or_unknown_source_status"]+=1
                            continue
                        if str(row.get("taxonRank") or "").strip().upper()!="SPECIES":
                            result["source_quality_counters"]["n_not_species_taxon_rank"]+=1
                            continue
                        species=str(row.get("speciesKey") or "").strip()
                        # A non-numeric source species ID in a legacy
                        # numeric GRIIS query indicates taxonomy-version
                        # incompatibility, not verified plant absence.
                        if species and not species.isdecimal():
                            result["source_quality_counters"]["n_original_positive_source_rows_with_alphanumeric_COL_like_speciesKey"]+=1
                        if species not in keys:
                            result["source_quality_counters"]["n_country_GRIIS_nonmember_unknown_status"]+=1
                            continue
                        entry["n_source_raw_PRESENCE_rows_geographically_valid"]+=1
                        entry["n_source_distinct_GBIF_occurrence_ids"]+=1
                        seen_dataset[(period,cid)].add(dataset)
                        uncertainty=_num(row.get("coordinateUncertaintyInMeters"))
                        if uncertainty is None or uncertainty<=0:
                            entry["n_coord_uncertainty_missing"]+=1
                        else:
                            if uncertainty>MAX_UNCERTAINTY_PRECISE_M:
                                entry["n_coord_uncertainty_gt5000m"]+=1
                            if uncertainty>MAX_UNCERTAINTY_SENSITIVITY_M:
                                entry["n_coord_uncertainty_gt25000m"]+=1
                        if not str(row.get("recordedBy") or "").strip():
                            entry["n_recorder_name_missing"]+=1
                        event_id=str(row.get("eventID") or "").strip()
                        eventdate=str(row.get("eventDate") or "").strip()
                        if not event_id:
                            entry["n_eventID_missing"]+=1
                        # Same dataset event/species record counted once
                        # when a real eventID exists. Fallback to the
                        # unique GBIF occurrence ID when absent; no
                        # observer/coordinate pseudo-survey invented.
                        unit=("EVENT",dataset,event_id,eventdate,species,cid) if event_id else (
                            "OCCURRENCE",gbif_id,species,cid)
                        he=_digest(*unit)
                        if conn.execute("INSERT OR IGNORE INTO event_species(hash) VALUES (?)",(he,)).rowcount==0:
                            entry["n_event_species_raw_composites_deduplicated"]+=1
                            continue
                        count=entry["country_GRIIS_introduced_event_speciesKey_count_all_valid_geometry"]
                        count[species]=count.get(species,0)+1
                        if uncertainty is not None and uncertainty>0 and uncertainty<=MAX_UNCERTAINTY_SENSITIVITY_M:
                            q=entry["country_GRIIS_introduced_event_speciesKey_count_sensitivity_25000m"]
                            q[species]=q.get(species,0)+1
                        if uncertainty is not None and uncertainty>0 and uncertainty<=MAX_UNCERTAINTY_PRECISE_M:
                            q=entry["country_GRIIS_introduced_event_speciesKey_count_precise_5000m"]
                            q[species]=q.get(species,0)+1
                        if ix%PROGRESS_N==0:
                            conn.commit()
                conn.commit()
                for (period,cid),dset in seen_dataset.items():
                    result["periods"][period][cid]["n_distinct_dataset_keys"]=len(dset)
            finally:
                conn.close()
    if (not taxonomy["COL_XR_mapping_used"] and
        result["source_quality_counters"].get(
            "n_original_positive_source_rows_with_alphanumeric_COL_like_speciesKey",0)>0):
        raise ValueError("COL XR alpha-numeric speciesKey found in original GBIF source, but only old numeric GRIIS keys were supplied. Provide the verified country COL-XR bridge; zero matches is NOT ecological absence.")
    result["source_quality_counters"]=dict(result["source_quality_counters"])
    result["n_50km_region_source_cells"]=16
    result["n_year_periods"]=2
    result["n_original_GBIF_DWCA_rows_processed"]=result["source_quality_counters"].get("n_archive_rows_processed",0)
    result["n_analyzed_country_GRIIS_presence_records"]=sum(
        r["n_source_raw_PRESENCE_rows_geographically_valid"]
        for period in result["periods"].values() for r in period.values())
    result["n_analyzed_source_events_with_GRIIS_species_and_valid_coordinates"]=sum(
        sum(r["country_GRIIS_introduced_event_speciesKey_count_all_valid_geometry"].values())
        for period in result["periods"].values() for r in period.values())
    result["DOI_archived_event_data_actually_read"]=True
    result["no_train_railway_plant_causality_identified"]=True
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--block",required=True,choices=PILOT)
    p.add_argument("--dwca-zip",required=True)
    p.add_argument("--gbif-download-metadata",required=True)
    p.add_argument("--griis-accepted-keys",required=True)
    p.add_argument("--COL-XR-bridge",default=None,
        help="Verified v6.3 official six-country accepted/synonym COL XR mapping JSON. Required when source speciesKey is alphanumeric.")
    p.add_argument("--include-COL-synonyms",action="store_true",
        help="Secondary sensitivity: include accepted IDs from old GRIIS synonym matches, not primary strict mapping.")
    p.add_argument("--output",default="results/audit/source_DWCA_event_quality_v58.json")
    args=p.parse_args()
    metadata=json.loads(Path(args.gbif_download_metadata).read_text(encoding="utf-8"))
    checklists=json.loads(Path(args.griis_accepted_keys).read_text(encoding="utf-8"))
    bridge=json.loads(Path(args.COL_XR_bridge).read_text(encoding="utf-8")) if args.COL_XR_bridge else None
    result=analyze_archive(
        args.dwca_zip,args.block,checklists,metadata,
        col_taxonomy_bridge=bridge,
        include_synonym_resolution=args.include_COL_synonyms
    )
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "original_archive_read":True,
        "source_region":args.block,
        "n_GBIF_original_rows":result["n_original_GBIF_DWCA_rows_processed"],
        "n_country_GRIIS_source_presence_rows":result["n_analyzed_country_GRIIS_presence_records"],
        "n_country_GRIIS_event_species_records":result["n_analyzed_source_events_with_GRIIS_species_and_valid_coordinates"],
        "metadata_not_authenticated_against_live_GBIF_API":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
