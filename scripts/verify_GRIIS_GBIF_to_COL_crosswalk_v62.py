#!/usr/bin/env python3
"""v6.2 OFFICIAL legacy GBIF backbone -> COL26.5 XR GRIIS source-key crosswalk.

Fix a critical taxonomy compatibility bug BEFORE interpreting the actual
public August2026 GBIF Parquet original-event study as biological absence.
GBIF old Backbone 2023 speciesKey numeric source != COL Extended Release
alphanumeric keys in the public August 2026 interpreted DWCA preview.

Official precomputed crosswalk:
 https://download.checklistbank.org/col/gbif/gbif-col265.tsv.gz
 https://download.checklistbank.org/col/gbif/README.html
 May 2026 COL 26.5 XR; not guaranteed same taxonomy revision as the
 Aug 2026 source or current Nov 2026 COL XR. Preserve all failures.

The strict PRIMARY mapping accepts ONLY:
 - old GRIIS country-listed GBIF backbone accepted *species* key,
 - one unambiguous corresponding COL species status ACCEPTED or
   PROVISIONALLY ACCEPTED and verified kingdom PLANTAE.
 A separate SENSITIVITY maps COL synonyms via acceptedID; never
 silently treat synonyms as multiple independent species or impute
 unmapped country-listed taxa as absent.

Read the SAME authentic one-file US Midwest August 2026 original
Parquet source already fixed before outcome inspection in v6.0/v6.1.
Report country GRIIS row-level PRESENCE, same event/occurrence ID
dedup and coordinate precision using v6.1 original 50km source cells.
The source Parquet is NOT the full geographic/date/taxon population,
even when the crosswalk succeeds. No credentials or new DOI.
"""
from __future__ import annotations

import argparse,csv,gzip,io,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import requests

from analyze_griis_alien_rail_source_v25 import valid_griis
from audit_public_GBIF_original_GRIIS_rows_v61 import (
    BLOCK,SOURCE_FILE,UNDERLYING_SOURCE_DOI,BATCH_SIZE,
    audit_rows,source_query
)

OFFICIAL_CROSSWALK="https://download.checklistbank.org/col/gbif/gbif-col265.tsv.gz"
OFFICIAL_README="https://download.checklistbank.org/col/gbif/README.html"
COL_VERSION="26.5 XR / May 2026"
MAX_MAPPING_SOURCE_ROWS=15_000_000
STRICT_ACCEPTED=("ACCEPTED","PROVISIONALLY ACCEPTED")
REQUIRED_COLUMNS=(
    "gbif:ID","gbif:rank","gbif:status",
    "col:ID","col:rank","col:scientificName","col:status",
    "col:acceptedID","col:kingdom"
)


def map_country_griis_old_keys(row_iter,old_keys):
    """Read official GBIF-COL crosswalk; reject unmapped, non-species and non-Plantae.

    Never count one underlying COL species twice because two GBIF
    backbone names collapse to the same accepted taxon.
    """
    if not isinstance(old_keys,set) or not old_keys:
        raise ValueError("Old source country-GRIIS list required")
    result={
        "source_legacy_griis_accepted_backbone_plant_species_count":len(old_keys),
        "n_mapping_rows_scanned":0,
        "n_old_keys_found_in_public_mapping":0,
        "n_old_keys_unmapped_or_no_usable_COL":0,
        "n_old_backbone_species_mapping_to_accepted_COL_species":0,
        "n_old_backbone_species_mapping_to_COL_synonyms":0,
        "n_rejected_non_Plantae_or_non_species_or_missing":0,
        "n_ambiguous_duplicate_old_backbone_keys":0,
        "n_distinct_strict_COL_speciesIDs":0,
        "n_distinct_COL_speciesIDs_with_accepted_synonym_resolution":0,
        "n_many_to_one_old_species_collapsed_into_same_STRICT_COL_ID":0,
        "strict_accepted":set(),
        "synonym_accepted_additional":set(),
        "source_old_to_strict_COL":{},
        "source_old_to_synonym_resolution":{},
    }
    seen_old=set()
    for row in row_iter:
        result["n_mapping_rows_scanned"]+=1
        if result["n_mapping_rows_scanned"]>MAX_MAPPING_SOURCE_ROWS:
            raise ValueError("Official mapping row count exceeds safe source limit")
        if not isinstance(row,dict):
            raise ValueError("Malformed official mapping row")
        key=str(row.get("gbif:ID") or "").strip()
        if key not in old_keys:
            continue
        if key in seen_old:
            result["n_ambiguous_duplicate_old_backbone_keys"]+=1
            raise ValueError("Official mapping duplicate legacy source key")
        seen_old.add(key)
        oldrank=str(row.get("gbif:rank") or "").upper().strip()
        kingdom=str(row.get("col:kingdom") or "").upper().strip()
        colrank=str(row.get("col:rank") or "").upper().strip()
        status=str(row.get("col:status") or "").upper().strip()
        colid=str(row.get("col:ID") or "").strip()
        accepted=str(row.get("col:acceptedID") or "").strip()
        if oldrank!="SPECIES" or colrank!="SPECIES" or kingdom!="PLANTAE":
            result["n_rejected_non_Plantae_or_non_species_or_missing"]+=1
            continue
        if status in STRICT_ACCEPTED and colid:
            result["strict_accepted"].add(colid)
            result["source_old_to_strict_COL"][key]=colid
            result["n_old_backbone_species_mapping_to_accepted_COL_species"]+=1
        elif status=="SYNONYM" and accepted:
            result["synonym_accepted_additional"].add(accepted)
            result["source_old_to_synonym_resolution"][key]=accepted
            result["n_old_backbone_species_mapping_to_COL_synonyms"]+=1
        else:
            result["n_rejected_non_Plantae_or_non_species_or_missing"]+=1
    result["n_old_keys_found_in_public_mapping"]=len(seen_old)
    result["n_old_keys_unmapped_or_no_usable_COL"]=len(old_keys)-len(
        result["source_old_to_strict_COL"])-len(result["source_old_to_synonym_resolution"])
    result["n_distinct_strict_COL_speciesIDs"]=len(result["strict_accepted"])
    result["n_distinct_COL_speciesIDs_with_accepted_synonym_resolution"]=len(
        result["strict_accepted"]|result["synonym_accepted_additional"])
    result["n_many_to_one_old_species_collapsed_into_same_STRICT_COL_ID"]=(
        result["n_old_backbone_species_mapping_to_accepted_COL_species"]-
        result["n_distinct_strict_COL_speciesIDs"])
    if result["n_old_keys_unmapped_or_no_usable_COL"]<0:
        raise ValueError("Crosswalk mapping results have impossible negative missing source")
    return result


def fetch_official_mapping(country_old):
    if not country_old:
        raise ValueError("No old GRIIS checklist source")
    response=requests.get(OFFICIAL_CROSSWALK,stream=True,timeout=(30,240),
                          headers={"User-Agent":"tetsu-GBIF-GRIIS-public-COL-XR-verified-crosswalk-v6.2"})
    response.raise_for_status()
    # Use streaming gzip decompression: no publication or upload of
    # raw crosswalk file and no model-selected source names.
    with gzip.GzipFile(fileobj=response.raw) as gz:
        with io.TextIOWrapper(gz,encoding="utf-8-sig",newline="") as stream:
            # The official broad taxonomic crosswalk includes some
            # exceptionally long col:classification entries. Python's
            # default 131,072-char CSV cap is insufficient, despite the
            # schema being valid. Raise a finite defensive limit.
            csv.field_size_limit(10_000_000)
            reader=csv.DictReader(stream,delimiter="\t")
            if not reader.fieldnames or any(c not in reader.fieldnames for c in REQUIRED_COLUMNS):
                raise ValueError("Official GBIF-COL mapping TSV header schema mismatch")
            result=map_country_griis_old_keys(reader,country_old)
    if not result["strict_accepted"] and not result["synonym_accepted_additional"]:
        raise ValueError("No official crosswalk matches to old US introduced taxa")
    result["source_official_mapping_url"]=OFFICIAL_CROSSWALK
    result["source_official_README"]=OFFICIAL_README
    result["COL_source_version"]=COL_VERSION
    result["mapping_response_etag"]=response.headers.get("ETag")
    result["mapping_response_last_modified"]=response.headers.get("Last-Modified")
    result["mapping_downloaded_and_streamed_from_GBIF_official_source"]=True
    return result


def compare_original_public_event_rows(old_source,translated):
    import duckdb
    db=duckdb.connect(database=":memory:")
    db.execute("INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1';")
    db.execute("SET threads=2; SET memory_limit='3GB';")
    try:
        query,fields,frame,meta=source_query(db)
        cur=db.execute(query)
        def source():
            while True:
                chunk=cur.fetchmany(BATCH_SIZE)
                if not chunk:break
                yield from chunk
        combined=translated["strict_accepted"]|translated["synonym_accepted_additional"]
        rows,statistics=audit_rows(source(),fields,combined)
    finally:
        db.close()
    records=sum(sum(c["country_GRIIS_specieskeys_original_row_counts"].values())
                for d in rows.values() for c in d.values())
    strict_id_set=translated["strict_accepted"]
    outperiods={}
    for p,grid in rows.items():
        outperiods[p]={
            "n_source_50km_cells_with_crosswalk_COUNTRY_GRIIS_records":
                sum(c["n_original_source_rows_GRIIS_PRESENCE"]>0 for c in grid.values()),
            "n_original_event_rows_matching_STRICT_current_COL_crosswalk":
                sum(sum(v for k,v in c["country_GRIIS_specieskeys_original_row_counts"].items()
                        if k in strict_id_set) for c in grid.values()),
            "n_original_event_rows_matching_STRICT_plus_synonym_resolution":
                sum(c["n_original_source_rows_GRIIS_PRESENCE"] for c in grid.values()),
            "n_original_event_rows_STRICT_COL_ID_with_known_uncertainty_le5km":
                sum(sum(v for k,v in c["country_GRIIS_specieskeys_precise_5000m_counts"].items()
                        if k in strict_id_set) for c in grid.values()),
            "n_original_event_rows_STRICT_plus_synonym_with_uncertainty_le5km":
                sum(sum(c["country_GRIIS_specieskeys_precise_5000m_counts"].values())
                    for c in grid.values()),
            "n_source_cells_STRICT_COL_GRIIS_at_least5taxa_and20_records":
                sum(sum(v for k,v in c["country_GRIIS_specieskeys_original_row_counts"].items()
                        if k in strict_id_set)>=20 and len(
                        [k for k in c["country_GRIIS_specieskeys_original_row_counts"]
                         if k in strict_id_set])>=5 for c in grid.values()),
            "source_GRIIS_original_50km_event_species":grid
        }
    result={
        "utc":datetime.now(timezone.utc).isoformat(),
        "geography":BLOCK,
        "country":"US",
        "source_actual_GBIF_Aug2026_Parquet_object":SOURCE_FILE,
        "underlying_GBIF_original_download_DOI":UNDERLYING_SOURCE_DOI,
        "original_source_legacy_GRIIS_species_count":len(old_source),
        "actual_original_Parquet_rows_streamed":statistics.get(
            "raw_source_rows_in_original_bounding_rectangle",0),
        "actual_original_Parquet_rows_numeric_speciesKey":statistics.get(
            "source_numeric_speciesKey_records",0),
        "actual_original_Parquet_rows_non_numeric_speciesKey":statistics.get(
            "source_NONNUMERIC_speciesKey_records",0),
        "country_GRIIS_original_event_rows_after_crosswalk":records,
        "periods":outperiods,
        "original_source_ledger":statistics,
        "only_ONE_Parquet_file_from_development_snapshot":True,
        "not_GBIF_original_FULL_DOI_archive":True,
        "no_railway_seed_dispersal_or_active_freight_inference":True,
    }
    return result


def sanitized_crosswalk_metadata(data):
    """Never publish millions of old-new taxon name mappings as a giant artifact."""
    return {k:v for k,v in data.items() if k not in (
        "strict_accepted","synonym_accepted_additional",
        "source_old_to_strict_COL","source_old_to_synonym_resolution"
    )}


def run(checklists):
    legacy=valid_griis(checklists).get("US")
    if not legacy:
        raise ValueError("No real accepted US old country GRIIS species backbone keys")
    lookup=fetch_official_mapping(legacy)
    result=compare_original_public_event_rows(legacy,lookup)
    result["official_crosswalk"]=sanitized_crosswalk_metadata(lookup)
    result["taxonomic_interpretation_limits"]=[
        "GBIF old country GRIIS IDs from frozen 2023 Backbone are NOT equivalent to COL XR keys without a verified crosswalk.",
        "This is May 2026 COL26.5 XR mapping, while public original data preview is Aug2026 and may include later taxonomy; nonmatched taxa remain unknown.",
        "Only strict accepted old-GBIF-to-current-COL species mapping is primary; synonym acceptedID associations are a marked sensitivity, not proof of same species concept.",
        "Crosswalk may be many-to-one and does NOT give independent biological species records for all old GRIIS accepted names.",
        "One lexically preselected August2026 public A5 Plantae parquet object is not the complete original GBIF source for the 16 50km cells or both historical periods.",
        "Mapped plant source events alone do not establish railway dispersal or any real active historic freight train vector.",
    ]
    print(json.dumps({
        "complete":True,"old_griis_species_keys":len(legacy),
        "official_crosswalk":result["official_crosswalk"],
        "raw_source_rows":result["actual_original_Parquet_rows_streamed"],
        "raw_source_legacy_numeric_keys":result["actual_original_Parquet_rows_numeric_speciesKey"],
        "raw_source_COL_alphanumeric_keys":result["actual_original_Parquet_rows_non_numeric_speciesKey"],
        "Griis_original_PRESENCE_rows_after_crosswalk":result["country_GRIIS_original_event_rows_after_crosswalk"],
        "source_period_summaries":{p:{k:v for k,v in x.items()
                                     if k!="source_GRIIS_original_50km_event_species"}
                                  for p,x in result["periods"].items()},
        "not_full_population_or_active_train_effect":True
    },ensure_ascii=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--griis-keys",required=True)
    p.add_argument("--output",default="results/audit/GBIF_public_original_GRIIS_COL_crosswalk_US_v62.json")
    a=p.parse_args()
    checklists=json.loads(Path(a.griis_keys).read_text(encoding="utf-8"))
    result=run(checklists)
    out=Path(a.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")


if __name__=="__main__":
    main()
