#!/usr/bin/env python3
"""v6.9: verify previously unread 33 GBIF original A5 Parquet files in REAL
country-GRIIS introduced flora event counts for all original France
Burgundy 16×2500 km² equal-area 50km cells, 2018-2025.

NO outcome-driven source selection: the original 48 source object keys,
new 33 source objects, exact full-source 697577-row bbox denominator,
country and 50km polygon boundaries were independently audited in
v6.5 / v6.7 / v6.8 BEFORE GRIIS taxon source reanalysis.

Use the SAME official May2026 COL crosswalk, GRIIS 2023 old GBIF
accepted keys, deduplication, row-level PRESENT, event/date/species
and known coordinate uncertainty <=5km as the independently executed
v6.6 France 48-file source. Newly added files may only increase or
leave unchanged legitimate event/species counts, never decrease them.

The 81-file list is complete only within the v6.7 A5 polygon+centroid
candidate-ID UNION; v6.7 found a polygon-vs-centroid disagreement.
Therefore DO NOT call the result a complete regional GBIF occurrence
archive or a verified source-independent ecological railway effect.
The input GBIF A5 public Parquet is an EXPERIMENTAL development preview
derived from official August 2026 DWCA DOI 10.15468/dl.8yrbe7, not a
newly downloaded DOI. Raw source GBIF records remain on GBIF public
S3; only safe aggregated 50km source counts are stored here.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from hashlib import sha256
import json,re
from pathlib import Path

from audit_public_GBIF_original_FR_GRIIS_COL_v66 import (
    BLOCK,PREVIEW_ROOT,PUBLIC_S3_BUCKET,ORIGINAL_UNDERLYING_DOI,
    source_manifest,sql_source,audit_rows,resolve_country_species_key_set,
    MAX_BOUNDING_ROWS,BATCH_SIZE,STRICT_SOURCE_MIN_TAXA,STRICT_SOURCE_MIN_EVENTS
)

OLD_N_FILES=48
NEW_N_FILES=33
COMPLETE_CANDIDATES=81
OLD_BBOX_ROWS=671528
EXTRA_BBOX_ROWS=26049
COMBINED_BBOX_ROWS=OLD_BBOX_ROWS+EXTRA_BBOX_ROWS
HISTORICAL_PERIOD="2018_2025"


def candidate_keys(original48,all_A5meta,extra33):
    """Strict independent pre-plant taxon source selection provenance."""
    oldkeys,old=source_manifest(original48)
    if len(oldkeys)!=OLD_N_FILES or len(set(oldkeys))!=OLD_N_FILES:
        raise ValueError("Original 48 independent A5 source object identities missing")
    if old.get("n_actual_interpreted_event_source_rows_in_selected_files_and_bbox")!=OLD_BBOX_ROWS:
        raise ValueError("Old 48-source all-Plantae original bbox count drift")
    if all_A5meta.get("original_region")!=BLOCK or all_A5meta.get(
        "n_total_source_parquet_files_in_all_true_A5_polygon_intersecting_partitions")!=COMPLETE_CANDIDATES:
        raise ValueError("v6.7 country A5 polygon+centroid 81-object census missing")
    if all_A5meta.get("n_additional_source_files_not_covered_by_old_v65_query")!=NEW_N_FILES:
        raise ValueError("v6.7 original source file omission count drift")
    if extra33.get("block")!=BLOCK or extra33.get("source_country")!="FR" or extra33.get("complete") is not True:
        raise ValueError("v6.8 missing genuine France source count success")
    if extra33.get("source_query_errors") or extra33.get("n_source_files_unread_or_invalid")!=0:
        raise ValueError("v6.8 source objects must all be verified, not imputed as zero")
    if extra33.get("n_new_additional_source_Parquet_objects_successful")!=NEW_N_FILES or extra33.get(
        "2018to2025_raw_all_Plantae_original_rows_from_new_33_files")!=EXTRA_BBOX_ROWS or extra33.get(
        "2018to2025_combined_81_file_raw_all_Plantae_source_rows_in_same_bbox")!=COMBINED_BBOX_ROWS:
        raise ValueError("v6.8 independent 697577-row source denominator drift")
    extra=extra33.get("source_additional_object_count_diagnostics")
    if not isinstance(extra,list) or len(extra)!=NEW_N_FILES:
        raise ValueError("Missing complete source-verified additional Parquet object list")
    keys=[]
    for item in extra:
        key=item.get("source_Parquet_key")
        if not isinstance(key,str) or not key.startswith(PREVIEW_ROOT) or (
            not re.fullmatch(r"[A-Za-z0-9_./=-]+",key)) or not key.endswith(".parquet"):
            raise ValueError("Unapproved additional original GBIF Parquet source URI")
        keys.append(key)
        count=item.get("n_original_Plantae_rows_in_2018to2025_FR_200km_bbox")
        if isinstance(count,bool) or not isinstance(count,int) or count<0:
            raise ValueError("Missing real original GBIF source file row count")
    if len(set(keys))!=NEW_N_FILES or set(keys)&set(oldkeys):
        raise ValueError("Duplicated GBIF source file in original/added 81 partition")
    if sum(x["n_original_Plantae_rows_in_2018to2025_FR_200km_bbox"] for x in extra)!=EXTRA_BBOX_ROWS:
        raise ValueError("Additional 33 real source file row sum invalid")
    selected=sorted(oldkeys+keys)
    if len(selected)!=COMPLETE_CANDIDATES or len(set(selected))!=COMPLETE_CANDIDATES:
        raise ValueError("Real original 81 source file list is not unique")
    return selected,{
        "n_original_48_source_files":OLD_N_FILES,
        "n_all_additional_independently_scanned_source_files":NEW_N_FILES,
        "n_source_files_in_conservative_A5_geography_candidate_union":COMPLETE_CANDIDATES,
        "n_original_48_bbox_rows_independently_verified":OLD_BBOX_ROWS,
        "n_additional_33_bbox_rows_independently_verified":EXTRA_BBOX_ROWS,
        "n_all_81_bbox_rows_EXPECTED":COMBINED_BBOX_ROWS,
        "candidate_A5_polygon_vs_centroid_disagreement":
            all_A5meta.get("source_50km_polygon_plus_centroid_A5_grid_audit",{}).get(
                "polygon_cover_disagreed_with_centroid_lookup"),
        "all_81_source_keys_sha256":sha256("\n".join(selected).encode()).hexdigest(),
        "NOT_A_FULL_GBIF_RECORD_POPULATION_CENSUS":True,
    }


def compare_previous(old_event,periods):
    """Loss-of-species and precision regression gate on the IDENTICAL cells."""
    if old_event.get("geographic_block")!=BLOCK or old_event.get(
        "n_pre_outcome_verified_source_Parquet_objects")!=OLD_N_FILES or (
        old_event.get("source_only_period_analyzed")!=HISTORICAL_PERIOD):
        raise ValueError("Original 48-file COL GRIIS country-event source provenance missing")
    before=old_event.get("periods",{}).get(HISTORICAL_PERIOD)
    after=periods.get(HISTORICAL_PERIOD)
    if not before or not after or set(before["cells"])!=set(after["cells"]):
        raise ValueError("Original 16 source-cell footprint changed")
    oldcount=before["n_strict_country_GRIIS_distinct_event_species_counts"]
    newcount=after["n_strict_country_GRIIS_distinct_event_species_counts"]
    oldprec=before["n_strict_country_GRIIS_event_species_known_uncertainty_le5000m"]
    newprec=after["n_strict_country_GRIIS_event_species_known_uncertainty_le5000m"]
    if oldcount!=22435 or oldprec!=17528:
        raise ValueError("Original v6.6 published 48-file strict/precise figures changed")
    if newcount<oldcount or newprec<oldprec:
        raise ValueError("Adding real GBIF source records cannot lose country-GRIIS accepted event counts")
    n_new_taxa=0
    n_count_changed=0
    for cid,oldcell in before["cells"].items():
        newcell=after["cells"][cid]
        for key in (
            "COL_strict_country_GRIIS_speciesKey_event_counts",
            "COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts"
        ):
            earlier=oldcell[key]
            later=newcell[key]
            if any(later.get(sp,0)<n for sp,n in earlier.items()):
                raise ValueError("Adding 33 GRIIS original source files decreased a previously observed species/cell count")
            if key=="COL_strict_country_GRIIS_speciesKey_event_counts":
                n_new_taxa+=len(set(later)-set(earlier))
                n_count_changed+=sum(later.get(sp,0)!=count for sp,count in earlier.items())
    return {
        "old_source_GRIIS_strict_COL_unique_event_species":oldcount,
        "new_source_GRIIS_strict_COL_unique_event_species":newcount,
        "net_new_country_GRIIS_strict_COL_unique_event_species":newcount-oldcount,
        "old_source_known5km_uncertainty_GRIIS_events":oldprec,
        "new_source_known5km_uncertainty_GRIIS_events":newprec,
        "net_new_known5km_GRIIS_events":newprec-oldprec,
        "n_new_strict_COL_species_by_cell_occurrence_categories":n_new_taxa,
        "n_original_species_by_cell_categories_with_increased_event_count":n_count_changed,
        "not_proof_of_new_biological_plant_colonizations":True,
    }


def run(original48,fullA5,extra33,old_event,bridge,checklists):
    selected,source=candidate_keys(original48,fullA5,extra33)
    strict,taxmeta=resolve_country_species_key_set(checklists,BLOCK,bridge,False)
    extended,secondarymeta=resolve_country_species_key_set(checklists,BLOCK,bridge,True)
    if not strict or not strict.issubset(extended):
        raise ValueError("May2026 accepted COL GRIIS species namespace invalid")
    import duckdb
    db=duckdb.connect()
    db.execute("SET threads=2; SET memory_limit='3GB'; SET enable_progress_bar=false")
    db.execute("INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1'")
    try:
        sql,columns,metadata=sql_source(db,selected,expected_source_object_count=COMPLETE_CANDIDATES)
        cur=db.execute(sql)
        def stream():
            while True:
                batch=cur.fetchmany(BATCH_SIZE)
                if not batch:break
                yield from batch
        result_by_period,quality=audit_rows(
            stream(),columns,strict,extended,expected_late_bbox_count=COMBINED_BBOX_ROWS,
            max_rows=MAX_BOUNDING_ROWS
        )
    finally:
        db.close()
    # 2010-17 was neither source-file audited nor row-count certified
    # in v6.5/v6.8. Do NOT report its genuine plant records as zeros.
    result_by_period.pop("2010_2017",None)
    validation=compare_previous(old_event,result_by_period)
    summary=result_by_period[HISTORICAL_PERIOD]
    result={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_type":"REAL interpreted original GBIF August 2026 EXPERIMENTAL A5 Parquet rows, not API indexed facets",
        "source_underlying_GBIF_official_DWCA_DOI":"10.15468/dl.8yrbe7",
        "NO_new_GBIF_download_or_DOI_claim":True,
        "source_country":"FR","frozen_geographical_block":BLOCK,
        "n_original_equal_area_50km_cells":16,
        "source_period_analyzed":HISTORICAL_PERIOD,
        "source_2010to2017_period_NOT_SCANNED_or_absent":True,
        "source_pre_GRIIS_taxon_outcome_81_file_selection":source,
        "source_original_parquet_schema_and_bbox":metadata,
        "source_independently_verified_same_all_Plantae_bbox_row_count":
            quality.get("n_original_late_period_bounding_envelope_rows")==COMBINED_BBOX_ROWS,
        "source_quality":quality,
        "strict_country_GRIIS_taxon_namespace":taxmeta,
        "extended_country_GRIIS_synonym_taxon_namespace":secondarymeta,
        "change_from_source_48_to_81_files":validation,
        "n_country_GRIIS_strict_COL_event_species":summary["n_strict_country_GRIIS_distinct_event_species_counts"],
        "n_country_GRIIS_strict_COL_event_species_known_uncertainty_le5000m":
            summary["n_strict_country_GRIIS_event_species_known_uncertainty_le5000m"],
        "n_country_GRIIS_strict_event_species_key_cells_with_min5species_min20events":
            summary["n_source_cells_strict_GRIIS_at_least5_species_and20_event_species"],
        "n_country_GRIIS_strict_5km_event_species_key_cells_with_min5species_min20events":
            summary["n_source_cells_strict_GRIIS_at_least5_species_and20_5km_precise_event_species"],
        "by_period":{HISTORICAL_PERIOD:summary},
        "source_archive_NOT_complete_regional_original_events_until_A5_partition_gap_resolved":True,
        "NO_RAIL_effect_or_train_dispersal_result":True,
        "limitations":[
            "The selected 81 publicly listed source files are fully scanned ONLY within the observed A5 polygon+centroid union; v6.7 found polygon vs centroid A5 source disagreement and completeness beyond the union remains unverified.",
            "Experimental August2026 GBIF A5 Parquet source is a derivative of the official underlying DOI, not a new DOI created by tetsu.",
            "Only 2018-2025 and original 200km bbox all-Plantae source rows are independently count verified; do not infer 2010-2017 source periods are absent.",
            "Current source GBIF speciesKey COL taxonomy crosswalk may be incomplete for 409/1618 original French country-GRIIS accepted taxa; unmatched taxa are UNKNOWN not biological absence.",
            "User/observation GBIF identifiers, human recorder names, original precision coordinates and event IDs remain unexported; duplicate provider and observer biases may remain.",
            "Country GRIIS accepted plant keys and event precision provide source coverage but not a randomized botanical survey or plant seed dispersal along active railway lines.",
            "No statistically independent railway ecology effect is estimated by this source-quality stage."
        ]
    }
    print(json.dumps({
        "complete":True,"block":BLOCK,
        "source_parquet_files":len(selected),
        "independent_all_plants_bbox_source_rows":COMBINED_BBOX_ROWS,
        "actual_all_plants_bbox_source_rows":quality.get("n_original_late_period_bounding_envelope_rows"),
        "strict_country_GRIIS_event_species":result["n_country_GRIIS_strict_COL_event_species"],
        "strict_GRIIS_known5km_events":result["n_country_GRIIS_strict_COL_event_species_known_uncertainty_le5000m"],
        "net_added_strict_GRIIS_events_v66_to_v69":validation["net_new_country_GRIIS_strict_COL_unique_event_species"],
        "net_added_known5km_GRIIS_events":validation["net_new_known5km_GRIIS_events"],
        "n_50km_GRIIS_cells_min5_20":result["n_country_GRIIS_strict_event_species_key_cells_with_min5species_min20events"],
        "n_50km_precise5km_GRIIS_cells_min5_20":result["n_country_GRIIS_strict_5km_event_species_key_cells_with_min5species_min20events"],
        "A5_polygon_vs_centroid_gap_unresolved":
          source["candidate_A5_polygon_vs_centroid_disagreement"],
        "no_railway_seed_transport_inferred":True
    },ensure_ascii=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-48-source",required=True)
    p.add_argument("--metadata-81-source",required=True)
    p.add_argument("--new-33-source",required=True)
    p.add_argument("--original-48-country-GRIIS-event-source",required=True)
    p.add_argument("--sixcountry-COL-taxonomy-bridge",required=True)
    p.add_argument("--griis-original-country-keys",required=True)
    p.add_argument("--output",default="results/audit/GBIF_FR_Burgundy_original_81file_GRIIS_COL_v69.json")
    a=p.parse_args()
    load=lambda q:json.loads(Path(q).read_text(encoding="utf-8"))
    results=run(*[load(p) for p in (
        a.original_48_source,a.metadata_81_source,a.new_33_source,
        a.original_48_country_GRIIS_event_source,a.sixcountry_COL_taxonomy_bridge,
        a.griis_original_country_keys)])
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding="utf-8")


if __name__=="__main__":
    main()
