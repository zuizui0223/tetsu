#!/usr/bin/env python3
"""v5.5 PREPARE-ONLY official GBIF DOI-bearing DWCA download predicates.

Ten original OUTCOME-FROZEN geographic blocks across eight countries,
each with EXACTLY 16 existing EPSG:6933 squares of 50×50 km (2,500 km²)
polygon filters. Use OR(within(original WGS84 WKT)) for all 16 source
cells, plus country, Plantae kingdomKey 6, species rank, 2010–2025
inclusive, coordinate=TRUE and no flagged geospatial issue.

Never use outdated GET occurrenceStatus=PRESENT / incompatible DETECTED
search filter. The full original DWCA archive must be inspected for
record-level detection/absence semantics, observer/dataset duplicates,
coordinate uncertainty and source taxonomy before biological inference.

This script NEVER posts an API request and NEVER reads, stores, creates
or prints GBIF credentials or personal data. GBIF authenticated
occurrence downloads require separate user-owned credentials. Once
actually requested and completed by the user, report the assigned
download DOI and individual dataset citations. Never claim manifests
are downloaded data or DOI publications.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path

from shapely import wkt

from audit_gbif_equalarea_macroblocks_v29 import PILOT,BLOCKS,make_grid,validate_frame

DOWNLOAD_URL="https://api.gbif.org/v1/occurrence/download/request"
HISTORICAL_START=2010
HISTORICAL_END=2025
MAX_WITHIN_SOURCE_POINTS=10000


def predicate_for(block):
    if block not in PILOT:
        raise ValueError("Source block must be one of 10 original frozen macro windows")
    country=BLOCKS[block][0]
    cells=make_grid(block)
    if len(cells)!=16 or any(c["country"]!=country or c["epsg6933_equal_area_km2"]!=2500. for c in cells):
        raise ValueError("Original 4x4 source geographic frame changed")
    within=[]
    used_ids=set()
    total_points=0
    for cell in cells:
        cid=cell["cell_id"]
        if cid in used_ids:raise ValueError("Duplicated source macro cell")
        used_ids.add(cid)
        geometry=cell["gbif_geometry_wgs84"]
        poly=wkt.loads(geometry)
        if poly.geom_type!="Polygon" or poly.is_empty or not poly.is_valid or not poly.exterior.is_ccw:
            raise ValueError("Source polygon must be valid anti-clockwise WGS84 WKT, not inverse world")
        n=len(poly.exterior.coords)
        total_points+=n
        within.append({"type":"within","geometry":geometry})
    if total_points>MAX_WITHIN_SOURCE_POINTS:
        raise ValueError("GBIF within geometry vertex budget exceeded")
    preds=[
        {"type":"equals","key":"COUNTRY","value":country},
        {"type":"equals","key":"KINGDOM_KEY","value":"6"},
        {"type":"equals","key":"TAXON_RANK","value":"SPECIES"},
        {"type":"equals","key":"HAS_COORDINATE","value":"true"},
        {"type":"equals","key":"HAS_GEOSPATIAL_ISSUE","value":"false"},
        {"type":"greaterThanOrEquals","key":"YEAR","value":str(HISTORICAL_START)},
        {"type":"lessThanOrEquals","key":"YEAR","value":str(HISTORICAL_END)},
        {"type":"or","predicates":within},
    ]
    for p in preds:
        if p.get("key")=="OCCURRENCE_STATUS":
            raise ValueError("Current source occurrenceStatus indexing API incompatibility")
    return {
        "sendNotification":False,
        "format":"DWCA",
        "description":"tetsu source-only pre-frozen 50km railway-plant macro block "+block+
            "; ORIGINAL GBIF indexed Plantae rows 2010-2025 for source QA",
        "predicate":{"type":"and","predicates":preds},
    },total_points


def prepare(outputdir):
    validate_frame()
    out=Path(outputdir);out.mkdir(parents=True,exist_ok=True)
    manifest={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "GBIF_auth_required_before_user_submits_any_download":True,
        "NO_GBIF_download_REQUEST_posted":True,
        "NO_DOI_CREATED_by_these_request_manifest_files":True,
        "GBIF_official_download_API":DOWNLOAD_URL,
        "GBIF_official_download_guide":"https://techdocs.gbif.org/en/data-use/api-downloads",
        "original_country_blocks":list(PILOT),
        "source_years_inclusive":[HISTORICAL_START,HISTORICAL_END],
        "source_output_is_PLANTAE_not_verified_introduced_species_only":True,
        "GBIF_legacy_Present_filter_deliberately_omitted":True,
        "true_country_GRIIS_annotate_ONLY_after_raw_file_is_acquired":True,
        "source_region_download_specs":{}
    }
    for block in PILOT:
        data,vertices=predicate_for(block)
        path=out/(block+"_official_DWCA_request.json")
        encoded=(json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+"\n").encode()
        path.write_bytes(encoded)
        manifest["source_region_download_specs"][block]={
            "country":BLOCKS[block][0],
            "n_exact_50km_geographic_cell_WKT_filters":16,
            "nominal_total_projected_source_area_km2":40000,
            "within_total_WKT_coordinate_vertices":vertices,
            "source_predicate_file":path.name,
            "source_request_sha256":sha256(encoded).hexdigest(),
            "download_request_submitted":False,
            "DOI_unknown_until_GBIF_completes_user_authenticated_DOWNLOAD":True
        }
        print(json.dumps({"block":block,"country":BLOCKS[block][0],
            "n_geographic_cells":16,
            "source_WKT_vertices":vertices,
            "prepare_only_no_API_post":True
        },ensure_ascii=False),flush=True)
    reg=out/"manifest_index.json"
    reg.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"complete":True,"manifest_files":len(PILOT),
        "user_authenticated_source_download_NOT_requested":True,
        "manifest_directory":str(out)},ensure_ascii=False),flush=True)
    return manifest


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--outputdir",default="results/audit/GBIF_official_DWCA_request_manifests_v55")
    a=p.parse_args()
    prepare(a.outputdir)


if __name__=="__main__":
    main()
