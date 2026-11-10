#!/usr/bin/env python3
"""v3.1 actual GRIIS-introduced speciesKey incidence from 50km equal-area cells.

Outcome-blind source windows are the pre-frozen macro block grids (v2.9),
not reselected after observing botanical or railway outcomes. Only pairs
of periods with >=100 indexed source plant records and >=3 provider
datasets in BOTH periods are queried for fully paginated GBIF speciesKey
facets. Regional GRIIS-listed plant keys are validated source-country
introduced taxa; all other keys have UNKNOWN status, never native.

These are GBIF indexed record aggregates, not DOI-backed original event
rows or verified population occupancy. Railway/connectivity effect is NOT
estimated by this source acquisition script.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import requests
import time

from audit_gbif_equalarea_macroblocks_v29 import (
    BLOCKS, PILOT, PERIODS, make_grid
)
from audit_griis_country_backbone_keys_v24 import DATASETS
from analyze_griis_alien_rail_source_v25 import valid_griis
from harvest_gbif_cell_species_facets_v20 import (
    GBIF, parse_page, source_json, FACET_PAGE, MAX_FACET_PAGES
)

MIN_OBSERVED_ALIEN_KEYS=5
MIN_INDEXED_ALIEN_RECORDS=20


def species_params(cell,period,offset):
    if period not in PERIODS or offset<0 or offset%FACET_PAGE:
        raise ValueError("Unfrozen macro species facet request")
    return {
        "country":cell["country"],
        "geometry":cell["gbif_geometry_wgs84"],
        "kingdomKey":6,"taxonRank":"SPECIES",
        "occurrenceStatus":"PRESENT",
        "year":f"{period[0]},{period[1]}",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "limit":0,"facet":"speciesKey",
        "facetLimit":FACET_PAGE,"facetOffset":offset,
        "speciesKey.facetLimit":FACET_PAGE,"speciesKey.facetOffset":offset,
    }


def full_macro_species(client,cell,period,alienkeys,request=source_json,sleep=time.sleep):
    if not isinstance(alienkeys,set) or not alienkeys:
        raise ValueError("Country GRIIS introduced species list is unavailable")
    records=None
    seen={}
    pages=0
    for page in range(MAX_FACET_PAGES):
        _,items=parse_page(request(client,species_params(cell,period,page*FACET_PAGE)),"speciesKey")
        payload=request  # intentionally never used as source data
        # Fetch each actual page once. The source response count is
        # returned by parse_page; do not run a second request here.
        # A source count is instead computed from first parsed page below.
        raise AssertionError("Internal caller expects full source count")
    raise AssertionError("unreachable")


def source_species_paginated(client,cell,period,alienkeys,request=source_json,sleep=time.sleep):
    total=None;seen={};page_count=0
    if not alienkeys:
        raise ValueError("Missing GRIIS species key set")
    for page in range(MAX_FACET_PAGES):
        count,items=parse_page(
            request(client,species_params(cell,period,page*FACET_PAGE)),
            "speciesKey"
        )
        if total is None:total=count
        elif count!=total:
            raise ValueError("GBIF indexed source record count changed during paging")
        for key,n in items:
            if key in seen:
                raise ValueError("Overlapping GBIF species facet pages; do not invent complete country flora")
            seen[key]=n
        page_count+=1
        if len(items)<FACET_PAGE:
            break
        sleep(.1)
    else:
        raise ValueError("Species facet exceeds capped pagination; source incomplete")
    if total is None or sum(seen.values())>total:
        raise ValueError("GBIF species facet sum exceeds indexed source")
    aliens={k:n for k,n in seen.items() if k in alienkeys}
    record_n=sum(aliens.values())
    return {
        "indexed_all_plantae_source_records":total,
        "n_all_observed_species_key_categories":len(seen),
        "n_observed_GRIIS_introduced_plant_species_keys":len(aliens),
        "n_indexed_GRIIS_introduced_plant_records":record_n,
        "observed_alien_specieskey_counts":aliens,
        "n_species_facet_source_pages":page_count,
        "full_facet_source_pagination_completed":True,
        "alien_observed_minimum_5_species_and_20_records":(
            len(aliens)>=MIN_OBSERVED_ALIEN_KEYS and record_n>=MIN_INDEXED_ALIEN_RECORDS
        ),
        "nonlisted_species_status_unknown_not_native":True,
    }


def run(blocks,source,checklists,client=None,request=source_json,sleep=time.sleep):
    if not source.get("complete") or source.get("source_errors"):
        raise ValueError("Source 50km GBIF record opportunity incomplete")
    if set(blocks)-set(source.get("blocks",{})):
        raise ValueError("New 50km blocks absent from pre-outcome source frame")
    if not blocks or any(b not in BLOCKS for b in blocks) or len(set(blocks))!=len(blocks):
        raise ValueError("Invalid frozen macro sampling blocks")
    allowed=valid_griis(checklists)
    own=client is None
    if client is None:
        client=requests.Session()
        client.headers.update({"User-Agent":"tetsu-macro-50km-GRIIS-alien-species-facet-v3.1"})
    out={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "stage":"source-identified introduced-only speciesKey GBIF faceted records in preselected equal-area 50km macro grids",
        "requested_blocks":list(blocks),
        "periods":[list(p) for p in PERIODS],
        "regions":{},"source_errors":{},
        "not_a_railway_ecological_effect":True,
    }
    for block in blocks:
        country=BLOCKS[block][0]
        keys=allowed.get(country)
        if keys is None:
            out["source_errors"][block]="Country-specific source GRIIS checklist unavailable or below quality minimum"
            continue
        old=source["blocks"][block]["cells"]
        new={}
        for cell in make_grid(block):
            cid=cell["cell_id"]
            eligibility=old.get(cid)
            if not eligibility or eligibility["country"]!=country:
                raise ValueError("Changed macro source cell identity/country")
            rows=eligibility["period_coverage"]
            eligible=all(rows.get(f"{p[0]}_{p[1]}",{}).get("source_eligible") is True for p in PERIODS)
            outcome={"cell_id":cid,"country":country,
                     "source_pre_outcome_both_period_recording_eligible":eligible,
                     "periods":{}}
            if eligible:
                for period in PERIODS:
                    code=f"{period[0]}_{period[1]}"
                    try:
                        item=source_species_paginated(client,cell,period,keys,request=request,sleep=sleep)
                        prior=rows[code]["indexed_present_plant_record_count"]
                        item["prior_2026_GBIF_indexed_source_record_count"]=prior
                        item["source_count_index_drift"]=item["indexed_all_plantae_source_records"]-prior
                        outcome["periods"][code]=item
                    except Exception as exc:
                        out["source_errors"][cid+"_"+code]=f"{type(exc).__name__}: {str(exc)[:190]}"
                        outcome["periods"][code]={"source_error":True}
                    sleep(.05)
            new[cid]=outcome
        successful=sum(
            bool(v["periods"]) and len(v["periods"])==len(PERIODS) and
            all(r.get("full_facet_source_pagination_completed") for r in v["periods"].values())
            for v in new.values()
        )
        usable_intro=sum(
            len(v["periods"])==len(PERIODS) and
            all(r.get("alien_observed_minimum_5_species_and_20_records") is True
                for r in v["periods"].values())
            for v in new.values()
        )
        out["regions"][block]={
            "country":country,
            "n_equal_area_source_cells":16,
            "n_preoutcome_both_period_source_eligible_cells":sum(x["source_pre_outcome_both_period_recording_eligible"] for x in new.values()),
            "n_complete_introduced_species_facets_two_period_cells":successful,
            "n_introduced_source_minimum_both_period_cells":usable_intro,
            "n_country_GRIIIS_backbone_keys":len(keys),
            "cells":new
        }
        print(json.dumps({"block":block,"country":country,
            "preoutcome_source_cells":out["regions"][block]["n_preoutcome_both_period_source_eligible_cells"],
            "complete_species_cells":successful,
            "introduced_eligible_both_period_cells":usable_intro,
            "source_errors":len(out["source_errors"])
        },ensure_ascii=False),flush=True)
    out["complete"]=not out["source_errors"] and len(out["regions"])==len(blocks)
    out["limits"]=[
        "Full speciesKey facets are still indexed GBIF record aggregates, not independent survey plots or exact native plant absences.",
        "Country GRIIS lists are source checklists, incomplete and different vintages; non-membership is unknown, not native.",
        "Cells crossing national boundaries/ocean are not yet filtered on terrestrial land fraction.",
        "Record opportunity >100 and >=3 dataset categories is a source minimum, not equal observer effort or deduplicated event data.",
        "No railway or road network, climate, ports or plant-trait exposure enters this source-only species acquisition.",
        "No actual train seed transport, local dispersal, true colonization or railway ecological association inferred."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gbif-coverage",required=True)
    p.add_argument("--griis-keys",required=True)
    p.add_argument("--blocks",nargs="+",choices=list(BLOCKS),
                   default=["jp_tohoku","us_midwest","fr_burgundy"])
    p.add_argument("--output",default="results/audit/griis_macro_50km_three_regions_v31.json")
    a=p.parse_args()
    old=json.loads(Path(a.gbif_coverage).read_text(encoding="utf-8"))
    griis=json.loads(Path(a.griis_keys).read_text(encoding="utf-8"))
    result=run(a.blocks,old,griis)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":result["complete"],"n_regions":len(result["regions"]),
                      "source_errors":result["source_errors"]},ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("50km GRIIS species-source acquisition incomplete")


if __name__=="__main__":
    main()
