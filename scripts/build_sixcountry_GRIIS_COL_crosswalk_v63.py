#!/usr/bin/env python3
"""v6.3: single official GBIF Backbone -> COL XR taxon crosswalk for SIX
original macroecological country GRIIS introduced-plant species lists.

Source all country old accepted GBIF Backbone speciesKey sets from v2.4
GRIIS/GBIF checklist. One streaming pass over GBIF official May2026
COL26.5 XR taxonomic crosswalk source (download/checklistbank).
This is independent of ANY 50km rail or species-similarity outcomes.

Retain country-specific strict accepted old->COL species IDs and
separate synonym acceptedID crosswalk. Many old IDs can collapse into
same COL species, so do NOT count original keys as independent taxa.

This is a *taxonomy key bridge*, NOT source plant occurrence data or
evidence of seeds, tracks, observed alien establishments, or even
complete species matching across Aug2026 COL versions. The old 2023
Backbone and newer COL26.5 May/2026-Aug taxonomy are mismatched releases;
unmapped are UNKNOWN, NEVER 'native' or 'absent'.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import gzip,io,json
from pathlib import Path
import requests

from analyze_griis_alien_rail_source_v25 import valid_griis
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX
from verify_GRIIS_GBIF_to_COL_crosswalk_v62 import (
    map_country_griis_old_keys,sanitized_crosswalk_metadata,
    OFFICIAL_CROSSWALK,OFFICIAL_README,COL_VERSION,REQUIRED_COLUMNS
)

MAX_SOURCE_ROWS=30_000_000
MAX_SELECTED_ROWS=100_000


def source_selected_old_keys(checklists,blocks=BLOCKS_SIX):
    source=valid_griis(checklists)
    countries={}
    for block in blocks:
        if block not in BLOCKS:
            raise ValueError("Nonfrozen country source block")
        country=BLOCKS[block][0]
        if country in countries:raise ValueError("Country block duplicate not independent country")
        aliases=source.get(country)
        if not isinstance(aliases,set) or len(aliases)<100:
            raise ValueError("Country GRIIS accepted taxonomic source insufficient")
        countries[country]=aliases
    if len(countries)!=len(blocks):
        raise ValueError("Country source coverage not fixed")
    return countries


def parsed_selected_official_tsv(lines,union_old):
    """Process official TSV safely, retaining only selected old IDs.

    GBIF's classification column contains quote-like data: csv.DictReader
    can merge rows incorrectly and hit huge CSV-field limits. Official
    source is TAB separated; filter first physical field BEFORE parsing
    any of the giant classification strings.
    """
    header=next(lines).rstrip("\r\n").split("\t")
    if not header or header[0]!="gbif:ID" or any(
            name not in header for name in REQUIRED_COLUMNS):
        raise ValueError("Official source crosswalk columns unexpected")
    selected_indices={field:header.index(field) for field in REQUIRED_COLUMNS}
    matches=[];scanned=0
    for line in lines:
        scanned+=1
        if scanned>MAX_SOURCE_ROWS:
            raise ValueError("Official source too many taxonomy rows")
        key=line.split("\t",1)[0]
        if key not in union_old:continue
        pieces=line.rstrip("\r\n").split("\t")
        if len(pieces)!=len(header):
            raise ValueError("Malformed selected official crosswalk source row")
        matches.append({name:pieces[idx] for name,idx in selected_indices.items()})
        if len(matches)>MAX_SELECTED_ROWS:
            raise ValueError("Too many selected GRIIS-source crosswalk records")
    if scanned<1 or not matches:
        raise ValueError("Official taxonomic source crosswalk had no country GRIIS matches")
    return matches,scanned


def analyze(countries,selected_rows,n_scanned,provenance):
    if len(countries)!=len(BLOCKS_SIX) or not selected_rows:
        raise ValueError("Incorrect six-country original source study frame")
    output={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "source":"Official May2026 COL26.5 XR - GBIF old Backbone crosswalk",
        "official_mapping_url":OFFICIAL_CROSSWALK,
        "official_mapping_readme":OFFICIAL_README,
        "source_taxonomy_version":COL_VERSION,
        "source_full_crosswalk_rows_scanned":n_scanned,
        "n_unique_crosswalk_rows_from_union_of_country_GRIIS_old_keys":len(selected_rows),
        "source_response_provenance":provenance,
        "original_country_macro_blocks":list(BLOCKS_SIX),
        "countries":{},"summary":{},
        "not_source_original_GBIF_occurrences_or_rail_transport_evidence":True
    }
    for country,old in countries.items():
        rows=[r for r in selected_rows if r["gbif:ID"] in old]
        mapped=map_country_griis_old_keys(rows,old)
        # A post-2023 old GBIF accepted key may be absent from the
        # official 2023->2026 mapping. Count as UNMAPPED, not as a
        # biological species absence or a failed taxonomic download.
        if mapped["n_old_keys_found_in_public_mapping"]==0:
            raise ValueError(f"No old country GRIIS keys mapped at all in {country}")
        result=sanitized_crosswalk_metadata(mapped)
        result["country"]=country
        result["n_source_country_GRIIS_accepted_old_specieskeys"]=len(old)
        result["old_GRIIS_to_COL_strict_accepted_map"]=mapped["source_old_to_strict_COL"]
        result["old_GRIIS_to_COL_synonym_acceptedID_map"]=mapped["source_old_to_synonym_resolution"]
        # Number of current accepted COL IDs across accepted+synonym
        # (not total old GBIF GRIIS keys).
        result["n_distinct_source_GRIIS_current_STRICT_COL_species_IDs"]=len(mapped["strict_accepted"])
        result["n_distinct_source_GRIIS_current_STRICT_plus_synonym_COL_species_IDs"]=len(
            mapped["strict_accepted"]|mapped["synonym_accepted_additional"])
        if any(not v for v in result["old_GRIIS_to_COL_strict_accepted_map"].values()):
            raise ValueError("Empty official mapped accepted ID")
        output["countries"][country]=result
        print(json.dumps({
            "country":country,
            "old_country_GRIIS_plant_keys":len(old),
            "strict_COL_accepted_species_mapped":result["n_old_backbone_species_mapping_to_accepted_COL_species"],
            "synonym_accepted_ID_resolution":result["n_old_backbone_species_mapping_to_COL_synonyms"],
            "unmapped_or_species_concept_rejected":result["n_old_keys_unmapped_or_no_usable_COL"],
            "distinct_strict_plus_synonym_COL_IDs":result["n_distinct_source_GRIIS_current_STRICT_plus_synonym_COL_species_IDs"],
        },ensure_ascii=False),flush=True)
    output["summary"]={
        "n_independent_country_taxonomic_sources":len(output["countries"]),
        "n_old_country_GRIIS_specieskey_entries_across_country_lists":sum(len(x) for x in countries.values()),
        "n_old_country_specieskeys_direct_strict_COL_match":sum(
            r["n_old_backbone_species_mapping_to_accepted_COL_species"]
            for r in output["countries"].values()),
        "n_old_country_specieskeys_synonym_resolved":sum(
            r["n_old_backbone_species_mapping_to_COL_synonyms"]
            for r in output["countries"].values()),
        "n_old_country_specieskeys_not_mapped_or_rejected":sum(
            r["n_old_keys_unmapped_or_no_usable_COL"]
            for r in output["countries"].values()),
    }
    if sum(output["summary"][k] for k in (
        "n_old_country_specieskeys_direct_strict_COL_match",
        "n_old_country_specieskeys_synonym_resolved",
        "n_old_country_specieskeys_not_mapped_or_rejected"))!=output["summary"]["n_old_country_GRIIS_specieskey_entries_across_country_lists"]:
        raise ValueError("Country taxonomy migration source count nonconservation")
    output["limitations"]=[
        "Official mapping is May2026 COL26.5 XR; actual August2026 public GBIF DWCA preview taxonomic version may differ and requires row-level matching checks.",
        "Old country GRIIS accepted species keys may map to new accepted or synonym targets or no concept; no unlisted/nonmapped taxon is classified native.",
        "Same COL accepted ID reused from multiple old keys is ONE concept, not independent species records. Checklist completeness differs among countries.",
        "Country-matched taxonomic key arrays are a necessary source join, not original GBIF occurrence event data or verified plant establishment.",
        "No actual train/freight service or seed transfer is included, and no causal railway flora inference is possible from taxonomy mapping alone.",
    ]
    return output


def run(checklists,source=requests.get):
    countries=source_selected_old_keys(checklists)
    response=source(
        OFFICIAL_CROSSWALK,stream=True,timeout=(30,240),
        headers={"User-Agent":"tetsu-sixcountry-GRIIS-to-COL26.5XR-taxonomy-crosswalk-v6.3"}
    )
    response.raise_for_status()
    with gzip.GzipFile(fileobj=response.raw) as compressed:
        with io.TextIOWrapper(compressed,encoding="utf-8-sig",newline="") as lines:
            rows,n=parsed_selected_official_tsv(lines,set().union(*countries.values()))
    meta={
        "official_HTTP_status":response.status_code,
        "official_ETag":response.headers.get("ETag"),
        "official_LastModified":response.headers.get("Last-Modified"),
        "the_official_data_NOT_a_live_GBIF_species_match_verification":True,
    }
    return analyze(countries,rows,n,meta)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--griis-accepted-keys",required=True)
    parser.add_argument("--output",default="results/audit/sixcountry_GRIIS_GBIF_to_COL_XR_crosswalk_v63.json")
    args=parser.parse_args()
    checklists=json.loads(Path(args.griis_accepted_keys).read_text(encoding="utf-8"))
    results=run(checklists)
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,"summary":results["summary"],
        "six_country_source_mapping_saved":True,
        "no_train_seed_effect":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
