#!/usr/bin/env python3
"""v28.4 author ORIGINAL Pladias source semantic fingerprint audit.

Pladias downloadTraitData xlsx may regenerate volatile workbook ZIP bytes
on each request; v14.3 original dispersal/life-form XLSX SHA 135fff... and
607638... no longer match 2026-11-18 downloads, and even freshly
downloaded XLSX bytes changed on the next source fetch. NEVER silently
claim immutable XLSX source authenticity.

Canonicalize only the VALUE/Taxon mapping used by the original
strict source reader: unambiguous author Sádlo dispersion archetype and
Raunkiaer LIFE form parsed from exact original source data rows. Do not
include raw source byte-stream packaging timestamps/ZIP metadata.
SHA256 over UTF-8 canonical sorted, unambiguous strict taxon→trait map,
not an invented trait. Check original published 3147 and 4599 rows,
original source category distribution; all missing/ambiguous remain
missing. A stale or altered categorical source mapping must fail closed
once the genuine semantic digest is frozen in the original CI.

No original full Pladias taxon-trait mappings or source workbook
redistributed; report aggregate counts and content digest only.
Always cite original Pladias users and respect Pladias consortium board
permissions for substantial external dissemination.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import parse_traits


def fingerprint(dispersal:dict[str,str],life:dict[str,str])->str:
    if not dispersal or not life:
        raise ValueError("Empty original Pladias trait source values")
    payload=json.dumps({
        "source_trait_categories_version":"Pladias_Sadlo_strategy_Raunkiaer_lifeform_v1",
        "unambiguous_dispersal_taxon_to_archetype":sorted(dispersal.items()),
        "unambiguous_lifeform_taxon_to_value":sorted(life.items())
    },sort_keys=True,separators=(",",":"),ensure_ascii=False).encode("utf8")
    return hashlib.sha256(payload).hexdigest()


def run(disp:Path,life:Path,expected_content_sha:str|None=None)->dict:
    d,l,source=parse_traits(disp,life)
    if (source["source_Pladias_original_diaspore_rows"]!=3147 or
        source["source_Pladias_original_lifeform_rows"]!=4599 or
        len(d)!=3147 or len(l)!=4161):
        raise ValueError("Original Pladias trait mappings/row counts differ from independently audited v14.3 source")
    digest=fingerprint(d,l)
    if expected_content_sha is not None and digest!=expected_content_sha:
        raise ValueError("Pinned original Pladias 2026-11-18 semantic taxon→trait map SHA drift")
    return {
        "source_Pladias_updated_public_download_byte_SHA": {
            "dispersal":source["source_strategy_SHA256"],
            "lifeform":source["source_lifeform_SHA256"]
        },
        "source_Pladias_canonical_taxon_trait_semantic_SHA256":digest,
        "source_original_dispersal_rows":3147,
        "source_original_lifeform_rows":4599,
        "n_source_strict_unique_dispersal_taxon_traits":len(d),
        "n_source_strict_unique_lifeform_taxon_traits":len(l),
        "source_expected_trait_category_counts":source["source_dispersal_code_counts"],
        "canonical_payload_carries_no_workbook_ZIP_file_metadata":True,
        "compares_original_strict_taxon_trait_values_not_byte_packaging":True,
        "no_raw_source_taxon_to_trait_table_exported":True,
        "original_Pladias_rights_board_approval_may_be_needed_for_manuscript":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--dispersal",required=True,type=Path)
    p.add_argument("--lifeform",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--expect-semantic-sha256",default=None)
    args=p.parse_args()
    result=run(args.dispersal,args.lifeform,args.expect_semantic_sha256)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(result,ensure_ascii=False))


if __name__=="__main__":main()
