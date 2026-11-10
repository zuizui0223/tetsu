#!/usr/bin/env python3
"""v3.6: 8-country, 10-block 50km GRIIS introduced species source census.

The 10 blocks were geographically frozen in v2.9 before rail or
botanical response. Apply an independently measured Natural Earth
country polygon fraction >=80%, and prior GBIF two-period minimum
record/dataset gate BEFORE querying species facets. Requests use
country-appropriate GRIIS species-key membership; nonmembers UNKNOWN.

A country polygon surface mask is NOT a full land cover/ecoregion mask.
These are indexed GBIF positive taxon record aggregates, NOT original
occurrences, observed absences, colonization dates, train transport,
or an ecological railway effect.
"""
from __future__ import annotations
import argparse
import copy
from datetime import datetime,timezone
import json
from pathlib import Path

from audit_gbif_equalarea_macroblocks_v29 import PILOT,BLOCKS
from harvest_griis_macro_50km_facets_v31 import run as harvest_alien


def apply_country_quality(coverage,country_mask,blocks):
    if not coverage.get("complete") or coverage.get("source_errors"):
        raise ValueError("Incomplete original GBIF source")
    if not country_mask.get("complete") or country_mask.get("errors"):
        raise ValueError("Incomplete independent Natural Earth country mask")
    if not blocks or len(set(blocks))!=len(blocks) or any(x not in PILOT for x in blocks):
        raise ValueError("Non-frozen 10-region country macro study frame")
    original=coverage.get("blocks",{})
    masks=country_mask.get("regions",{})
    if not set(blocks).issubset(original) or not set(blocks).issubset(masks):
        raise ValueError("Country mask or GBIF source missing predefined region")
    out=copy.deepcopy(coverage)
    trace={}
    for block in blocks:
        country=BLOCKS[block][0]
        old=original[block]["cells"]
        src=masks[block]["cells"]
        if set(old)!=set(src) or masks[block].get("country")!=country:
            raise ValueError("Macro source geographic cell / country mismatch")
        pre=0;post=0;masked=0
        for cid,original_cell in old.items():
            if original_cell["country"]!=country:
                raise ValueError("Country source mismatch in GBIF cell")
            land=src[cid].get("country_fraction_min0p80")
            if not isinstance(land,bool):
                raise ValueError("Natural Earth quality flag missing/invalid")
            for p,v in original_cell["period_coverage"].items():
                if not isinstance(v.get("source_eligible"),bool):
                    raise ValueError("Prior source eligibility quality flag missing")
                pre+=int(v["source_eligible"])
                post+=int(v["source_eligible"] and land)
                masked+=int(v["source_eligible"] and not land)
                out["blocks"][block]["cells"][cid]["period_coverage"][p]["source_eligible"]=(
                    v["source_eligible"] and land
                )
        trace[block]={
            "country":country,
            "n_two_period_source_bins_before_country_mask":pre,
            "n_two_period_source_bins_after_country_mask":post,
            "n_source_bins_excluded_for_country_fraction":masked,
            "n_source_cells_country_80pct_and_both_period_eligible":sum(
                c["country_and_prior_GBIF_source_qualified"] for c in src.values()
            ),
        }
    return out,trace


def run(coverage,country_mask,griis,blocks):
    derived,trace=apply_country_quality(coverage,country_mask,blocks)
    output=harvest_alien(blocks,derived,griis)
    output["land_country_quality"]={
        "source":"Natural Earth 10m countries independent of rail/flower outcomes",
        "source_sha256":country_mask.get("source_sha256"),
        "minimum_country_polygon_fraction":country_mask.get("country_polygon_fraction_threshold"),
        "source_quality_trace":trace,
    }
    output["study_kind"]="8-country 10-block GRIIS-introduced source census conditional on independent land-jurisdiction and GBIF recording opportunity"
    output["source_regional_replicates_not_independent_species_records"]=True
    output["mask_applied_before_species_facets_fetched"]=True
    output["summary"]={
        "n_country_blocks":len(blocks),
        "n_predefined_50km_cells":len(blocks)*16,
        "n_country_and_two_period_source_qualified_cells":sum(
            v["n_source_cells_country_80pct_and_both_period_eligible"]
            for v in trace.values()),
        "n_two_period_complete_GRIIS_specieskey_cells":sum(
            x["n_complete_introduced_species_facets_two_period_cells"]
            for x in output["regions"].values()),
        "n_both_period_introduced_source_eligible_cells":sum(
            x["n_introduced_source_minimum_both_period_cells"]
            for x in output["regions"].values()),
    }
    return output


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--gbif-coverage",required=True)
    a.add_argument("--country-mask",required=True)
    a.add_argument("--griis-keys",required=True)
    a.add_argument("--output",default="results/audit/griis_50km_10blocks_country_mask_v36.json")
    args=a.parse_args()
    files=[json.loads(Path(f).read_text(encoding="utf-8")) for f in (
        args.gbif_coverage,args.country_mask,args.griis_keys)]
    out=run(*files,list(PILOT))
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],
                      "errors":out["source_errors"]},ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Incomplete country/GBIF/GRIIS 10-block source; do not infer railway biological effect")
if __name__=="__main__":
    main()
