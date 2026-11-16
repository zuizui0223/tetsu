#!/usr/bin/env python3
"""v5.6 actual GBIF species-by-dataset source audit in original 50km cells.

For each existing prequalified six-country GRIIS-introduced-flora 50km
cell and BOTH historic periods, query official GBIF 2026 indexed
SPECIES_KEY facets twice with the IDENTICAL source geographic,
taxonomic, time and geospatial constraints:

  A: all datasetKey labels, no legacy occurrenceStatus=PRESENT GET filter
  B: the single most abundant ALL-PLANT datasetKey from the independent
     source v5.0 metadata audit, with all other constraints unchanged.

Subtract B from A for a *count-by-species complement*. This is NOT a
true introduced-only recording-effort control (the largest dataset was
selected from ALL Plantae). Not every species observed in multiple GBIF
dataset UUIDs represents independently confirmed plants. It still
permits a strict observational test of dependence on one source.

All results are indexed GBIF source-count categories rather than
DOI-backed original occurrence events. We do not claim verified
DETECTED at row level: non-detection fraction is separately bounded
using OCCURRENCE_STATUS FACET on the very same source requests. No
actual train-borne propagules, local trackside vegetation or native
species absence inferred. No occurrence IDs, human names or exact
observation locations are saved.

Six source blocks were frozen BEFORE this data-acquisition audit, not
chosen from source-effect signs. Two periods; unchanged source cells.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time

import requests
from audit_gbif_equalarea_macroblocks_v29 import make_grid, BLOCKS
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX, PERIODS
from analyze_griis_alien_rail_source_v25 import valid_griis
from harvest_griis_macro_50km_facets_v31 import species_params
from harvest_gbif_cell_species_facets_v20 import parse_page,source_json,FACET_PAGE,MAX_FACET_PAGES
from audit_gbif_detected_original_rows_v48b import parse_status_facet

MAX_SOURCE_OLD_DRIFT=.10
MAX_UNKNOWN_OR_NOT_DETECTED=.02
MIN_COUNTRY_SOURCE_ALIEN_KEYS=5
MIN_COUNTRY_SOURCE_ALIEN_RECORDS=20


def query_species(cell, period, dataset_key=None, offset=0):
    if period not in PERIODS or offset<0 or offset%FACET_PAGE:
        raise ValueError("Original source period/facet paging changed")
    p=species_params(cell,tuple(map(int,period.split("_"))),offset)
    removed=p.pop("occurrenceStatus",None)
    if removed!="PRESENT":
        raise ValueError("Legacy status contract source changed unexpectedly")
    # Current GBIF occurrence-status facet categories DETECTED /
    # NOT_DETECTED are NOT accepted as GET occurrenceStatus filters.
    p["facet"]=["speciesKey","occurrenceStatus"]
    p["occurrenceStatus.facetLimit"]=10
    if dataset_key is not None:
        if not isinstance(dataset_key,str) or not dataset_key.strip():
            raise ValueError("Missing true source datasetKey")
        p["datasetKey"]=dataset_key
    return p


def whole_source_species(client,cell,period,dataset_key=None,request=source_json,sleep=time.sleep):
    counts={}
    indexed=None
    statuses=None
    pages=0
    for i in range(MAX_FACET_PAGES):
        response=request(client,query_species(cell,period,dataset_key,i*FACET_PAGE))
        total,items=parse_page(response,"speciesKey")
        if indexed is None:
            indexed=total
            status=parse_status_facet(response)
            if status["indexed_source_count"]!=indexed:
                raise ValueError("Source species facet / DETECTED status counts disagree")
            statuses=status
        elif total!=indexed:
            raise ValueError("GBIF indexed counts drift during speciesKey paging")
        for key,n in items:
            if key in counts:
                raise ValueError("Duplicate indexed speciesKey across GBIF facet source pages")
            if not n or n<0:
                continue
            counts[key]=n
        pages+=1
        if len(items)<FACET_PAGE:
            break
        sleep(.06)
    else:
        raise ValueError("Full speciesKey pagination exceeds bounded GBIF indexed source budget")
    if indexed is None or sum(counts.values())>indexed:
        raise ValueError("GBIF speciesKey source count incomplete/invalid")
    unknown=statuses["unclassified_count"]
    non=statuses["statuses"].get("NOT_DETECTED",0)
    if unknown+non>indexed:
        raise ValueError("Invalid current GBIF indexed source detection status")
    return {
        "source_indexed_all_Plantae_records":indexed,
        "source_status_DETECTED_indexed_count":statuses["statuses"].get("DETECTED",0),
        "source_status_NOT_DETECTED_indexed_count":non,
        "source_status_unknown_indexed_count":unknown,
        "source_non_DETECTED_or_unknown_fraction":round((non+unknown)/indexed,7) if indexed>0 else None,
        "n_source_speciesKey_facet_pages":pages,
        "source_speciesKey_full_pagination":True,
        "speciesKey_count_by_all_taxa":counts,
    }


def split_top_and_other(full,top,griis_keys,old):
    total=full["source_indexed_all_Plantae_records"]
    n_top=top["source_indexed_all_Plantae_records"]
    prior=old["indexed_all_plantae_source_records"]
    if not isinstance(prior,int) or prior<=0 or total<=0 or n_top<0 or n_top>total:
        raise ValueError("GBIF source occurrence denominators impossible")
    if abs(total-prior)/prior>MAX_SOURCE_OLD_DRIFT:
        raise ValueError("Historical v3.6 to current indexed source count drift >10%")
    if full["source_non_DETECTED_or_unknown_fraction"] is None or (
        full["source_non_DETECTED_or_unknown_fraction"]>MAX_UNKNOWN_OR_NOT_DETECTED):
        raise ValueError("Indexed status outside verified DETECTED source gate")
    if n_top and (top["source_non_DETECTED_or_unknown_fraction"] is None or
                  top["source_non_DETECTED_or_unknown_fraction"]>MAX_UNKNOWN_OR_NOT_DETECTED):
        raise ValueError("Largest dataset has too many NOT_DETECTED/unknown status records")
    full_counts=full["speciesKey_count_by_all_taxa"]
    top_counts=top["speciesKey_count_by_all_taxa"]
    bad=[(k,n,full_counts.get(k,0)) for k,n in top_counts.items() if n>full_counts.get(k,0)]
    if bad:
        raise ValueError("Largest source dataset speciesKey counts exceed SAME-query total")
    country_all={k:n for k,n in full_counts.items() if k in griis_keys and n>0}
    country_top={k:n for k,n in top_counts.items() if k in griis_keys and n>0}
    country_other={k:n-country_top.get(k,0) for k,n in country_all.items()
                   if n-country_top.get(k,0)>0}
    fully_top=set(country_top)
    other=set(country_other)
    all_alien=set(country_all)
    if fully_top|other !=all_alien:
        raise ValueError("Alien GRIIS keys lost during dataset source partition")
    nsource_aliens=len(country_all)
    totals={"all":sum(country_all.values()),
            "top":sum(country_top.values()),
            "other":sum(country_other.values())}
    if totals["all"]!=totals["top"]+totals["other"]:
        raise ValueError("Source top/other taxon counts cannot partition full source")
    return {
        "full_source_current_country_GRIIS_aliens":country_all,
        "largest_dataset_country_GRIIS_aliens":country_top,
        "remaining_datasets_country_GRIIS_aliens":country_other,
        "n_alien_species_keys_all_indexed_sources":nsource_aliens,
        "n_species_keys_seen_only_in_top_dataset":len(fully_top-other),
        "n_species_keys_seen_only_in_other_datasets":len(other-fully_top),
        "n_species_keys_seen_in_top_AND_other_datasets":len(other&fully_top),
        "fraction_source_alien_species_keys_appearing_in_multiple_dataset_groups":round(
            len(other&fully_top)/nsource_aliens,7) if nsource_aliens else None,
        "n_indexed_country_GRIIS_alien_records_total":totals["all"],
        "n_indexed_country_GRIIS_alien_records_top_dataset":totals["top"],
        "n_indexed_country_GRIIS_alien_records_other_datasets":totals["other"],
        "fraction_source_GRIIS_introduced_records_from_top_dataset":round(
            totals["top"]/totals["all"],7) if totals["all"] else None,
        "source_alien_min5_species_and20records_all":len(country_all)>=5 and totals["all"]>=20,
        "source_alien_min5_species_and20records_without_top_dataset":len(country_other)>=5 and totals["other"]>=20,
        "source_alien_min5_species_and20records_top_dataset_only":len(country_top)>=5 and totals["top"]>=20,
        "dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys":True,
        "original_GRIIS_source_taxon_count":old.get("n_observed_GRIIS_introduced_plant_species_keys"),
    }


def run_one(block,flora,source_datasets,checklists,request=source_json,sleep=time.sleep):
    if block not in BLOCKS_SIX or flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("Original geographical GRIIS source not complete")
    if source_datasets.get("complete") is not True or source_datasets.get("errors"):
        raise ValueError("Original six-country GBIF source composition not complete")
    if source_datasets.get("block")!=block or source_datasets.get("country")!=BLOCKS[block][0]:
        raise ValueError("Source datasets country/block mismatch")
    if flora["regions"][block]["country"]!=BLOCKS[block][0]:
        raise ValueError("Original country source mismatched")
    gr=valid_griis(checklists)
    country_keys=gr.get(BLOCKS[block][0])
    if not country_keys:
        raise ValueError("Country GRIIS accepted plant taxon checklist unusable")
    geom={r["cell_id"]:r for r in make_grid(block)}
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-GBIF-GRIIS-introduced-dominant-dataset-leaveout-v5.6"})
    output={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_country":BLOCKS[block][0],
        "source_geo_block":block,
        "scope":"SpeciesKey source count split top ONE all-plant GBIF dataset vs remaining source datasets, in country GRIIS introduced flora",
        "source_GRIIS_taxa_from_independent_country_accepted_checklist":True,
        "full_source_status_NOT_DETECTED_marginal_checked":True,
        "source_cells":{},"source_errors":{},
    }
    for cid in sorted(geom):
        existing=flora["regions"][block]["cells"][cid]
        source=source_datasets["cells"][cid]
        eligible=(existing["source_pre_outcome_both_period_recording_eligible"] is True
                  and all(existing["periods"].get(p,{}).get(
                      "alien_observed_minimum_5_species_and_20_records") is True for p in PERIODS))
        if source["original_GRIIS_two_period_source_eligible"]!=eligible:
            raise ValueError("GBIF dataset source and original country GRIIS source eligibility disagree")
        item={"country":BLOCKS[block][0],"original_two_period_country_GRIIS_eligible":eligible,
              "periods":{}}
        if eligible:
            for period in PERIODS:
                try:
                    prior=existing["periods"][period]
                    ds=source["periods"][period]
                    if not ds.get("source_profile_valid_for_recording_bias_covariates"):
                        raise ValueError("Original datasetKey source concentration QA failed")
                    top_id=ds["top_datasetKey_source_id"]
                    all_src=whole_source_species(client,geom[cid],period,request=request,sleep=sleep)
                    top_src=whole_source_species(client,geom[cid],period,dataset_key=top_id,request=request,sleep=sleep)
                    n_ref=ds["source_indexed_all_Plantae_records"]
                    n_top_ref=ds["dataset_key_record_counts"][top_id]
                    if abs(all_src["source_indexed_all_Plantae_records"]-n_ref)>max(10,.02*n_ref):
                        raise ValueError("Total indexed count mismatches independent earlier dataset facet source >2%")
                    if abs(top_src["source_indexed_all_Plantae_records"]-n_top_ref)>max(10,.02*n_top_ref):
                        raise ValueError("Top dataset indexed count mismatch >2%")
                    split=split_top_and_other(all_src,top_src,country_keys,prior)
                    item["periods"][period]={
                        "top_dataset_selected_from_all_Plantae_source_UUID":top_id,
                        "all_source_indexed_Plantae_records":all_src["source_indexed_all_Plantae_records"],
                        "top_source_indexed_Plantae_records":top_src["source_indexed_all_Plantae_records"],
                        "all_source_DETECTED_marginal_count":all_src["source_status_DETECTED_indexed_count"],
                        "all_source_NOT_DETECTED_marginal_count":all_src["source_status_NOT_DETECTED_indexed_count"],
                        "n_speciesKey_source_facet_pages_all":all_src["n_source_speciesKey_facet_pages"],
                        "n_speciesKey_source_facet_pages_top":top_src["n_source_speciesKey_facet_pages"],
                        **split,
                    }
                except Exception as exc:
                    output["source_errors"][cid+"_"+period]=f"{type(exc).__name__}: {str(exc)[:250]}"
                    item["periods"][period]={"source_error":True}
                sleep(.08)
        output["source_cells"][cid]=item
    selected=[x for x in output["source_cells"].values() if x["original_two_period_country_GRIIS_eligible"]]
    both_complete=sum(
        all(x["periods"].get(p,{}).get("dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys")
            for p in PERIODS) for x in selected
    )
    output["source_original_GRIIS_eligible_cells"]=len(selected)
    output["n_both_period_full_total_top_and_other_source_cells"]=both_complete
    output["n_both_period_without_top_dataset_min5species_20records"]=sum(
        all(x["periods"].get(p,{}).get("source_alien_min5_species_and20records_without_top_dataset")
            for p in PERIODS) for x in selected
    )
    output["complete"]=not output["source_errors"] and both_complete==len(selected)
    from statistics import median
    vals=[v["periods"]["2018_2025"]["fraction_source_GRIIS_introduced_records_from_top_dataset"]
          for v in selected if "fraction_source_GRIIS_introduced_records_from_top_dataset" in v["periods"]["2018_2025"] and
          v["periods"]["2018_2025"]["fraction_source_GRIIS_introduced_records_from_top_dataset"] is not None]
    corroboration=[v["periods"]["2018_2025"]["fraction_source_alien_species_keys_appearing_in_multiple_dataset_groups"]
          for v in selected if "fraction_source_alien_species_keys_appearing_in_multiple_dataset_groups" in v["periods"]["2018_2025"] and
          v["periods"]["2018_2025"]["fraction_source_alien_species_keys_appearing_in_multiple_dataset_groups"] is not None]
    output["late_period_median_fraction_GRIIS_recorded_introduced_records_from_TOP_dataset"]=round(median(vals),6) if vals else None
    output["late_period_median_fraction_source_GRIIS_species_keys_in_top_AND_other_datasets"]=round(median(corroboration),6) if corroboration else None
    print(json.dumps({
        "block":block,"country":BLOCKS[block][0],
        "original_GRIIS_two_period_cells":len(selected),
        "n_complete_top_other_source_partition_cells":both_complete,
        "n_both_period_other_datasets_min5_20":output["n_both_period_without_top_dataset_min5species_20records"],
        "median_late_GRIIS_record_share_from_top_dataset":output["late_period_median_fraction_GRIIS_recorded_introduced_records_from_TOP_dataset"],
        "median_late_GRIIS_species_keys_corrob_by_two_dataset_groups":output["late_period_median_fraction_source_GRIIS_species_keys_in_top_AND_other_datasets"],
        "source_errors":len(output["source_errors"])
    },ensure_ascii=False),flush=True)
    output["limits"]=[
        "A dataset UUID is a recording-provider identifier, not independent survey sites or observers; top-other source groups do not imply independently observed biological colonization.",
        "Dominant dataset is chosen from all-Plantae source records, not from outcome-introduced taxa; no outcome-driven exclusion selection.",
        "Current no-status GBIF speciesKey facets may have marginal NOT_DETECTED; the source fraction is audited, not accurately apportioned to each alien species.",
        "Original historical GRIIS accepted speciesKey registry, 2026 spatial network and indexed GBIF source taxonomies vary over time.",
        "No DOI-backed original occurrence event/download, coordinate uncertainty, observer effort, active train/freight lines, seed traps or actual ecological absences observed.",
    ]
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--block",choices=BLOCKS_SIX,required=True)
    p.add_argument("--flora",required=True)
    p.add_argument("--source-datasets",required=True)
    p.add_argument("--griis-accepted-keys",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    src=[json.loads(Path(x).read_text(encoding="utf-8")) for x in
         (a.flora,a.source_datasets,a.griis_accepted_keys)]
    result=run_one(a.block,*src)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":result["complete"],
                      "source_block":a.block,
                      "source_errors":result["source_errors"]},ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("Dominant GBIF dataset source source-taxon audit failed; no railway/alien inference")


if __name__=="__main__":
    main()
