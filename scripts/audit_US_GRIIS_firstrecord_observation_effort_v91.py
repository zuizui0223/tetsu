#!/usr/bin/env python3
"""v9.1: source-record observation effort audit before interpreting 2022
FAF6 crop freight and 2018-21/2023-25 U.S.-country-GRIIS taxon
'first recorded' state speciesKey categories. Strictly no true
biological colonization, local introduction, train-carriage or stress
adaptation. This is a 48-Census-state SOURCE-DATA PROCESS audit.
"""
from __future__ import annotations
import argparse,json,math
from datetime import datetime,timezone
from pathlib import Path
from statistics import median
import numpy as np
from analyze_macro_50km_rail_road_aliens_v34 import _rank

PERIOD_YEARS={"before_2022":4,"after_2022":3}
MIN_STATES=40
QUALITY_RULE="Both 2018-21 and 2023-25 GBIF source windows independently eligible"
PREDICTORS=("post_all_Plantae_record_count",
            "post_before_all_Plantae_annual_record_ratio",
            "post_GRIIS_2026_country_speciesKey_richness",
            "post_before_GRIIS_annual_record_ratio")

def source_state_rows(temporal):
    if temporal.get("complete") is not True or temporal.get("source_periods")!={
        "before_2022":"2018,2021","after_2022":"2023,2025"
    } or temporal.get("source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows") is not True:
        raise ValueError("Actual pre/post original GBIF speciesKey source windows do not match 2022 freight audit")
    rows=[]
    missing=[]
    for code,info in sorted(temporal["states"].items()):
        if not info["two_year_source"]["two_period_quality"]:
            missing.append(code)
            continue
        early=info["source_windows"]["before_2022"]
        late=info["source_windows"]["after_2022"]
        if not (early.get("source_eligible") is True and late.get("source_eligible") is True):
            raise ValueError("Original quality provenance inconsistent for included state")
        ea=early["n_all_Plantae_indexed_records"]
        la=late["n_all_Plantae_indexed_records"]
        eg=early["n_US_country_GRIIS_indexed_records"]
        lg=late["n_US_country_GRIIS_indexed_records"]
        if any(type(n)!=int or n<=0 for n in (ea,la,eg,lg)):
            raise ValueError("Missing original 2026 GBIF source indexed record opportunity")
        early_only=info["two_year_source"]["n_early_only_source_taxon_keys"]
        late_only=info["two_year_source"]["n_late_only_first_recorded_source_taxon_keys"]
        late_species=late["n_US_country_GRIIS_observed_speciesKey_categories"]
        if (type(early_only)!=int or type(late_only)!=int or
            type(late_species)!=int or min(early_only,late_only)<0 or late_only>late_species):
            raise ValueError("GBIF source taxon speciesKey first-recorded category counts malformed")
        rows.append({
            "Census_FIPS_state_code":code,
            "n_pre_all_Plantae_indexed_records":ea,
            "n_post_all_Plantae_indexed_records":la,
            "n_pre_country_GRIIS_indexed_records":eg,
            "n_post_country_GRIIS_indexed_records":lg,
            "post_before_all_Plantae_annual_record_ratio":(la/3)/(ea/4),
            "post_before_GRIIS_annual_record_ratio":(lg/3)/(eg/4),
            "post_all_Plantae_record_count":la,
            "post_GRIIS_2026_country_speciesKey_richness":late_species,
            "post_only_first_recorded_taxon_categories_NOT_immigrations":late_only,
            "pre_only_recorded_taxon_categories_NOT_extinctions":early_only,
            "post_only_fraction_of_post_recorded_taxon_categories":late_only/late_species,
        })
    if len(rows)<MIN_STATES or len(rows)+len(missing)!=51:
        raise ValueError("Insufficient fixed Census 51-state real GBIF quality frame")
    return rows,missing


def spearman_source(x,y):
    if len(x)<MIN_STATES or len(set(x))<5 or len(set(y))<5:
        return None
    xr=_rank(np.asarray(x,float));yr=_rank(np.asarray(y,float))
    xz=xr-np.mean(xr);yz=yr-np.mean(yr)
    norm=float(np.linalg.norm(xz)*np.linalg.norm(yz))
    return round(float(xz@yz/norm),6) if norm>1e-8 else None


def descriptive(rows,missing):
    rec=lambda key:[r[key] for r in rows]
    all_ratio=rec("post_before_all_Plantae_annual_record_ratio")
    griis_ratio=rec("post_before_GRIIS_annual_record_ratio")
    first=rec("post_only_first_recorded_taxon_categories_NOT_immigrations")
    last=rec("pre_only_recorded_taxon_categories_NOT_extinctions")
    fraction=rec("post_only_fraction_of_post_recorded_taxon_categories")
    results={
        "complete":True,
        "source_analysis_year":"2026 GBIF indexed records attributed to 2018-21 / 2023-25 event years",
        "source_pre2018_21_window_years":4,
        "source_post2023_25_window_years":3,
        "source_observation_process_NOT_ecological_response":True,
        "n_qualified_Census_US_states_for_both_year_windows":len(rows),
        "n_excluded_state_codes_no_plant_absences_imputed":len(missing),
        "Census_state_codes_source_excluded":missing,
        "annual_all_Plantae_source_record_ratio_post_over_pre":{
            "min":round(min(all_ratio),6),"median":round(median(all_ratio),6),
            "max":round(max(all_ratio),6),
            "n_states_more_post_records_per_year":sum(v>1 for v in all_ratio)},
        "annual_US_country_GRIIS_source_record_ratio_post_over_pre":{
            "min":round(min(griis_ratio),6),"median":round(median(griis_ratio),6),
            "max":round(max(griis_ratio),6),
            "n_states_more_post_records_per_year":sum(v>1 for v in griis_ratio)},
        "n_states_more_post_only_recorded_taxon_keys_than_pre_only":{
            "post_only_greater":sum(a>b for a,b in zip(first,last)),
            "pre_only_greater":sum(a<b for a,b in zip(first,last)),
            "equal":sum(a==b for a,b in zip(first,last))},
        "post_only_recorded_speciesKey_categories_per_state":{
            "min":min(first),"median":median(first),"max":max(first)},
        "post_only_recorded_taxa_fraction_of_POST_recorded_taxa":{
            "min":round(min(fraction),6),"median":round(median(fraction),6),
            "max":round(max(fraction),6)},
        "source_correlations_no_causal_P_or_independent_GBIF_survey":{},
        "state_FIPS_source_audit_rows":rows,
        "NO_invasion_seed_railroad_mowing_wind_adaptation_inferred":True,
    }
    for target in ("post_only_first_recorded_taxon_categories_NOT_immigrations",
                   "post_only_fraction_of_post_recorded_taxon_categories"):
        results["source_correlations_no_causal_P_or_independent_GBIF_survey"][target]={
            name:spearman_source(rec(name),rec(target)) for name in PREDICTORS}
    results["limits"]=[
        "GBIF 2026 indexed source counts are not independently collected field observations; repeated specimen and publisher records may exist.",
        "GBIF taxon FIRST RECORDED in post-2022 event years does not indicate first biological arrival or colonization.",
        "Different 4-year early and 3-year late bins, retroactive submissions, observer and dataset recruitment may all change recording opportunities; annualized counts do not fix them.",
        "USA-country GRIIS listing is not state-specific introduced status; stateProvince is a text filter not polygon clipping.",
        "State FIPS are spatially linked and repeatedly supplied by shared GBIF publishers; correlations are descriptive only, no independent state randomization p-values or causal test.",
        "Positive correlation of late-only taxon category counts with late recording effort may arise from both true rich states and increased observation; cannot apportion causes.",
        "No actual freight car cargo propagules, transport route, established population, mowing, wind or genetic adaptation measured."
    ]
    return results


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--temporal85",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    data=json.loads(Path(a.temporal85).read_text(encoding="utf-8"))
    rows,missing=source_state_rows(data)
    results=descriptive(rows,missing)
    results["utc"]=datetime.now(timezone.utc).isoformat()
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "state_count":results["n_qualified_Census_US_states_for_both_year_windows"],
        "all_Plantae_annual_growth":results["annual_all_Plantae_source_record_ratio_post_over_pre"],
        "GRIIS_annual_growth":results["annual_US_country_GRIIS_source_record_ratio_post_over_pre"],
        "post_vs_pre_taxa":results["n_states_more_post_only_recorded_taxon_keys_than_pre_only"],
        "source_spearman":results["source_correlations_no_causal_P_or_independent_GBIF_survey"],
        "NOT_newly_immigrated_seed_species":True},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
