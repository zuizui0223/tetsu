#!/usr/bin/env python3
"""v9.3: all prequalified 48 US Census state source subsets, independently
verify iNaturalist research-grade GBIF dataset UUID's share of ALL
Plantae indexed records and USA-country-GRIIS-positive speciesKey
records both 2018–21 and 2023–25.

GBIF source key 50c9509d-22c7-4a22-a47d-8c48425ef4a7 is officially
iNaturalist Research-grade Observations (GBIF public technical
constants). This is a RECORDING-BIAS macro audit, not a rail freight
seed-dispersal or ecological adaptation analysis.

An all-source GRIIS indexed record minus iNaturalist source GRIIS
indexed record count is EXACT for positive record numbers when
filters and source index are stable; taxa missing from iNaturalist
alone cannot be called biologically absent or native.
"""
from __future__ import annotations
import argparse,json,math
from datetime import datetime,timezone
from pathlib import Path
from statistics import median
import requests
from harvest_US_state_GRIIS_GBIF_plant_source_v81 import (
    CENSUS_URL,parse_census_state_table)
from audit_US_GBIF_matched_provider_temporal_v92 import (
    provider_facet,provider_GRIIS_facets)
from analyze_griis_alien_rail_source_v25 import valid_griis

INAT_GBIF_DATASET_UUID="50c9509d-22c7-4a22-a47d-8c48425ef4a7"
SOURCE="iNaturalist Research-grade Observations (DOI:10.15468/ab3s5x)"
WINDOWS=("before_2022","after_2022")
MIN_STATES=40
MIN_OTHER_GRIIS_INDEXED_RECORDS_PER_WINDOW=20

def source_frame(temporal,census,accepted):
    if (temporal.get("complete") is not True or len(census)!=51 or
        temporal.get("source_periods")!={"before_2022":"2018,2021","after_2022":"2023,2025"} or
        temporal.get("source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows") is not True or
        len(temporal.get("states",{}))!=51 or len(accepted)<100):
        raise ValueError("Original pre/post GBIF source universe or US checklist changed")
    chosen={c:info for c,info in temporal["states"].items()
            if info["two_year_source"]["two_period_quality"] is True}
    if len(chosen)!=48:
        raise ValueError("Original 48 qualified US state source list changed")
    if any(census[k]["name"]!=v["state_name"] for k,v in chosen.items()):
        raise ValueError("Census FIPS state source identity drift")
    return chosen

def analyse_inat_state(baseline,provider_stages):
    dates={}
    for win in WINDOWS:
        allsource=baseline["source_windows"][win]
        stage=provider_stages[win]
        ns=stage["dataset_source"]
        subset=stage["INAT_scoped"]
        indexed=allsource["n_all_Plantae_indexed_records"]
        gi=allsource["n_US_country_GRIIS_indexed_records"]
        ni=subset["original_indexed_all_Plantae_records"]
        ngi=subset["n_US_country_GRIIS_indexed_records"]
        if not(0<ni<=indexed and 0<=ngi<=gi):
            raise ValueError("Provider source subset violates independently completed all-provider record counts")
        alltax=allsource["US_country_GRIIS_positive_speciesKey_record_counts"]
        itax=subset["US_country_GRIIS_positive_speciesKey_record_counts"]
        if not set(itax).issubset(alltax) or any(itax[k]>alltax[k] for k in itax):
            raise ValueError("Original iNaturalist speciesKey counts not subset of original 2026 GBIF indexed source")
        dates[win]={
            "iNaturalist_original_all_Plantae_record_source_count":ni,
            "all_provider_original_all_Plantae_record_source_count":indexed,
            "iNaturalist_fraction_of_all_Plantae_indexed_records":round(ni/indexed,6),
            "iNaturalist_source_country_GRIIS_record_count":ngi,
            "all_provider_source_country_GRIIS_record_count":gi,
            "iNaturalist_fraction_of_country_GRIIS_indexed_records":round(ngi/gi,6) if gi else None,
            "NON_iNaturalist_ALLPlantae_indexed_source_records":indexed-ni,
            "NON_iNaturalist_country_GRIIS_indexed_source_records":gi-ngi,
            "NON_iNaturalist_country_GRIIS_taxa_only_not_in_iNaturalist_lowerbound":
                len(set(alltax)-set(itax)),
            "iNaturalist_source_country_GRIIS_taxon_keys":itax,
            "n_country_GRIIS_species_key_categories_covered_by_INAT":len(itax),
            "provider_top100_covered_total_source_record_fraction":
                ns["top100_source_count_coverage"],
        }
    a,b=(dates[w] for w in WINDOWS)
    before=baseline["source_windows"]["before_2022"]["US_country_GRIIS_positive_speciesKey_record_counts"]
    after=baseline["source_windows"]["after_2022"]["US_country_GRIIS_positive_speciesKey_record_counts"]
    an=a["iNaturalist_source_country_GRIIS_taxon_keys"]
    bn=b["iNaturalist_source_country_GRIIS_taxon_keys"]
    original_late=set(after)-set(before)
    source_late=set(bn)-set(an)
    if not original_late:
        # No newly recorded source taxon in that state; no fraction.
        retained=None
    else:
        retained=round(len(original_late&set(bn))/len(original_late),6)
    return {
        "two_GBIF_year_source_periods":dates,
        "iNaturalist_increase_in_indexed_source_share":(
           b["iNaturalist_fraction_of_all_Plantae_indexed_records"]>
           a["iNaturalist_fraction_of_all_Plantae_indexed_records"]),
        "all_provider_late_only_country_GRIIS_recorded_taxon_categories":len(original_late),
        "iNaturalist_late_only_country_GRIIS_recorded_taxon_categories":len(source_late),
        "all_provider_late_only_categories_also_recorded_late_in_INAT":len(original_late&set(bn)),
        "fraction_all_provider_late_only_categories_also_in_iNaturalist_late_source":retained,
        "iNaturalist_only_late_categories_already_recorded_before_by_other_provider":
            len(source_late&set(before)),
        "NON_iNaturalist_20_GRIIS_source_records_available_in_BOTH_windows":
            min(a["NON_iNaturalist_country_GRIIS_indexed_source_records"],
                b["NON_iNaturalist_country_GRIIS_indexed_source_records"])>=MIN_OTHER_GRIIS_INDEXED_RECORDS_PER_WINDOW,
        "iNaturalist_dataset_key_same_in_both_periods_NOT_same_observer_or_detection":True,
        "NOT_actual_establishment_or_new_alien_species":True,
    }

def run(sess,temporal,census,accepted):
    chosen=source_frame(temporal,census,accepted)
    out={
        "complete":True,"source_utc":datetime.now(timezone.utc).isoformat(),
        "study":"Original 2026 GBIF 48-state pre/post iNaturalist research grade observation dominance and recorded new taxon source audit",
        "iNaturalist_original_GBIF_dataset_UUID":INAT_GBIF_DATASET_UUID,
        "dataset_identity":SOURCE,
        "source_original_Census_FIPS_US_state_frame":51,
        "source_original_GBIF_qualified_states":len(chosen),
        "source_record_periods":temporal["source_periods"],
        "all_original_48_state_names_fixed_before_examining_source_signal":True,
        "states":{},"errors":{},
        "NO_riding_train_seed_causality_mowing_gust_or_adaptation":True
    }
    for code,original in sorted(chosen.items()):
        item={"name":original["state_name"],"complete":False,"provider_source_windows":{}}
        try:
            for window in WINDOWS:
                n=original["source_windows"][window]["n_all_Plantae_indexed_records"]
                ds=provider_facet(sess,original["state_name"],window,n)
                k=ds["top_dataset_UUID_counts"].get(INAT_GBIF_DATASET_UUID)
                if k is None or k<50:
                    # The UUID was not in a rank-limited top100 source
                    # facet; do not assert true zero iNaturalist records.
                    raise ValueError("iNaturalist dataset UUID missing from all-Plantae top100 source at >=50 records")
                inner=provider_GRIIS_facets(sess,original["state_name"],window,
                                            INAT_GBIF_DATASET_UUID,k,accepted)
                item["provider_source_windows"][window]={
                    "dataset_source":ds,"INAT_scoped":inner}
            item["source_result"]=analyse_inat_state(original,item["provider_source_windows"])
            item["complete"]=True
        except Exception as exc:
            item["error"]=type(exc).__name__+": "+str(exc)[:180]
            out["errors"][code]=item["error"]
        out["states"][code]=item
        r=item.get("source_result",{})
        print(json.dumps({
            "FIPS":code,
            "complete":item["complete"],
            "iNat_share_early":r.get("two_GBIF_year_source_periods",{}).get("before_2022",{}).get("iNaturalist_fraction_of_all_Plantae_indexed_records"),
            "iNat_share_late":r.get("two_GBIF_year_source_periods",{}).get("after_2022",{}).get("iNaturalist_fraction_of_all_Plantae_indexed_records"),
            "allsource_late_only_GRIIS_taxa":r.get("all_provider_late_only_country_GRIIS_recorded_taxon_categories"),
            "fraction_first_recorded_taxa_also_in_iNat":r.get("fraction_all_provider_late_only_categories_also_in_iNaturalist_late_source"),
            "non_iNat_20_20_GRIIS_indexed_source_records":r.get("NON_iNaturalist_20_GRIIS_source_records_available_in_BOTH_windows"),
            "n_errors":len(out["errors"])},ensure_ascii=False),flush=True)
    out["complete"]=len(out["states"])==len(chosen) and not out["errors"]
    valid=[x["source_result"] for x in out["states"].values() if x.get("complete")]
    if len(valid)<MIN_STATES:
        raise ValueError("Fewer than 40 48-state original provider-scoped GBIF source records were verified")
    before=[v["two_GBIF_year_source_periods"]["before_2022"] for v in valid]
    after=[v["two_GBIF_year_source_periods"]["after_2022"] for v in valid]
    fractions=[v["fraction_all_provider_late_only_categories_also_in_iNaturalist_late_source"]
               for v in valid
               if v["fraction_all_provider_late_only_categories_also_in_iNaturalist_late_source"] is not None]
    out["summary"]={
        "n_original_48state_GBIF_provider_scoped_source_completed":len(valid),
        "median_iNaturalist_share_all_Plantae_before":round(median(x[
            "iNaturalist_fraction_of_all_Plantae_indexed_records"] for x in before),6),
        "median_iNaturalist_share_all_Plantae_after":round(median(x[
            "iNaturalist_fraction_of_all_Plantae_indexed_records"] for x in after),6),
        "min_max_iNaturalist_share_all_Plantae_after":[min(x[
            "iNaturalist_fraction_of_all_Plantae_indexed_records"] for x in after),
            max(x["iNaturalist_fraction_of_all_Plantae_indexed_records"] for x in after)],
        "n_states_iNaturalist_share_increased":sum(v[
            "iNaturalist_increase_in_indexed_source_share"] for v in valid),
        "median_fraction_of_allsource_late_only_GRIIS_taxa_also_in_iNaturalist_post":
            round(median(fractions),6) if fractions else None,
        "n_states_with_non_iNaturalist_at_least_20_GRIIS_records_each_year_bin":
            sum(v["NON_iNaturalist_20_GRIIS_source_records_available_in_BOTH_windows"]
                for v in valid),
        "n_states_with_provider_only_late_categories_previously_recorded_by_OTHER_data_source":
            sum(v["iNaturalist_only_late_categories_already_recorded_before_by_other_provider"]>0
                for v in valid),
        "n_GBIF_provider_filter_errors":len(out["errors"])
    }
    out["limitations"]=[
        "A SINGLE dominant iNaturalist dataset UUID is used to document SOURCE RECORDING CONCENTRATION, NOT to infer observer identity, detection bias magnitude or independent ecological replication.",
        "US-country-GRIIS alien checklist statuses are not state-specific, so indexed later-only plant taxa cannot be called newly introduced in that state.",
        "The 2026 GBIF search source is historical event-year indexed data, which can change through retrospective uploads and taxonomic backbone interpretation.",
        "Records from other data providers are counted exactly by subtraction only for indexed RECORD COUNTS, NOT for unique taxon subsets without separate full source provider facets.",
        "Sufficient non-iNaturalist record counts alone do not prove >=5 unique taxa, independent sites/observers, original verified event records or ecological colonization dates.",
        "No actual rail freight car seed content, freight route use, rail-side mowing/airflow exposure, true plant survival, phenotype plasticity or genetic adaptation measured."
    ]
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("griis","baseline85","output"):p.add_argument("--"+k,required=True)
    args=p.parse_args()
    load=lambda k:json.loads(Path(getattr(args,k)).read_text(encoding="utf-8"))
    resp=requests.get(CENSUS_URL,timeout=(20,90))
    resp.raise_for_status()
    states=parse_census_state_table(resp.text)
    accepted=valid_griis(load("griis")).get("US")
    with requests.Session() as session:
        session.headers.update({"User-Agent":"tetsu-US48-GRIIStaxa-iNaturalist-dataset-source-audit-v93"})
        output=run(session,load("baseline85"),states,accepted)
    dst=Path(args.output)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":output["complete"],"summary":output["summary"],
        "error_state_codes":list(output["errors"]),"not_seed_colonization":True},
        ensure_ascii=False),flush=True)
    if not output["complete"]:
        raise SystemExit("Not all original 48 states had reconciled same iNaturalist source data")

if __name__=="__main__":
    main()
