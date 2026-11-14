#!/usr/bin/env python3
"""v9.6: Source-only conservative UPPER BOUND on same-observer same-site
country-GRIIS botanical RECORD support across original GBIF 48 US states.

Input only de-identified, independently verified v9.5 iNaturalist
original GBIF occurrence audit (nonrandom offsets 0 and 1000, 200
records/state/period). No original individual IDs, precise coordinates,
names or occurrence rows are retained or recreated.

For a later record to be an eligible *same-observer AND same-coarse-site
GRIIS record*, it MUST BOTH (1) share observer+site with an early-source
record, AND (2) be a GRIIS-recorded species in a late source site that
has an earlier record. The actual count is unidentified because v9.5
does not retain the joint taxon×observer×tile microdata.

Thus min(count1,count2) is a mathematically valid UPPER BOUND; if this
necessary (not sufficient) bound < K=20, a K20 repeated-observer-site
flora comparison cannot exist IN THIS bounded nonrandom source probe.

K=20 established in v8.6–v9.4. The >=40 states source macro gate is
already frozen in earlier 48-state study. No ecological rail freight
model is fitted unless source necessity is met for >=40 states.

Not actual railway seed transport, state-alien status, observer detection
equivalence, plant establishment, mowing, wind or adaptive evolution.
"""
from __future__ import annotations
import argparse,json
from datetime import datetime,timezone
from pathlib import Path
from statistics import median

K=20
MIN_MACRO_STATES=40
SOURCE_STATES=48
SOURCE_KM=(10,50)


def support_by_state(v):
    if not isinstance(v,dict) or "overlap" not in v:
        raise ValueError("V9.5 original aggregated source state not complete")
    periods=("source_pre_2018_21","source_post_2023_25")
    if any(v.get(p,{}).get("n_original_rows")!=200 for p in periods):
        raise ValueError("V9.5 source state sampled record frame changed")
    out={}
    for km in SOURCE_KM:
        tile=v["overlap"]["coarse_site_overlap"][str(km)+"km"]
        repeated=tile["n_late_records_with_same_observer_same_coarse_site_as_early"]
        source_griis=tile["n_late_GRIIS_records_in_repeated_geographical_cells"]
        site_late=tile["n_late_records_in_early_sampled_geographical_cells"]
        n_all=tile["n_late_records_with_original_valid_coordinates"]
        if any(type(n)!=int or n<0 or n>200 for n in (repeated,source_griis,site_late,n_all)):
            raise ValueError("V9.5 original source overlap count not a bounded nonnegative integer")
        if repeated>site_late or source_griis>site_late:
            raise ValueError("Observer+site or GRIIS repeated-tile records exceed original repeat site sample")
        ceiling=min(repeated,source_griis)
        out[str(km)+"km"]={
            "n_late_all_Plantae_same_observer_and_coarse_site_as_early":repeated,
            "n_late_GRIIS_country_listed_source_records_at_earlier_sampled_sites":
                source_griis,
            "n_late_all_Plantae_source_records_at_earlier_sampled_sites":site_late,
            "n_late_valid_original_coordinate_records":n_all,
            "upper_bound_GRIIS_late_records_same_observer_AND_same_site":ceiling,
            "necessary_not_sufficient_K20_GRIIS_source_record_support":ceiling>=K,
            "actual_joint_GRIIS_taxon_observer_site_count_is_NOT_identified":True,
            "n_valid_species_keys_repeat_observer_and_same_site_UNAVAILABLE":None,
        }
    return out


def run(source):
    if (source.get("complete") is not True or
        len(source.get("states",{}))!=SOURCE_STATES or
        source.get("no_raw_observer_names_ids_record_ids_or_latlon_saved") is not True or
        source.get("nonrandom_source_probe_DO_NOT_EXTRAPOLATE_TO_ALL_INAT_USERS") is not True or
        source.get("bounded_deterministic_search_offsets") not in ([0,1000],(0,1000))):
        raise ValueError("Original v9.5 observer and same-site aggregate source provenance incomplete")
    by_state={s:support_by_state(v) for s,v in source["states"].items()}
    out={
        "complete":True,
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_input":"v9.5 actual iNaturalist GBIF occurrence/search aggregate 48-state nonrandom source probe",
        "source_study":"Strict K20 same-observer and same coarse site GRIIS record availability UPPER BOUND, no inferred botanical change",
        "source_page_design":"2 deterministic offsets 0/1000 with 100 original records/offset per US state × before/after period",
        "source_frozen_US_state_count":SOURCE_STATES,
        "min_original_GRIIS_indexed_source_records_per_group":K,
        "minimum_original_qualified_macro_geographic_state_sources":MIN_MACRO_STATES,
        "by_state_source_support":by_state,
        "summary_by_coarse_site_km":{},
        "do_not_create_freight_plant_effect_when_source_gate_fails":True,
        "NO_transport_or_stress_or_evolution_inference":True,
    }
    for km in SOURCE_KM:
        label=str(km)+"km"
        entries=[values[label] for values in by_state.values()]
        observer=[x["n_late_all_Plantae_same_observer_and_coarse_site_as_early"] for x in entries]
        taxa=[x["n_late_GRIIS_country_listed_source_records_at_earlier_sampled_sites"] for x in entries]
        both=[x["upper_bound_GRIIS_late_records_same_observer_AND_same_site"] for x in entries]
        sites=[x["n_late_all_Plantae_source_records_at_earlier_sampled_sites"] for x in entries]
        necessary=sum(x>=K for x in both)
        out["summary_by_coarse_site_km"][label]={
            "n_states_original_audit_complete":len(entries),
            "sum_late_original_records_matching_early_observer_and_coarse_site":sum(observer),
            "median_state_late_original_records_matching_early_observer_and_site":median(observer),
            "median_state_late_GRIIS_source_records_at_early_sampled_site":median(taxa),
            "median_state_late_allplant_records_at_early_sampled_site":median(sites),
            "n_states_with_20_allplant_late_records_same_observer_and_site":sum(x>=K for x in observer),
            "n_states_with_20_GRIIS_late_records_at_early_sampled_site":sum(x>=K for x in taxa),
            "n_states_for_which_K20_GRIIS_same_observer_same_site_is_EVEN_POSSIBLE_UPPER_BOUND":necessary,
            "n_states_definitely_LACKING_K20_GRIIS_same_observer_same_site_support_in_this_NONRANDOM_PROBE":len(both)-necessary,
            "has_at_least_40_eligible_states":necessary>=MIN_MACRO_STATES,
            "no_known_source_same_observer_site_specific_GRIIS_species_taxon_lists":True,
        }
    if any(x["has_at_least_40_eligible_states"] for x in out["summary_by_coarse_site_km"].values()):
        raise ValueError("Source upper-bound unexpectedly reaches 40 states; requires true joint source taxon verification")
    out["macro_observer_site_matched_2022_FAF6_GRIIS_plant_effect_estimable"]=False
    out["limitations"]=[
        "This is a conservative mathematical necessary-support upper bound, not an observed count of same-observer+site GRIIS plant records or actual repeated species.",
        "Two ORIGINAL GBIF search offsets per source state×time period are NONRANDOM and sparse relative to tens of thousands of actual records; source probe insufficiency is NOT proof full iNaturalist/GBIF data cannot support controlled surveys.",
        "Same observed observerID/name may not identify same individual reliably across historical GBIF data, and coarse 10/50km tiles cannot identify actual vegetation survey sites.",
        "A source repeat observer-site record may not have any matched GRIIS taxon or repeat detection, even when the upper bound reaches 20.",
        "Source year bins are 4 versus 3 years and sample source publisher order is uncontrolled; GBIF first recorded never guarantees biological immigration or known prior absence.",
        "Country GRIIS introduction list is not state-specific alien status; 2022 FAF is estimated economic cargo OD not actual rail route or viable propagule carriage.",
        "No observed train slipstream/gust, railway mowing, seed mass or plant plasticity/genetics; no causal railway ecological interpretation is possible.",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original95",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    data=json.loads(Path(a.original95).read_text(encoding="utf-8"))
    out=run(data)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "K20_observer_and_site_necessary_upper_bound_geographic_support":{
            km:{"max_possible_states":v[
                "n_states_for_which_K20_GRIIS_same_observer_same_site_is_EVEN_POSSIBLE_UPPER_BOUND"],
                "median_late_allplant_matched_observer_site_records":v[
                "median_state_late_original_records_matching_early_observer_and_site"],
                "sum_late_source_matched_observer_site_records":v[
                "sum_late_original_records_matching_early_observer_and_coarse_site"]}
            for km,v in out["summary_by_coarse_site_km"].items()},
        "macro_K20_matched_observer_site_freight_effect_estimable":
            out["macro_observer_site_matched_2022_FAF6_GRIIS_plant_effect_estimable"],
        "no_causal_seed_or_adaptation_claim":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
