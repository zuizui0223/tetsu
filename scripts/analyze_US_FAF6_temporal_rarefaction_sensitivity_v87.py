#!/usr/bin/env python3
"""v8.7: 2022 FAF crop-directed RAIL vs TRUCK × 2018–21 origin
recorded USA-country GRIIS taxa × 2023–25 destination source-record
FIRST-RECORDED similarity at K=10/20/40/80/160 indexed-source
record rarefaction depths on exactly the SAME quality-state/pair cohort.

Also computes REVERSED freight OD orientation on identical botanical
source outcomes (a descriptive direction control, not a random null
or proof of train propagule carriage).

No 2022 biological arrivals observed; no state-native/alien assignment,
no track routing, mowing, gust wind or heritable adaptation.
"""
from __future__ import annotations
import argparse,json,math
from datetime import datetime,timezone
from pathlib import Path
import requests
from analyze_macro_rail_assemblages_v22 import incidence_probabilities
from analyze_FAF6_GRIIS_US_state_OD_alignment_v83 import (
    CENSUS_TIGER_URL,parse_census_centroids,original_crop_OD)
from analyze_US_FAF6_directional_source_taxa_v86 import (
    COHORTS,RESPONSES,source_states,construct_directions,cohort,partial_rank)

DEPTHS=(10,20,40,80,160)
RESPONSE="expected_20record_source_linked_first_record_fraction"
MIN_CANDIDATES=10


def state_probs_all_depths(temporal,eligible):
    out={}
    for state in eligible:
        source=temporal["states"][state]["source_windows"]
        early=source["before_2022"]["US_country_GRIIS_positive_speciesKey_record_counts"]
        late=source["after_2022"]["US_country_GRIIS_positive_speciesKey_record_counts"]
        if not early or not late:
            raise ValueError("Missing real temporal GBIF country-GRIIS taxa source")
        out[state]={}
        for depth in DEPTHS:
            if min(sum(early.values()),sum(late.values()))<depth:
                raise ValueError("Common 48-state source cohort cannot support frozen common K160")
            out[state][depth]={
              "before":incidence_probabilities(early,k=depth),
              "after":incidence_probabilities(late,k=depth),
            }
    return out


def reweight_rows(rows,eligible,prob,k,rail_OD=None):
    if k not in DEPTHS:raise ValueError("Unpredeclared source record count depth")
    rebuilt=[]
    for row in rows:
        o,d=row["o"],row["d"]
        preo=eligible[o]["early"]
        pred=eligible[d]["early"]
        cand=preo-pred
        if len(cand)<MIN_CANDIDATES:
            raise ValueError("Candidate group failed unchanged v8.6 source cohort minimum")
        pa=prob[o][k]["before"]
        pb=prob[d][k]["after"]
        den=sum(pa.get(t,0) for t in cand)
        if den<=0:raise ValueError("Undefined source-record normalized candidate first-record fraction")
        value=sum(pa.get(t,0)*pb.get(t,0) for t in cand)/den
        if not 0<=value<=1:
            raise ValueError("GBIF source record expected inclusion ratio not in [0,1]")
        updated=dict(row)
        updated["responses"]=dict(row["responses"])
        updated["responses"][RESPONSE]=value
        if rail_OD is not None:
            reverse=rail_OD.get((d,o),0.0)
            updated["x_rail"]=math.log1p(reverse)
            # Keep the same independently defined cohort membership and
            # geographic/source covariates; this is a reversed-predictor
            # diagnostic, not a reverse-date botanical colonization model.
        rebuilt.append(updated)
    return rebuilt


def source_results(rows):
    out={}
    for cohort_label in COHORTS:
        sample=cohort(rows,cohort_label)
        supported=len(sample)>=30 and len({r["o"] for r in sample})>=10 and len(
            {r["d"] for r in sample})>=10
        out[cohort_label]={
            "n_same_source_directed_pairs":len(sample),
            "n_unique_origins":len({r["o"] for r in sample}),
            "n_unique_destinations":len({r["d"] for r in sample}),
            "rail_rho_full_controls":partial_rank(sample,RESPONSE,"RAIL") if supported else None,
            "truck_rho_full_controls":partial_rank(sample,RESPONSE,"TRUCK") if supported else None,
            "supported":supported,
            "NOT_independent_state_pair_replicates":True,
        }
    return out


def run(temporal,freight,flora82,prior86,census_html):
    if (prior86.get("complete") is not True or
        prior86.get("n_both_GBIF_source_quality_US_states")!=48 or
        prior86.get("n_evaluable_directed_origin_destination_source_pairs")!=2255):
        raise ValueError("Original independently reproduced v8.6 source cohort changed")
    geo=parse_census_centroids(census_html,flora82["states"])
    rows,eligible,missing=construct_directions(temporal,freight,geo)
    if len(rows)!=2255 or len(eligible)!=48 or len(missing)!=3:
        raise ValueError("Fixed 48-state recorded-source orientation or candidate cohort changed")
    prob=state_probs_all_depths(temporal,eligible)
    results={}
    for depth in DEPTHS:
        actual=reweight_rows(rows,eligible,prob,depth)
        results[f"K{depth}"]=source_results(actual)
    for label in COHORTS:
        ref=prior86["all_results"][label]["results"][RESPONSE]
        now=results["K20"][label]
        for old,new in (("rail_full","rail_rho_full_controls"),("truck_full","truck_rho_full_controls")):
            if ref[old]!=now[new]:
                raise ValueError("Independent verified v8.6 K20 freight-source rho not exactly reproduced")
    direction_OD=original_crop_OD(freight)["RAIL"]
    reverse=reweight_rows(rows,eligible,prob,20,rail_OD=direction_OD)
    placebo=source_results(reverse)
    out={
        "complete":True,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "predeclared_record_rarefaction_depths":list(DEPTHS),
        "fixed_48_GBIF_quality_states_and_2255_source_direction_pairs":True,
        "source_origin_years":"2018–2021",
        "source_destination_first_positive_record_years":"2023–2025",
        "source_freight_estimated_base_year":2022,
        "all_depth_rail_and_truck_associations":results,
        "reversed_direction_2022_RAIL_OD_same_botanical_response":placebo,
        "K20_reproduced_verified_v86_rail_truck_coefficients":True,
        "no_detected_plant_absence_or_immigration_or_railway_adaptation":True,
        "interpretation_limits":[
           "This is indexed GBIF RECORDS, not actual immigrant plants, establishment, introductions or presence-absence-corrected species.",
           "Varying indexed record rarefaction depth does not remove observer/provider/source and geography/vegetation biases.",
           "2022 FAF cargo OD are economic estimates rather than verified railway route use or seed-carrying shipments.",
           "Reversed OD predictor is an imperfect direction comparator because geography, state freight strengths and bidirectional economic flows remain structured.",
           "All 2255 state-pair source comparisons share the same 48 state specimens and source lists; no dyadic p-value or independent ecosystem replication.",
           "GRIIS US-country introduced status may not hold for individual states; 2023-2025 first record does not establish a newly invasive plant.",
           "There is no railway mowing, train slipstream gust exposure, propagule shedding, phenotype plasticity, or genetic adaptation data.",
        ]}
    return out


def main():
    p=argparse.ArgumentParser()
    for x in ("temporal","freight","flora82","original86","output"):
        p.add_argument("--"+x,required=True)
    a=p.parse_args()
    read=lambda name:json.loads(Path(getattr(a,name)).read_text(encoding="utf-8"))
    r=requests.get(CENSUS_TIGER_URL,timeout=(20,90),
                   headers={"User-Agent":"tetsu-v87-2022-FAF6-source-record-depth-sensitivity"})
    r.raise_for_status()
    result=run(read("temporal"),read("freight"),read("flora82"),read("original86"),r.text)
    dst=Path(a.output)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "rarefaction_depths":result["predeclared_record_rarefaction_depths"],
        "full_controls":{k:v["all_temporal_source_state_directions"]["rail_rho_full_controls"]
                         for k,v in result["all_depth_rail_and_truck_associations"].items()},
        "rail_positive":{k:v["positive_RAIL_crop_OD"]["rail_rho_full_controls"]
                         for k,v in result["all_depth_rail_and_truck_associations"].items()},
        "reversed_RAIL_OD":result["reversed_direction_2022_RAIL_OD_same_botanical_response"],
        "NO_transport_or_evolutionary_evidence":True},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
