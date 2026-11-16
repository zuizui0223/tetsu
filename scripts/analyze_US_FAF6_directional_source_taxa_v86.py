#!/usr/bin/env python3
"""v8.6: 2022 FAA/BTS FAF6 state crop rail vs truck DIRECTIONAL economic
origin→destination tonnage and 2018–21 pre vs 2023–25 post GBIF
USA-country-GRIIS taxon FIRST-RECORDED source candidates.

Only source-index 'newly recorded after' conditional on origin-state
prior recorded source taxa, NEVER proven new botanical colonization,
state-specific alien species, railroad seed movement, state-specific
shipments, mowing/wind filtering or heritable adaptation.

No dyadic p-values: directed state pairs share all 48 state sources.
"""
from __future__ import annotations
import argparse,json,math
from collections import defaultdict
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import requests,numpy as np

from analyze_FAF6_GRIIS_US_state_OD_alignment_v83 import (
    CENSUS_TIGER_URL,parse_census_centroids,geodesic_km,original_crop_OD,_rank)
from analyze_macro_rail_assemblages_v22 import incidence_probabilities

K=20
MIN_BOTH_STATES=40
MIN_CANDIDATE_KEYS=10
MIN_PAIRS_ALL=100
MIN_PAIRS_RAIL=30
COHORTS=("all_temporal_source_state_directions","positive_RAIL_crop_OD",
         "positive_RAIL_and_TRUCK_crop_OD")
RESPONSES=("expected_20record_source_linked_first_record_fraction",
           "raw_source_linked_first_record_fraction")
PRIMARY_CONTROLS=(
    "log1p_truck_OD",
    "log_geographical_distance",
    "same_Census_region","same_Census_division",
    "log_before_origin_GBIF_indexed_plant_records",
    "log_before_destination_GBIF_indexed_plant_records",
    "log_after_destination_GBIF_indexed_plant_records",
    "log_origin_candidate_taxon_categories",
    "log_origin_state_land_area","log_destination_state_land_area",
    "log1p_crop_origin_rail_outgoing","log1p_crop_destination_rail_incoming",
    "log1p_crop_origin_truck_outgoing","log1p_crop_destination_truck_incoming",
)
REDUCED_CONTROLS=(
    "log1p_truck_OD","log_geographical_distance",
    "same_Census_region","same_Census_division",
    "log_before_origin_GBIF_indexed_plant_records",
    "log_before_destination_GBIF_indexed_plant_records",
    "log_after_destination_GBIF_indexed_plant_records",
    "log_origin_candidate_taxon_categories",
)
STATE_IDENTITY_IS_NOT_INDEPENDENT_PAIRWISE_REPLICATION=True

def source_states(temporal):
    if not temporal.get("complete") or temporal.get("n_qualified_source_states_gate_40") is not True:
        raise ValueError("Pre/post 2022 GBIF 51-state temporal-source quality insufficient")
    if (temporal.get("source_periods")!={"before_2022":"2018,2021","after_2022":"2023,2025"} or
        temporal.get("source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows") is not True or
        len(temporal.get("states",{}))!=51):
        raise ValueError("Temporal GBIF recorded source years or geography modified")
    states={}
    unavailable={}
    for code,row in sorted(temporal["states"].items()):
        v=row["source_windows"]
        early=v["before_2022"];late=v["after_2022"]
        if row["two_year_source"]["two_period_quality"] is not True:
            unavailable[code]={
                "early_reason":early.get("source_eligibility_failures"),
                "late_reason":late.get("source_eligibility_failures")}
            continue
        e=early["US_country_GRIIS_positive_speciesKey_record_counts"]
        l=late["US_country_GRIIS_positive_speciesKey_record_counts"]
        ne=early["n_all_Plantae_indexed_records"]
        nl=late["n_all_Plantae_indexed_records"]
        if min(ne,nl)<100 or min(sum(e.values()),sum(l.values()))<K:
            raise ValueError("Supposedly qualified source has too few records")
        states[code]={
          "early":set(e),"late":set(l),
          "early_prob20":incidence_probabilities(e,k=K),
          "late_prob20":incidence_probabilities(l,k=K),
          "effort_early":ne,"effort_late":nl}
    if len(states)<MIN_BOTH_STATES:
        raise ValueError("Insufficient 40 fixed states in both recorded-year source bins")
    return states,unavailable


def source_weights(OD):
    out=defaultdict(float);into=defaultdict(float)
    for (o,d),w in OD.items():
        if o!=d:
            out[o]+=w;into[d]+=w
    return out,into


def recorded_source_candidate(p,q):
    """At-risk source keys = origin pre recorded MINUS destination pre
    recorded, NOT true species absence at destination in nature."""
    candidates=p["early"]-q["early"]
    if len(candidates)<MIN_CANDIDATE_KEYS:
        return None
    raw=sum(t in q["late"] for t in candidates)/len(candidates)
    den=sum(p["early_prob20"].get(t,0) for t in candidates)
    if den<=0:return None
    rare=sum(p["early_prob20"].get(t,0)*q["late_prob20"].get(t,0)
             for t in candidates)/den
    if not (0<=raw<=1 and 0<=rare<=1):
        raise ValueError("Invalid GBIF taxon source inclusion probability")
    return {
      "expected_20record_source_linked_first_record_fraction":rare,
      "raw_source_linked_first_record_fraction":raw,
      "n_origin_before_recorded_taxa_not_recorded_destination_before":len(candidates),
      "n_origin_taxa_first_recorded_destination_later":sum(
          t in q["late"] for t in candidates),
      "NOT_true_source_to_dest_propagule_arrival":True}


def construct_directions(temporal,freight,geo):
    source,missing=source_states(temporal)
    od=original_crop_OD(freight)
    rail=od["RAIL"];truck=od["TRUCK"]
    ro,ri=source_weights(rail);to,ti=source_weights(truck)
    rows=[]
    for o in sorted(source):
        for d in sorted(source):
            if o==d:continue
            p,q=source[o],source[d]
            linked=recorded_source_candidate(p,q)
            if linked is None:continue
            distance=geodesic_km(geo[o],geo[d])
            if distance<=0:raise ValueError("Duplicate Census state centroids invalid")
            r=rail.get((o,d),0);t=truck.get((o,d),0)
            cov={
                "log1p_truck_OD":math.log1p(t),
                "log_geographical_distance":math.log(distance),
                "same_Census_region":int(geo[o]["region"]==geo[d]["region"]),
                "same_Census_division":int(geo[o]["division"]==geo[d]["division"]),
                "log_before_origin_GBIF_indexed_plant_records":math.log(p["effort_early"]),
                "log_before_destination_GBIF_indexed_plant_records":math.log(q["effort_early"]),
                "log_after_destination_GBIF_indexed_plant_records":math.log(q["effort_late"]),
                "log_origin_candidate_taxon_categories":
                   math.log(linked["n_origin_before_recorded_taxa_not_recorded_destination_before"]),
                "log_origin_state_land_area":math.log(geo[o]["land_sq_m"]),
                "log_destination_state_land_area":math.log(geo[d]["land_sq_m"]),
                "log1p_crop_origin_rail_outgoing":math.log1p(ro[o]),
                "log1p_crop_destination_rail_incoming":math.log1p(ri[d]),
                "log1p_crop_origin_truck_outgoing":math.log1p(to[o]),
                "log1p_crop_destination_truck_incoming":math.log1p(ti[d]),
            }
            if any(not math.isfinite(z) for z in cov.values()):
                raise ValueError("Invalid indexed GBIF source effort/Census/FAF freight covariate")
            rows.append({
                "o":o,"d":d,"x_rail":math.log1p(r),"x_truck":math.log1p(t),
                "rail_2022_estimated_thousand_short_tons":r,
                "truck_2022_estimated_thousand_short_tons":t,
                "responses":{x:linked[x] for x in RESPONSES},
                "candidate_taxon_categories":linked[
                    "n_origin_before_recorded_taxa_not_recorded_destination_before"],
                "late_only_intersect_origin_candidates":linked[
                    "n_origin_taxa_first_recorded_destination_later"],
                "covars":cov,
            })
    return rows,source,missing


def partial_rank(rows,response,mode="RAIL",control_names=PRIMARY_CONTROLS):
    if len(rows)<MIN_PAIRS_RAIL or response not in RESPONSES or mode not in ("RAIL","TRUCK"):
        return None
    x=np.asarray([r["x_rail" if mode=="RAIL" else "x_truck"] for r in rows],float)
    y=np.asarray([r["responses"][response] for r in rows],float)
    controls=list(control_names)
    if mode=="TRUCK":
        controls=[z for z in controls if z!="log1p_truck_OD"]
        controls=["x_rail"]+controls
    z=np.asarray([[r["x_rail"] if field=="x_rail" else r["covars"][field]
                    for field in controls] for r in rows],float)
    if (min(len(np.unique(x)),len(np.unique(y)))<5 or
        not all(np.all(np.isfinite(v)) for v in (x,y,z))):
        return None
    A=np.column_stack((np.ones(len(rows)),*(_rank(z[:,i]) for i in range(z.shape[1]))))
    if np.linalg.matrix_rank(A)<A.shape[1] or len(rows)-A.shape[1]<8:
        return None
    rx=_rank(x)-A@np.linalg.lstsq(A,_rank(x),rcond=None)[0]
    ry=_rank(y)-A@np.linalg.lstsq(A,_rank(y),rcond=None)[0]
    den=float(np.linalg.norm(rx)*np.linalg.norm(ry))
    return round(float(rx@ry/den),6) if den>1e-10 else None


def cohort(rows,name):
    if name==COHORTS[0]:return list(rows)
    if name==COHORTS[1]:
        return [x for x in rows if x["rail_2022_estimated_thousand_short_tons"]>0]
    if name==COHORTS[2]:
        return [x for x in rows
                if min(x["rail_2022_estimated_thousand_short_tons"],
                       x["truck_2022_estimated_thousand_short_tons"])>0]
    raise ValueError("Unpredefined directional freight source cohort")


def signed_results(rows):
    out={}
    for name in COHORTS:
        subset=cohort(rows,name)
        n_orig=len({r["o"] for r in subset})
        n_dest=len({r["d"] for r in subset})
        out[name]={
            "n_directed_state_pairs_sharing_origins_and_destinations_NOT_INDEPENDENT":len(subset),
            "n_distinct_source_origins":n_orig,
            "n_distinct_source_destinations":n_dest,
            "n_directed_source_candidate_taxon_comparisons":sum(r["candidate_taxon_categories"]
                                                               for r in subset),
            "model_supported":len(subset)>=(
                MIN_PAIRS_ALL if name==COHORTS[0] else MIN_PAIRS_RAIL)
                and min(n_orig,n_dest)>=10,
            "results":{},
        }
        for metric in RESPONSES:
            eligible=out[name]["model_supported"]
            results={
                "rail_full":partial_rank(subset,metric,"RAIL") if eligible else None,
                "truck_full":partial_rank(subset,metric,"TRUCK") if eligible else None,
                "rail_reduced":partial_rank(subset,metric,"RAIL",REDUCED_CONTROLS)
                    if eligible else None,
            }
            out[name]["results"][metric]=results
    return out


def node_deletion(rows,field,level="destination"):
    if level not in ("origin","destination"):
        raise ValueError("Only source-origin or receiving-destination state-jackknife")
    key="o" if level=="origin" else "d"
    values={}
    for code in sorted({r[key] for r in rows}):
        subset=[r for r in rows if r[key]!=code]
        values[code]=partial_rank(subset,field)
    valid=[v for v in values.values() if v is not None]
    return {
        "n_distinct_deleted_source_state_nodes":len(values),
        "n_supported_refits":len(valid),
        "rho_min":min(valid) if valid else None,
        "rho_max":max(valid) if valid else None,
        "n_positive":sum(v>0 for v in valid),
        "n_negative":sum(v<0 for v in valid),
        "NOT_independent_country_replication_or_inferential_p_value":True,
    }


def run(temporal,freight,flora82,census_html,census_bytes):
    if flora82.get("complete") is not True or set(temporal.get("states",{}))!=set(flora82.get("states",{})):
        raise ValueError("Original v8.2 all-51-state botanical source frame drift")
    geo=parse_census_centroids(census_html,flora82["states"])
    rows,source,missing=construct_directions(temporal,freight,geo)
    summaries=signed_results(rows)
    out={
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "2018_21_US_country_GRIIS_positive_records_pre_2022":True,
        "2023_25_US_country_GRIIS_positive_records_post_2022":True,
        "2022_freight_OD_ESTIMATED_not_seeds":True,
        "Census_source_sha256":sha256(census_bytes).hexdigest(),
        "n_both_GBIF_source_quality_US_states":len(source),
        "n_unsupported_state_source_ids":len(missing),
        "n_evaluable_directed_origin_destination_source_pairs":len(rows),
        "all_results":summaries,
        "jackknife_on_all_directional_pairs":{},
        "NOT_species_introduction_or_train_carried_seeds_or_mowing_wind_adaptation":True,
    }
    for node in ("origin","destination"):
        out["jackknife_on_all_directional_pairs"][node]=node_deletion(
            cohort(rows,COHORTS[0]),RESPONSES[0],node)
    out["limits"]=[
        "GBIF indexed species recorded in 2023–25 but not 2018–21 is NOT newly immigrated, newly naturalized, newly established or known absent in earlier years.",
        "Source recorded at origin before 2022 DOES NOT establish the taxon was loaded on a train or shipped within a 2022 SCTG cargo; origin/destination cargo is economic OD.",
        "2022 FAF6 is one base-year freight estimate, and 2022 plant records are omitted; post-2022 source-year bin is not a validated response to 2022 freight.",
        "State source GRIIS is an USA-wide taxon list, not confirmed locally introduced status; no origin/source-to-destination taxonomic transfer established.",
        "Early window four years versus late window three years causes uneven chance of source records. K20 rarefaction of GBIF indexed records cannot fix observer/provider/plot biases.",
        "StateProvince location is an interpreted source text field, not an independently polygon-clipped occurrence.",
        "Both directional origin/destination pairs reuse states and candidate species lists. No dyadic p-values or independent cross-national replications.",
        "Rail/road active freight pathways, seed mass, actual falling seeds, disturbance, grass cutting, wind gusts, plasticity and genetic adaptation are unmeasured."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for key in ("temporal","freight","flora82","output"):
        p.add_argument("--"+key,required=True)
    a=p.parse_args()
    load=lambda k:json.loads(Path(getattr(a,k)).read_text(encoding="utf-8"))
    req=requests.get(CENSUS_TIGER_URL,timeout=(20,90),
                     headers={"User-Agent":"tetsu-US-FAF6-temporal-state-candidates-v86"})
    req.raise_for_status()
    result=run(load("temporal"),load("freight"),load("flora82"),
               req.text,req.content)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "states":result["n_both_GBIF_source_quality_US_states"],
        "directed_pairs":result["n_evaluable_directed_origin_destination_source_pairs"],
        "model":result["all_results"],
        "leave_source":result["jackknife_on_all_directional_pairs"],
        "NOT_verified_seed_transfer_or_ecological_adaptation":True},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
