#!/usr/bin/env python3
"""v8.3: USA Census 2025 state centroids × FAF6 2022 domestic agricultural
RAIL/TRUCK state-OD flows × 2018–25 GBIF country-GRIIS plant recorded taxa.

Source-only *cross-sectional* associational MACRO diagnostic, NOT
observed train-borne weed propagules, local state-level introduction,
rail-specific colonization, railway management or adaptive evolution.

51 original FAF economic source state codes; 49 original GBIF quality
eligible states only. Within US geographic state pairs SHARE states,
no pairwise p-values, no independent country replications. Predeclare
K=20 source-record rarefied similarity as primary and raw Jaccard as
diagnostic. Hold physical/geography/record opportunity controls fixed.
"""
from __future__ import annotations
import argparse,io,json,math
from collections import defaultdict
from datetime import datetime,timezone
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from statistics import median
import numpy as np
import pandas as pd
import requests

from analyze_macro_50km_rail_road_aliens_v34 import _rank
from analyze_macro_rail_assemblages_v22 import (
    incidence_probabilities,expected_rarefied_overlap)

CENSUS_TIGER_URL="https://tigerweb.geo.census.gov/tigerwebmain/Files/acs26/tigerweb_acs26_state_2025_acs25_us.html"
REQUIRED_CENSUS_COLUMNS={"STATE","NAME","STUSAB","CENTLAT","CENTLON","AREALAND","REGION","DIVISION"}
METRICS=("rarefied20_expected_inclusion_overlap","raw_recorded_GRIIS_Jaccard")
COHORTS=("all_eligible_state_pairs","RAIL_OD_positive_only","RAIL_and_TRUCK_OD_both_positive")
PRIMARY_CONTROLS=(
    "log1p_truck_OD",
    "log_state_centroid_geo_distance",
    "log_allplant_record_geomean",
    "log_allplant_record_ratio",
    "log_state_land_area_geomean",
    "log_state_land_area_ratio",
    "same_census_region",
    "same_census_division",
    "log1p_rail_node_strength_geomean",
    "log1p_truck_node_strength_geomean",
)
MIN_STATES=40
MIN_STATE_PAIRS=80
MIN_NONZERO_STATE_PAIRS=30
K=20


def parse_census_centroids(html,source_states):
    tables=pd.read_html(io.StringIO(html))
    candidates=[table for table in tables if REQUIRED_CENSUS_COLUMNS.issubset(
        {str(x).upper() for x in table.columns})]
    if len(candidates)!=1:
        raise ValueError("Exact one official 2025 Census TIGER US states centroid table not found")
    df=candidates[0]
    result={}
    for row in df.to_dict(orient="records"):
        code=str(row["STATE"]).strip().zfill(2)
        if code not in source_states:
            continue
        name=str(row["NAME"]).strip()
        if name!=source_states[code]["official_Census_state_name"]:
            raise ValueError("Official FIPS Census state name changed between 2025 TIGER and v8.2")
        if code in result:
            raise ValueError("Duplicate Census FIPS code in TIGER state geometry table")
        lat,lon=float(row["CENTLAT"]),float(row["CENTLON"])
        area=float(row["AREALAND"])
        region=str(row["REGION"]).strip()
        division=str(row["DIVISION"]).strip()
        if not (-90<=lat<=90 and -180<=lon<=180 and area>0):
            raise ValueError("Census geographic source lat/lon/land area invalid")
        result[code]={"lat":lat,"lon":lon,"land_sq_m":area,
                      "region":region,"division":division}
    if set(result)!=set(source_states) or len(result)!=51:
        raise ValueError("Census geometric state source does not match FAF6 51 states")
    return result


def geodesic_km(a,b):
    lat1,lat2=math.radians(a["lat"]),math.radians(b["lat"])
    dlat=lat2-lat1
    dlon=math.radians(b["lon"]-a["lon"])
    h=math.sin(dlat/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return 6371.0088*2*math.asin(math.sqrt(min(1.,max(0.,h))))


def validated_botanical_states(source):
    if source.get("complete") is not True or source.get("date_period")!="2018_2025":
        raise ValueError("Original 51-state GBIF GRIIS recent source is not complete")
    if source.get("n_requested_state_fips_codes")!=51 or len(source.get("states",{}))!=51:
        raise ValueError("Original Census FIPS state source geographic scope mismatch")
    result={}
    discarded={}
    for state,row in source["states"].items():
        v=row["GBIF_2018_2025_GRIIS_country_taxa_source"]
        if v.get("source_eligible") is not True:
            discarded[state]=v.get("reason_not_qualified")
            continue
        counts=v.get("source_GRIIS_speciesKey_indexed_counts")
        effort=v.get("source_all_Plantae_indexed_records")
        if not isinstance(counts,dict) or len(counts)<5 or any(
            not isinstance(k,str) or type(n)!=int or n<=0 for k,n in counts.items()):
            raise ValueError("Complete GBIF GRIIS country source speciesKey counts invalid")
        if type(effort)!=int or effort<100:
            raise ValueError("Indexed GBIF ALL-Plantae source opportunity invalid")
        if sum(counts.values())<K:
            raise ValueError("Cannot rarefy source group with fewer than 20 records")
        result[state]={"counts":counts,"effort":effort,
                       "probs":incidence_probabilities(counts,k=K),
                       "country_GRIIS_not_state_introduction":True}
    if len(result)<MIN_STATES:
        raise ValueError("Fewer than 40 quality-eligible independently defined US state sources")
    return result,discarded


def original_crop_OD(source):
    if source.get("complete") is not True or source.get("source_hash_verified_vs_prior_v79") is not True or (
        source.get("source_trace",{}).get("previous_verified_v79_total_and_modes_reconciled") is not True):
        raise ValueError("Original authenticated 2022 FAF6 source OD network not verified")
    if source.get("source_year")!=2022:
        raise ValueError("Unfrozen freight source base year")
    by_mode={}
    for mode in ("RAIL","TRUCK"):
        rows=source["raw_official_state_OD_2022_crop_estimated_tonnage"][mode]
        d={}
        for row in rows:
            pair=(row["o"],row["d"])
            v=row["estimated_thousand_short_tons"]
            if pair in d or type(v) not in (int,float) or not math.isfinite(v) or v<=0:
                raise ValueError("Malformed or duplicate actual domestic crop FAF6 economic OD")
            d[pair]=float(v)
        if not d:raise ValueError("Missing official FAF6 agricultural source OD")
        by_mode[mode]=d
    return by_mode


def node_strength_nonself(od):
    strength=defaultdict(float)
    for (o,d),w in od.items():
        if o!=d:
            strength[o]+=w
            strength[d]+=w
    return strength


def state_pair_rows(botanical,geo,fa):
    by_mode=original_crop_OD(fa)
    r=by_mode["RAIL"];t=by_mode["TRUCK"]
    rs=node_strength_nonself(r)
    ts=node_strength_nonself(t)
    rows=[]
    for o,d in combinations(sorted(botanical),2):
        a,b=botanical[o],botanical[d]
        g=geodesic_km(geo[o],geo[d])
        if g<=0:raise ValueError("Two distinct US state centroids collapsed to same point")
        rail=r.get((o,d),0.)+r.get((d,o),0.)
        truck=t.get((o,d),0.)+t.get((d,o),0.)
        akey=set(a["counts"]);bkey=set(b["counts"])
        raw=len(akey&bkey)/len(akey|bkey)
        rare=expected_rarefied_overlap(a["probs"],b["probs"])
        if rare is None:raise ValueError("Expected US state introduced taxon rarity source missing")
        ea,eb=a["effort"],b["effort"]
        aa,ab=geo[o]["land_sq_m"],geo[d]["land_sq_m"]
        if min(ea,eb,aa,ab)<=0:raise ValueError("Missing source area or recording effort")
        responses={"rarefied20_expected_inclusion_overlap":rare,
                   "raw_recorded_GRIIS_Jaccard":raw}
        covars={
            "log1p_truck_OD":math.log1p(truck),
            "log_state_centroid_geo_distance":math.log(g),
            "log_allplant_record_geomean":(math.log(ea)+math.log(eb))/2,
            "log_allplant_record_ratio":abs(math.log(ea/eb)),
            "log_state_land_area_geomean":(math.log(aa)+math.log(ab))/2,
            "log_state_land_area_ratio":abs(math.log(aa/ab)),
            "same_census_region":int(geo[o]["region"]==geo[d]["region"]),
            "same_census_division":int(geo[o]["division"]==geo[d]["division"]),
            "log1p_rail_node_strength_geomean":(math.log1p(rs[o])+math.log1p(rs[d]))/2,
            "log1p_truck_node_strength_geomean":(math.log1p(ts[o])+math.log1p(ts[d]))/2,
        }
        if any(not math.isfinite(v) for v in [*covars.values(),rare,raw]):
            raise ValueError("Nonfinite source plant/network/geography comparison")
        rows.append({"a":o,"b":d,"x_rail":math.log1p(rail),"x_truck":math.log1p(truck),
                     "rail_estimated_thousand_short_tons":round(rail,6),
                     "truck_estimated_thousand_short_tons":round(truck,6),
                     "responses":responses,"covars":covars})
    return rows


def partial_rank_rail(rows,response,controls=PRIMARY_CONTROLS,mode="RAIL"):
    if len(rows)<MIN_STATE_PAIRS or response not in METRICS:
        return None
    y=np.asarray([v["responses"][response] for v in rows],float)
    x=np.asarray([v["x_rail" if mode=="RAIL" else "x_truck"] for v in rows],float)
    use_controls=[c for c in controls if not (mode=="TRUCK" and c=="log1p_truck_OD")]
    if mode=="TRUCK":use_controls=["x_rail"]+use_controls
    z=np.asarray([[v["x_rail"] if c=="x_rail" else v["covars"][c]
                   for c in use_controls] for v in rows],float)
    if len(np.unique(x))<5 or len(np.unique(y))<5 or not (
        np.all(np.isfinite(x)) and np.all(np.isfinite(y)) and np.all(np.isfinite(z))):
        return None
    A=np.column_stack([np.ones(len(rows)),*[_rank(z[:,i]) for i in range(z.shape[1])]])
    if np.linalg.matrix_rank(A)<A.shape[1] or len(rows)-A.shape[1]<10:
        return None
    rx=_rank(x)-A@np.linalg.lstsq(A,_rank(x),rcond=None)[0]
    ry=_rank(y)-A@np.linalg.lstsq(A,_rank(y),rcond=None)[0]
    norm=float(np.linalg.norm(rx)*np.linalg.norm(ry))
    if norm<1e-8:return None
    value=float(rx@ry/norm)
    if not -1.0001<=value<=1.0001:
        raise ValueError("Out-of-range residual rank coefficient")
    return round(value,6)


def analysis(rows):
    out={}
    for cohort in COHORTS:
        source=([r for r in rows] if cohort=="all_eligible_state_pairs"
                else [r for r in rows if r["rail_estimated_thousand_short_tons"]>0
                      and (cohort!="RAIL_and_TRUCK_OD_both_positive"
                           or r["truck_estimated_thousand_short_tons"]>0)])
        uniq=len({s for x in source for s in (x["a"],x["b"])})
        out[cohort]={
            "n_pairwise_source_state_comparisons_NOT_INDEPENDENT":len(source),
            "n_unique_states_in_pair_subset":uniq,
            "mode_qualifying_cohort_is_independent_of_plant_similarity":True,
            "responses":{}}
        for metric in METRICS:
            rho=partial_rank_rail(source,metric) if len(source)>=(
                MIN_STATE_PAIRS if cohort=="all_eligible_state_pairs" else MIN_NONZERO_STATE_PAIRS) else None
            if rho is None:
                out[cohort]["responses"][metric]={
                    "rail_partial_rho":None,"n_pairs":len(source),"source_support":False}
                continue
            out[cohort]["responses"][metric]={
                "rail_partial_rho":rho,"n_pairs":len(source),"source_support":True,
                "interpretation":"positive means higher modeled crop rail OD coincides with more similar recorded GRIIS-listed state taxon categories, conditional on truck/geometry/effort; NOT seed transport"}
    # Leave-one-state-out sensitivity for main state-pair cohort is a
    # perturbation diagnostic, never independent geographic replication.
    main=out["all_eligible_state_pairs"]["responses"]["rarefied20_expected_inclusion_overlap"]["rail_partial_rho"]
    leave={}
    for state in sorted({s for x in rows for s in (x["a"],x["b"])}):
        remaining=[r for r in rows if state not in (r["a"],r["b"])]
        leave[state]=partial_rank_rail(remaining,"rarefied20_expected_inclusion_overlap")
    valid=[v for v in leave.values() if v is not None]
    out["leave_one_state_out_primary_rarefied20"]={
        "n_deleted_US_state_source_ids":len(leave),
        "n_estimable_withheld_state_refits":len(valid),
        "full_primary_rho":main,
        "min_rho":min(valid) if valid else None,
        "max_rho":max(valid) if valid else None,
        "n_positive":sum(v>0 for v in valid),
        "n_negative":sum(v<0 for v in valid),
        "NOT_49_independent_country_replicates":True,
    }
    return out


def run(fa,botanical,census_response):
    botanical_eligible,excluded=validated_botanical_states(botanical)
    geo=parse_census_centroids(census_response.text,botanical["states"])
    pairs=state_pair_rows(botanical_eligible,geo,fa)
    if len(pairs)!=len(botanical_eligible)*(len(botanical_eligible)-1)//2:
        raise ValueError("One of the fixed eligible US state source pairs missing")
    return {
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "original_FA6_2022_OD_source_sha256_from_v80":fa.get("source_hash_verified_vs_prior_v79"),
        "Census_official_2025_TIGER_state_table_sha256":sha256(census_response.content).hexdigest(),
        "original_GBIF_indexed_flora_period":"2018_2025",
        "FAF6_freight_estimate_base_year":2022,
        "Census_state_FIPS_geo_universe":51,
        "n_eligible_country_GRIIS_indexed_US_state_sources":len(botanical_eligible),
        "noneligible_state_codes_and_GBIF_status_reasons":excluded,
        "n_physical_or_economic_OD_state_pairs_NOT_independent":len(pairs),
        "trait_and_introduction_status_NOT_verified_for_individual_state":True,
        "analysis":analysis(pairs),
        "NO_freight_weed_seed_causality_mowing_wind_adaptation":True,
        "limitations":[
            "Cross-sectional 2022 estimated FAF6 commodity OD vs 2018–2025 GBIF recorded plants: NO before/after seed transport chronology or local establishment observed.",
            "Country GRIIS introduction status is not STATE-specific; 2026 GBIF stateProvince indexed search is a text filter and not verified state polygon geolocation.",
            "2022 economic OD state freight weights are NOT actual track routes, viable transported weed seeds or independent freight delivery experiments.",
            "More similarity in state pairs can arise from geography, climate, flora source datasets, land use, agriculture, urbanization and nonrandom botanical recording.",
            "FRA crossing train counts and state-level freight OD weights must not be interchanged; the source flow has no measured mowing, slipstream wind, disturbance or adaptive evolution.",
            "Pairs reuse states and cannot be assigned 1000+ independent sample p-values. State jackknife is sensitivity, not a set of independent countries.",
            "Raw Jaccard and K20 rarefaction have different source effort sensitivity; no favorable metric selection should be treated as universal rail effect.",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    for key in ("freight-v80","flora-v82","output"):
        ap.add_argument("--"+key,required=True)
    a=ap.parse_args();read=lambda name:json.loads(Path(name).read_text(encoding="utf-8"))
    r=requests.get(CENSUS_TIGER_URL,timeout=(20,90),
                    headers={"User-Agent":"tetsu-2025-Census-state-geography-FAF-GRIIS-v83"})
    r.raise_for_status()
    output=run(read(a.freight_v80),read(a.flora_v82),r)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(output,indent=2,ensure_ascii=False),encoding="utf-8")
    a=output["analysis"]
    print(json.dumps({"complete":True,
       "n_source_states":output["n_eligible_country_GRIIS_indexed_US_state_sources"],
       "n_state_pairs":output["n_physical_or_economic_OD_state_pairs_NOT_independent"],
       "cohorts":{k:v for k,v in a.items() if k in COHORTS},
       "leave_one_state_out":a["leave_one_state_out_primary_rarefied20"],
       "NO_causal_seed_dispersal_or_mowing_wind_evolution":True},
       ensure_ascii=False),flush=True)


if __name__=="__main__":main()
