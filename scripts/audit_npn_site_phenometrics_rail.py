#!/usr/bin/env python3
"""Observed site phenometrics × railway proximity screening, 2023–2025.

USA-NPN Site Phenometrics are derived site×species annual first/last yes
summaries, NOT true flowering onset and NOT environmental railway effects.
Screen all species using phenophase_id=501 (verified as Open flowers in actual
Liriodendron and Carnegiea 2024 status records) *before* selecting any taxa.

No individual IDs, site IDs, coordinates, or absolute onset dates exported.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
from pathlib import Path

import requests

from audit_npn_registered_sites_rail import STATES, registered_sites
from audit_npn_rail_proximity import get_rail_features
from audit_npn_active_rail_sites import validate_station
from audit_usa_npn_status import parse_stream

NPN_SITE_API="https://services.usanpn.org/npn_portal/observations/getSiteLevelData.ndjson"
PHASE_ID=501
YEARS=(2023,2024,2025)


def get_site_phenometrics(state,year):
    if state not in STATES or year not in YEARS:
        raise ValueError("state/year not in predeclared screen")
    form={
        "request_src":"tetsu_rail_repeated_phenology_feasibility",
        "climate_data":"0",
        "state[1]":state,
        "phenophase_id[1]":str(PHASE_ID),
        "start_date":f"{year}-01-01",
        "end_date":f"{year}-12-31",
        "num_days_quality_filter":"30",
    }
    with requests.post(
        NPN_SITE_API, data=form, stream=True, timeout=(20,200),
        headers={"User-Agent":"tetsu-USANPN-openflower-rail-eligibility/0.8"}
    ) as r:
        r.raise_for_status()
        rows, sha, nbytes = parse_stream(r)
    if rows and not any("site_id" in r or "station_id" in r for r in rows[:3]):
        raise ValueError(f"NPN site-level source omitted ID: {sorted(rows[0])}")
    # A year / state / phenophase filter silently ignored is worse than an error.
    unexpected_phase=[
        r.get("phenophase_description") for r in rows
        if r.get("phenophase_description") not in (None,"Open flowers")
    ]
    if unexpected_phase:
        raise ValueError(
            f"NPN site phenophase filter failed: {unexpected_phase[:4]}"
        )
    unexpected_year=[
        r.get("year") for r in rows
        if r.get("year") not in (None,"",year,str(year))
    ]
    if unexpected_year:
        raise ValueError(f"NPN year filter failed: {unexpected_year[:3]}")
    return rows,sha,nbytes


def build_registered_rail_frame(state):
    near_rail, fra_info=get_rail_features(state)
    site_rail_bin={}
    invalid=collections.Counter()
    for rec in registered_sites(state):
        x=validate_station(rec)
        if x is None:
            invalid["missing_or_bad_id_or_coords"]+=1
            continue
        site,lon,lat=x
        if site in site_rail_bin:
            invalid["duplicate_site"]+=1
            continue
        distance=near_rail(lon,lat)
        site_rail_bin[site]=(
            "near" if distance<=250 else
            "far" if distance>=1000 else
            "middle_excluded"
        )
    return site_rail_bin,{"rail_n_features":fra_info["n_fra_features"],
                          "n_near":sum(x=="near" for x in site_rail_bin.values()),
                          "n_far":sum(x=="far" for x in site_rail_bin.values()),
                          "invalid":dict(invalid)}


def summarize_state_year(rows,rail_bins,state,year):
    # site×species rows represent observed first YES (possibly averaged
    # across monitored individuals); count each unique site/species once.
    quality=collections.Counter()
    sites=collections.defaultdict(set)
    species=collections.defaultdict(lambda:collections.defaultdict(set))
    site_rows=set()
    first_fields=set()
    for row in rows:
        try:
            site_value=row.get("site_id",row.get("station_id"))
            site_id=int(site_value)
            species_id=int(row.get("species_id"))
            if species_id<=0: raise ValueError()
        except (TypeError, ValueError):
            quality["invalid_site_or_species_id"]+=1
            continue
        b=rail_bins.get(site_id)
        if b is None:
            quality["site_not_in_registered_state_frame"]+=1
            continue
        if b=="middle_excluded":
            quality["intermediate_rail_distance"]+=1
            continue
        # The existence of a site phenometric row means that the phenophase
        # occurred in the source record. It does not establish systematic
        # NO observation coverage beyond the API quality filter.
        sites[b].add(site_id)
        species[species_id][b].add(site_id)
        k=(site_id,species_id)
        if k in site_rows:
            quality["duplicate_site_species_rows"]+=1
        site_rows.add(k)
    common=[]
    for taxon,pairs in species.items():
        a=len(pairs.get("near",set()))
        b=len(pairs.get("far",set()))
        if a and b:
            common.append({
                "npn_species_id":taxon,
                "n_near_sites_with_positive_site_phenometric":a,
                "n_far_sites_with_positive_site_phenometric":b,
            })
    common.sort(key=lambda r:(-min(
        r["n_near_sites_with_positive_site_phenometric"],
        r["n_far_sites_with_positive_site_phenometric"]
    ),-(r["n_near_sites_with_positive_site_phenometric"]+
        r["n_far_sites_with_positive_site_phenometric"]),
        r["npn_species_id"]))
    return {
        "state":state,"year":year,
        "raw_source_site_phenometric_rows":len(rows),
        "n_registered_sites_with_positive_Open_flowers_near":len(sites["near"]),
        "n_registered_sites_with_positive_Open_flowers_far":len(sites["far"]),
        "n_species_with_any_positive_site_phenometric":len(species),
        "n_species_with_positive_phenometric_both_near_far":len(common),
        "species_eligible_for_later_individual_status_lookups":common[:30],
        "quality":dict(quality),
        "definition":"Presence of an annual site-level Open flowers phenometric; MUST independently check actual 0/1 status visits at selected sites to estimate onset intervals.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/npn_site_phenometrics_rail_2023_25.json")
    a=p.parse_args()
    audit={
        "run_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "api":NPN_SITE_API,"years":YEARS,"states":{},
        "phase_filter":"phenophase_id[1]=501 (Open flowers)",
        "errors":{},
        "notes":[
            "Site Phenometrics selects presence of flowers, not observed biological absence; original status data needed for no-to-yes interval tests.",
            "Annual site first/last YES may be an average of multiple plants; it does not identify true first flowering or multiple reflowering bouts.",
            "Rail exposure is station-centroid distance to FRA 2026 map, not known 2023–25 train service, rail-side habitat, mowing or slipstream.",
            "All enrolled taxa were screened without reading their relative onset timing, which is never exported as a railway effect.",
            "No study site IDs, individual IDs or coordinates are released."
        ]
    }
    for state in STATES:
        try:
            rail,source=build_registered_rail_frame(state)
            cases=[]
            for year in YEARS:
                try:
                    rows,digest,nbytes=get_site_phenometrics(state,year)
                    result=summarize_state_year(rows,rail,state,year)
                    result["source_sha256"]=digest
                    result["source_bytes"]=nbytes
                    cases.append(result)
                    print(json.dumps({
                        "state":state,"year":year,
                        "source_rows":len(rows),
                        "sites_near":result["n_registered_sites_with_positive_Open_flowers_near"],
                        "sites_far":result["n_registered_sites_with_positive_Open_flowers_far"],
                        "shared_species":result["n_species_with_positive_phenometric_both_near_far"],
                        "species_top3":result["species_eligible_for_later_individual_status_lookups"][:3],
                    }),flush=True)
                except Exception as ex:
                    audit["errors"][f"{state}-{year}"]=f"{type(ex).__name__}: {str(ex)[:400]}"
            audit["states"][state]={"source_registered":source,"year_results":cases}
        except Exception as ex:
            audit["errors"][state]=f"{type(ex).__name__}: {str(ex)[:400]}"
    out=Path(a.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    n=sum(len(v["year_results"]) for v in audit["states"].values())
    print(json.dumps({"successful_state_year_queries":n,"errors":audit["errors"]}),flush=True)
    if n==0: raise SystemExit("No valid site phenometric API result, no inference")


if __name__=="__main__":
    main()
