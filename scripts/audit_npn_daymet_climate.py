#!/usr/bin/env python3
"""Outcome-blind Daymet climate audit for observed NPN railway-near/far sites.

The same taxa, years, site-distance definitions, and 75-km far-site matches
are frozen by the existing NPN data-opportunity pipeline. Daymet 1-km daily
temperature and precipitation are then fetched; no flowering-date outcome is
used in selecting a weather station, cutoff day or thermal threshold.

Outputs aggregate site-equal-weight weather medians, not station locations,
onset effects, mowing measurements or train slipstream.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import statistics
import time
from pathlib import Path

import requests

from audit_npn_candidate_flowering_status import (
    site_frame, site_phenometric_species, select_candidate, match_far
)
from audit_npn_active_rail_sites import validate_station
from audit_npn_registered_sites_rail import registered_sites, STATES
from audit_npn_site_phenometrics_rail import get_site_phenometrics, YEARS

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.daymet import parse_daymet, summarize_annual_climate, FIXED_DAYS

DAYMET = "https://daymet.ornl.gov/single-pixel/api/data"
VERSION = "tetsu-daymet-v0.9"
REQ = requests.Session()
REQ.headers.update({"User-Agent": "tetsu-rail-phenology-climate-audit/0.9 (https://github.com/zuizui0223/tetsu)"})
# Days 60 (MD Cercis before observed spring flowering) and 90 (AZ
# Prosopis before observed first flowers) were already fixed as climate
# cutoffs; expose them, rather than relying only on DOY120 which can
# occur AFTER first flowering for both species.
COLUMNS = (
    "gdd5_to_doy60", "gdd5_to_doy90", "gdd5_to_doy120",
    "mean_air_c_to_doy60", "mean_air_c_to_doy90",
    "mean_air_c_to_doy120",
    "precip_mm_to_doy60", "precip_mm_to_doy90",
    "precip_mm_to_doy120"
)


def collect_frames_and_cases():
    """Freeze site IDs by existing data availability and geography only."""
    frames={}
    coordinates={}
    errors={}
    for state in STATES:
        try:
            frames[state],_ = site_frame(state)
            coordinates[state]={}
            for row in registered_sites(state):
                parsed=validate_station(row)
                if parsed:
                    site, lon, lat=parsed
                    if site in frames[state]:
                        coordinates[state][site]=(lat,lon)
        except Exception as ex:
            errors[state]=f"{type(ex).__name__}: {str(ex)[:160]}"
    cases=[]
    for year in YEARS:
        try:
            rows,digest,_ = get_site_phenometrics(year)
        except Exception as ex:
            errors[f"source-{year}"]=f"{type(ex).__name__}: {str(ex)[:160]}"
            continue
        for state,frame in frames.items():
            species,_ = site_phenometric_species(rows,frame,state)
            top=select_candidate(species)
            if not top:
                continue
            species_id, _, _ = top[0]
            near=set(species[species_id]["near"])
            matches=match_far(near,species[species_id]["far"],frame)
            far={f for n,f,_ in matches}
            cases.append({
                "state":state, "year":year, "species_id":species_id,
                "source_sha256":digest,
                "_near":near, "_far":far,
                "site_count_near":len(near),
                "site_count_far_geo_matched":len(far),
                "n_pairs_geo":len(matches)
            })
    return cases, coordinates, errors


def download_daymet(lat, lon):
    args={
        "lat":f"{lat:.5f}","lon":f"{lon:.5f}",
        "vars":"tmax,tmin,prcp",
        "years":",".join(str(y) for y in YEARS)
    }
    for attempt in range(3):
        try:
            r=REQ.get(DAYMET,params=args,timeout=(20,120))
            if r.status_code in (429,500,502,503,504):
                time.sleep(2+3*attempt)
                continue
            r.raise_for_status()
            if len(r.content)>1_000_000:
                raise ValueError("unexpectedly large single-pixel Daymet response")
            series=parse_daymet(r.text,YEARS)
            result={
                year:summarize_annual_climate(series[year]) for year in YEARS
            }
            return result,hashlib.sha256(r.content).hexdigest()
        except (requests.RequestException,ValueError) as ex:
            if attempt==2: raise
            time.sleep(2+attempt)
    raise RuntimeError("could not fetch Daymet after retries")


def median_or_none(values):
    return round(statistics.median(values),3) if values else None


def run():
    cases, coordinates, frame_errors=collect_frames_and_cases()
    needed=collections.defaultdict(set)
    for case in cases:
        needed[case["state"]].update(case["_near"] | case["_far"])
    weather={}
    source_hashes=[]
    failures=[]
    for state in STATES:
        for sid in sorted(needed[state]):
            coords=coordinates.get(state,{}).get(sid)
            if not coords:
                failures.append(f"coordinate not verified state {state}")
                continue
            try:
                measures,digest=download_daymet(*coords)
                weather[(state,sid)]=measures
                source_hashes.append(digest)
                # Courtesy pause; Daymet 1-km spatial data are server-generated.
                time.sleep(0.35)
            except Exception as ex:
                failures.append(
                    f"Daymet {state}: {type(ex).__name__}: {str(ex)[:145]}"
                )
    outputs=[]
    for case in cases:
        near=case["_near"];far=case["_far"]
        state=case["state"];year=case["year"]
        near_ok=[weather[(state,x)][year] for x in near if (state,x) in weather]
        far_ok=[weather[(state,x)][year] for x in far if (state,x) in weather]
        entry={
            "state":state, "year":year, "NPN_species_id":case["species_id"],
            "n_selected_near_sites":case["site_count_near"],
            "n_selected_far_sites_within_75km":case["site_count_far_geo_matched"],
            "n_near_with_usable_daymet":len(near_ok),
            "n_far_with_usable_daymet":len(far_ok),
            "n_unique_geographic_links":case["n_pairs_geo"],
            "source_phenometric_sha256":case["source_sha256"],
            "fixed_cumulative_day_of_year_cutoffs":list(FIXED_DAYS),
            "fixed_gdd_base_c":5.0,
            "variables":{},
        }
        for key in COLUMNS:
            a=median_or_none([v[key] for v in near_ok])
            b=median_or_none([v[key] for v in far_ok])
            entry["variables"][key]={
                "site_median_near":a,"site_median_far":b,
                "near_minus_far":round(a-b,3) if a is not None and b is not None else None
            }
        outputs.append(entry)
        print(json.dumps({
            "case":f"{state}_{year}_spp{case['species_id']}",
            "source_eligible_near":len(near_ok),"source_eligible_far":len(far_ok),
            "GDD5_D60_near_minus_far":entry["variables"]["gdd5_to_doy60"]["near_minus_far"],
            "GDD5_D90_near_minus_far":entry["variables"]["gdd5_to_doy90"]["near_minus_far"],
            "GDD5_D120_near_minus_far":entry["variables"]["gdd5_to_doy120"]["near_minus_far"],
            "mean_air_D120_near_minus_far":entry["variables"]["mean_air_c_to_doy120"]["near_minus_far"]
        },ensure_ascii=False),flush=True)
    return {
        "run_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "type":"fixed-season temperature confounding audit, NOT railway phenology effect",
        "Daymet_source":DAYMET,
        "Daymet_DOI":"10.3334/ORNLDAAC/2361",
        "Daymet_pixel_resolution":"about 1 km; no microclimate, ballast, wind or cutting exposure",
        "source_requests_success":len(weather),
        "source_requests_expected":sum(len(v) for v in needed.values()),
        "source_checksum_manifest_sha256":hashlib.sha256(
            "".join(sorted(source_hashes)).encode()
        ).hexdigest(),
        "selection_errors":frame_errors,
        "daymet_source_failures":failures[:12],
        "results":outputs,
        "interpretation_limits":[
            "Site choices from NPN positive phenometric opportunity and geography; not selected by flowering date or effect direction.",
            "Only 2 independent rail-near sites per selected taxon-year; this is not an inferential comparison.",
            "No actual railway-vegetation habitat, mowing dates, wind measurements or establishment data.",
            "Daymet gridded macroclimate cannot represent track ballast temperatures or slipstream-induced gusts.",
            "Thermal cutoff days 60/90/120/150 and GDD base 5 C fixed before querying weather.",
            "No genetic adaptation, phenotypic plasticity or railway-related flowering shift is inferred."
        ]
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="results/audit/npn_daymet_gdd_control_v09.json")
    args=parser.parse_args()
    result=run()
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    if result["source_requests_success"]==0:
        raise SystemExit("Daymet unavailable; do not claim temperature-adjusted phenology")


if __name__=="__main__":
    main()
