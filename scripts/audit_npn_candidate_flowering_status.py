#!/usr/bin/env python3
"""Follow up repeated individual status for outcome-blind common-species sites.

Sources: USA-NPN annual site phenometrics (positive 'Open flowers' at a site),
USA-NPN individual 0/1 status, registered site coordinates, FRA 2026 rails.
Candidates are ranked by number of observed sites in BOTH rail-near/far bins,
never by a flowering-time difference or an effect sign.

Does not estimate train-induced change, mowing effects, or adaptation.
Only aggregate counts and verified taxon names are published.
"""
from __future__ import annotations
import argparse
import collections
import datetime as dt
import json
import math
import statistics
from pathlib import Path

import requests
from pyproj import Transformer

from audit_npn_registered_sites_rail import STATES, registered_sites
from audit_npn_rail_proximity import UTM_EPSG, get_rail_features
from audit_npn_active_rail_sites import validate_station
from audit_npn_site_phenometrics_rail import YEARS, get_site_phenometrics
from audit_usa_npn_status import ENDPOINT, parse_stream
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.npn_individuals import normalize_open_flowers, summarize_individuals

PHASE_ID=501
YEAR_CHOICES=(2023,2024,2025)
MAX_DISTANCE_KM=75
MAX_FAR_PER_NEAR=2
SPECIES_URL="https://services.usanpn.org/npn_portal/species/getSpeciesById.json"

SESSION=requests.Session()
SESSION.headers.update({"User-Agent":"tetsu-rail-phenology-candidate-verification/0.8"})


def site_frame(state):
    nearest, rail_meta=get_rail_features(state)
    utm=Transformer.from_crs(4326, UTM_EPSG[state], always_xy=True)
    sites={}
    for row in registered_sites(state):
        x=validate_station(row)
        if not x:
            continue
        sid,lon,lat=x
        if sid in sites:
            continue
        dist=nearest(lon,lat)
        xy=utm.transform(lon,lat)
        sites[sid]={
            "bin":"near" if dist<=250 else "far" if dist>=1000 else "excluded",
            "xy":xy
        }
    return sites, rail_meta


def site_phenometric_species(rows,sites,state):
    by_taxon=collections.defaultdict(lambda:collections.defaultdict(set))
    seen=collections.Counter()
    for r in rows:
        try:
            site=int(r.get("site_id", r.get("station_id")))
            taxon=int(r["species_id"])
        except (KeyError,TypeError,ValueError):
            seen["missing_ids"]+=1
            continue
        if site not in sites:
            continue
        if str(r.get("state") or "").strip() not in ("",state):
            seen["state_conflict"]+=1
            continue
        bin=sites[site]["bin"]
        if bin in ("near","far"):
            by_taxon[taxon][bin].add(site)
    return by_taxon,dict(seen)


def select_candidate(by_taxon):
    # At least 2 near and 3 far verified positive site-phenometric rows.
    # Pick one per state-year; rank only by observation opportunity.
    eligible=[]
    for taxon,groups in by_taxon.items():
        n=len(groups["near"]);f=len(groups["far"])
        if n>=2 and f>=3:
            eligible.append((taxon,n,f))
    return sorted(eligible,key=lambda r:(-min(r[1],r[2]),-(r[1]+r[2]),r[0]))[:1]


def match_far(near_ids,far_ids,frame):
    """Select unique far active phenometric sites within 75km of near sites.

    Selection is geography-only once the species occurrence opportunity is
    established. No flowering date or state is examined for control choice.
    """
    available=set(far_ids)
    matches=[]
    for sid in sorted(near_ids):
        ax,ay=frame[sid]["xy"]
        options=sorted(
            ((math.hypot(ax-frame[fid]["xy"][0],ay-frame[fid]["xy"][1])/1000, fid)
             for fid in available),
            key=lambda x:(x[0],x[1])
        )
        best=[(dist,fid) for dist,fid in options if dist<=MAX_DISTANCE_KM][:MAX_FAR_PER_NEAR]
        for d,fid in best:
            matches.append((sid,fid,d))
            available.remove(fid)
    return matches


def species_name(taxon):
    try:
        r=SESSION.get(SPECIES_URL,params={"species_id":taxon},timeout=(15,45))
        r.raise_for_status()
        obj=r.json()
        arr=obj if isinstance(obj,list) else [obj]
        found=[x for x in arr if isinstance(x,dict)]
        if len(found)!=1:
            return {"id":taxon,"name":None,"status":"unverified_response"}
        data=found[0]
        # Source schemas may use genus+species, scientific_name, or taxon_name.
        candidate=(data.get("scientific_name") or data.get("species_name") or
                   data.get("taxon_name"))
        if candidate is None and data.get("genus") and data.get("species"):
            candidate=f"{data['genus']} {data['species']}"
        if candidate is None:
            return {"id":taxon,"name":None,"status":"name_field_not_identified",
                    "returned_fields":sorted(data.keys())}
        return {"id":taxon,"name":str(candidate),"status":"api_lookup"}
    except (requests.RequestException,ValueError) as ex:
        return {"id":taxon,"name":None,"status":f"error:{type(ex).__name__}"}


def midpoint_doy_of_observed_onset_bracket(individual_summary):
    """Interval midpoint in day-of-year, NOT known true flowering onset."""
    earlier=individual_summary.get("first_yes_prior_no_date")
    later=individual_summary.get("first_yes_observed_date")
    if not earlier or not later:
        return None
    a=dt.date.fromisoformat(earlier)
    b=dt.date.fromisoformat(later)
    if not a<b or a.year!=b.year or (b-a).days>21:
        return None
    return ((a-dt.date(a.year,1,1)).days+1 +
            (b-dt.date(b.year,1,1)).days+1)/2.0


def status_for_station(site,year,taxon):
    form={
        "request_src":"tetsu_rail_phenology_candidate_verification",
        "climate_data":"0",
        "station_id[1]":str(site),
        "species_id[1]":str(taxon),
        "phenophase_id[1]":str(PHASE_ID),
        "start_date":f"{year}-01-01",
        "end_date":f"{year}-12-31",
    }
    with SESSION.post(ENDPOINT,data=form,stream=True,timeout=(20,140)) as resp:
        resp.raise_for_status()
        rows,_digest,_size=parse_stream(resp)
    if rows:
        if {int(x.get("site_id",-1)) for x in rows}!={int(site)}:
            raise ValueError("Ignored station filter; refuse source")
        if {int(x.get("species_id",-1)) for x in rows}!={int(taxon)}:
            raise ValueError("Ignored taxon filter; refuse source")
        if {str(x.get("observation_date",""))[:4] for x in rows}!={str(year)}:
            raise ValueError("Ignored year filter; refuse source")
        if any(x.get("phenophase_description")!="Open flowers" for x in rows):
            raise ValueError("Ignored Open flowers phenophase filter")
    clean,excluded=normalize_open_flowers(rows)
    individuals=summarize_individuals(clean)
    eligible=[x for x in individuals if x["n_dates"]>=4 and x["n_yes"]>0 and x["n_no"]>0]
    midpoints=[
        midpoint_doy_of_observed_onset_bracket(x) for x in eligible
    ]
    midpoints=[x for x in midpoints if x is not None]
    widths=[x["first_yes_bracket_days"] for x in eligible
            if midpoint_doy_of_observed_onset_bracket(x) is not None]
    return {
        "n_raw":len(rows),"n_verified_visit_dates":len(clean),
        "n_monitored_individuals":len(individuals),
        "n_eligible_individual_years":len(eligible),
        "n_interval_bracketed_onsets_leq21d":len(midpoints),
        "n_monitored_reflower_candidates":sum(x["has_two_no_interrupted_flowering"] for x in eligible),
        "site_median_bracket_midpoint_doy":statistics.median(midpoints) if midpoints else None,
        "site_median_bracket_width_days":statistics.median(widths) if widths else None,
    }


def run():
    report={
        "created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "phase":"Open flowers (NPN phenophase ID 501)",
        "source_NPN_site_phenometrics":"https://services.usanpn.org/npn_portal/observations/getSiteLevelData.ndjson",
        "source_NPN_raw_status":ENDPOINT,
        "states":{},"errors":{},
        "warnings":[
            "Source positive site-phenometrics select taxa by DATA AVAILABILITY, not by phenological effect.",
            "A matched far site is <=75km away within a state but land cover, road exposure, climate/GDD, cutting and actual train traffic are NOT controlled.",
            "Station rail proximity uses 2026 FRA rail network; even a 100m rail-near station is not proof of trackside plant microhabitat.",
            "No phenological effect contrast, test statistic or causal p-value is computed here.",
            "No precise station coordinates, station/individual IDs or dates published."
        ]
    }
    frames={}
    for state in STATES:
        try:
            frames[state],_ = site_frame(state)
            report["states"][state]=[]
        except Exception as ex:
            report["errors"][state]=f"{type(ex).__name__}: {str(ex)[:240]}"
    for year in YEAR_CHOICES:
        try:
            rows, digest, nbytes=get_site_phenometrics(year)
        except Exception as ex:
            report["errors"][str(year)]=f"{type(ex).__name__}: {str(ex)[:240]}"
            continue
        for state,frame in frames.items():
            species,diagnostics=site_phenometric_species(rows,frame,state)
            top=select_candidate(species)
            if not top:
                report["states"][state].append({
                    "year":year,"status":"no_taxon_with_2_near_and_3_far_phenometric_sites",
                    "source_national_sha256":digest
                })
                continue
            taxon,n_near,n_far=top[0]
            near_ids=species[taxon]["near"]
            far_ids=species[taxon]["far"]
            matches=match_far(near_ids,far_ids,frame)
            # Also report geographic overlap before querying individual
            # flowering time. No matching is based on the phenology outcome.
            selected_near=set(near_ids)
            selected_far={z[1] for z in matches}
            response_groups={"near":[],"far":[]}
            errors=[]
            for group,ids in (("near",selected_near),("far",selected_far)):
                for sid in sorted(ids):
                    try:
                        response_groups[group].append(status_for_station(sid,year,taxon))
                    except Exception as ex:
                        errors.append(f"{type(ex).__name__}: {str(ex)[:155]}")
            def reduce_group(values):
                # Site is the independent unit. Do NOT pool individual-day
                # observations to create artificial treatment replication.
                site_doys=[
                    v["site_median_bracket_midpoint_doy"] for v in values
                    if v["site_median_bracket_midpoint_doy"] is not None
                ]
                site_widths=[
                    v["site_median_bracket_width_days"] for v in values
                    if v["site_median_bracket_width_days"] is not None
                ]
                return {
                    "queried_successful_sites":len(values),
                    "n_sites_with_eligible_individual_yes_no":sum(x["n_eligible_individual_years"]>0 for x in values),
                    "n_eligible_individual_years":sum(x["n_eligible_individual_years"] for x in values),
                    "n_bracketed_first_yes":sum(x["n_interval_bracketed_onsets_leq21d"] for x in values),
                    "n_reflower_candidates":sum(x["n_monitored_reflower_candidates"] for x in values),
                    "n_raw_records":sum(x["n_raw"] for x in values),
                    "n_sites_with_onset_brackets":len(site_doys),
                    "site_equal_weight_median_bracket_midpoint_doy":round(statistics.median(site_doys),2) if site_doys else None,
                    "site_equal_weight_median_bracket_width_days":round(statistics.median(site_widths),2) if site_widths else None,
                }
            near_summary=reduce_group(response_groups["near"])
            far_summary=reduce_group(response_groups["far"])
            raw_shift=None
            if (near_summary["site_equal_weight_median_bracket_midpoint_doy"] is not None and
                far_summary["site_equal_weight_median_bracket_midpoint_doy"] is not None):
                raw_shift=round(
                    near_summary["site_equal_weight_median_bracket_midpoint_doy"] -
                    far_summary["site_equal_weight_median_bracket_midpoint_doy"],2
                )
            item={
                "year":year, "state":state,
                "selected_taxon":species_name(taxon),
                "positive_phenometric_near_sites":n_near,
                "positive_phenometric_far_sites_statewide":n_far,
                "far_sites_within_75km_geographically_matched":len(selected_far),
                "n_unique_75km_near_far_geographic_links":len(matches),
                "minimum_matched_site_separation_km":round(min((x[2] for x in matches),default=0),1) if matches else None,
                "near_raw_status_eligibility":near_summary,
                "matched_far_raw_status_eligibility":far_summary,
                "exploratory_near_minus_far_median_interval_midpoint_days":raw_shift,
                "raw_shift_caveat":"Conditioned on sites with positive flowers; only 2 near sites; no controls for site climate/urban/roads/observer frequency; NOT a railway effect.",
                "n_status_query_errors":len(errors),
                "status_query_error_examples":errors[:3],
                "national_site_phenometric_sha256":digest,
                "meaning":"source and geography eligibility only; NOT railway-associated flowering effect"
            }
            report["states"][state].append(item)
            print(json.dumps({
                "state":state,"year":year,"taxon":taxon,
                "phenometric_near":n_near,"far":n_far,
                "far_matched75km":len(selected_far),
                "near_eligible_sites":item["near_raw_status_eligibility"]["n_sites_with_eligible_individual_yes_no"],
                "far_eligible_sites":item["matched_far_raw_status_eligibility"]["n_sites_with_eligible_individual_yes_no"],
                "errors":len(errors)
            },ensure_ascii=False),flush=True)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/npn_candidate_flowering_verification.json")
    args=p.parse_args()
    report=run()
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    if not any(any(x.get("selected_taxon") for x in v) for v in report["states"].values()):
        raise SystemExit("No candidate taxon qualified before individual status; no inference")
    print(json.dumps({"n_valid_candidate_checks":sum(
        len([x for x in v if x.get("selected_taxon")]) for v in report["states"].values()),
        "errors":report["errors"]},ensure_ascii=False))


if __name__=="__main__":
    main()
