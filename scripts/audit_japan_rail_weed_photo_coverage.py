#!/usr/bin/env python3
"""Real Japanese herbaceous rail-near photo source-coverage audit.

Use historical Japanese railway weeds as the *predeclared taxa*, but obtain
modern independently sourced iNaturalist 2018-2025 photos. Source counts
are not abundance and railway proximity is not confirmed rail verge habitat.

Queries ALL photo observations by taxon/area when <= 2000 records, rather
than cherry-picking first-flowers or the nicest-looking plants.
No raw coordinates, observation IDs, or third-party photos are exported.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import math
from pathlib import Path

from audit_inat_rail_phenophases import (
    AREAS, N02, URL, RateLimitedClient, acceptable_observation,
    rail_index_for_bbox
)

# Taxa fixed from Abe et al. (2012) top 20, plus a modern accepted
# scientific synonym for Conyza. Grouped Setaria is a GENUS query.
TAXA = {
    "Equisetum arvense": "perennial",
    "Solidago altissima": "perennial",
    "Setaria": "annual_genus",
    "Miscanthus sinensis": "perennial",
    "Imperata cylindrica": "perennial",
    "Digitaria ciliaris": "annual",
    "Pueraria lobata": "perennial_vine",
    "Artemisia princeps": "perennial",
    "Erigeron canadensis": "annual_synonym_for_Conyza_canadensis",
}
YEARS = ("2018-01-01", "2025-12-31")
MAX_COMPLETE_RESULTS = 2000


def classify_rail_distance(distance_m: float) -> str:
    if not math.isfinite(distance_m) or distance_m < 0:
        raise ValueError("invalid rail distance")
    if distance_m <= 100:
        return "0-100m"
    if distance_m <= 250:
        return "100-250m"
    if distance_m <= 1000:
        return "250-1000m"
    return "gt_1000m"


# Published 2012 name Pueraria lobata is an accepted infraspecific taxon
# under Pueraria montana var. lobata according to Kew POWO (2026).
# Never automatically accept undifferentiated Pueraria montana species: it
# includes multiple varieties and would mix distinct source taxon concepts.
CURATED_ALIASES = {
    "Pueraria lobata": {
        "Pueraria lobata",
        "Pueraria montana var. lobata",
        "Pueraria montana var lobata",
    }
}


def matches_focal_name(ob, requested: str):
    raw=(ob.get("taxon") or {}).get("name") or ""
    if not raw:
        return False
    if requested == "Setaria":
        return raw == "Setaria" or raw.startswith("Setaria ")
    if requested in CURATED_ALIASES:
        return raw in CURATED_ALIASES[requested]
    return raw == requested or raw.startswith(requested + " ")


def queried_observations(client, area, taxon):
    swlat,swlng,nelat,nelng = AREAS[area]
    params={
        "taxon_name":taxon, "quality_grade":"research", "photos":"true",
        "d1":YEARS[0], "d2":YEARS[1],
        "swlat":swlat,"swlng":swlng,"nelat":nelat,"nelng":nelng,
        "per_page":200,"order_by":"id","order":"asc"
    }
    first=client.get(URL,params={**params,"page":1}).json()
    total=int(first.get("total_results") or 0)
    if total>MAX_COMPLETE_RESULTS:
        return total, None
    obs=list(first.get("results") or [])
    for p in range(2, math.ceil(total/200)+1):
        obs.extend(client.get(URL,params={**params,"page":p}).json().get("results") or [])
    if len(obs)!=total or len({x["id"] for x in obs})!=len(obs):
        raise RuntimeError("Paginaton drift or duplicate observation IDs")
    return total, obs


def summarize_taxon_observations(all_records, requested, dist_func):
    quality=collections.Counter()
    returned_unmatched_taxa=collections.Counter()
    rail_bins=collections.Counter()
    observers=collections.defaultdict(set)
    phases=collections.Counter()
    for ob in all_records:
        if not matches_focal_name(ob,requested):
            quality["returned_taxon_does_not_match_exact_group"]+=1
            returned_unmatched_taxa[str((ob.get("taxon") or {}).get("name") or "missing_name")]+=1
            continue
        valid,reason=acceptable_observation(ob)
        quality[reason]+=1
        if not valid:
            continue
        lng,lat=ob["geojson"]["coordinates"]
        category=classify_rail_distance(dist_func(lng,lat))
        rail_bins[category]+=1
        who=(ob.get("user") or {}).get("id")
        if who is not None:
            observers[category].add(str(who))
        seen={
            x.get("controlled_value_id") for x in (ob.get("annotations") or [])
            if x.get("controlled_attribute_id")==12
        }
        if 13 in seen and 21 in seen:
            phases["contradictory_13_and_21"]+=1
        elif 13 in seen:
            phases["visible_open_flower_annotation"]+=1
        elif 21 in seen:
            phases["no_visible_flower_or_fruit_annotation"]+=1
        else:
            phases["phase_not_annotated"]+=1
    return {
        "quality_gate":dict(quality),
        "unmatched_api_accepted_taxa":dict(returned_unmatched_taxa.most_common(15)),
        "n_valid_geo":sum(rail_bins.values()),
        "rail_bins":dict(rail_bins),
        "distinct_observers_per_rail_bin":{
            k:len(v) for k,v in observers.items()
        },
        "flower_visibility_annotation_counts":dict(phases),
    }


def run():
    client=RateLimitedClient(delay_seconds=1.1)
    raw=client.get(N02,timeout=(20,150)).content
    n02_sha=hashlib.sha256(raw).hexdigest()
    profiles=[]
    local_lines={}
    distance_functions={}
    for area in AREAS:
        function,n_lines=rail_index_for_bbox(raw,AREAS[area])
        distance_functions[area]=function
        local_lines[area]=n_lines

    for area in ("tokyo_bbox","osaka_bbox","sendai_bbox"):
        for taxon,life in TAXA.items():
            total,observations=queried_observations(client,area,taxon)
            base={
                "area":area,"focal_name":taxon,
                "source_life_history_for_planning":life,
                "api_total":total,
                "complete_source_census":observations is not None,
                "incomplete_source_reason":"more_than_2000_per_taxon_area" if observations is None else None
            }
            if observations is not None:
                base.update(summarize_taxon_observations(
                    observations,taxon,distance_functions[area]
                ))
            profiles.append(base)
            print(json.dumps({
                "area":area,"taxon":taxon,"n_total":total,
                "quality":base.get("n_valid_geo"),
                "rail_le_250":sum(base.get("rail_bins",{}).get(x,0)
                                   for x in ("0-100m","100-250m"))
            },ensure_ascii=False),flush=True)
    return {
        "created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "type":"modern citizen-photo data-source feasibility, not railway habitat or mowing effects",
        "source":URL,"rail_gis_source":N02,
        "n02_sha256":n02_sha,"rail_line_parts_indexed":local_lines,
        "study_area_definitions":AREAS,
        "photo_observation_years":YEARS,
        "focal_taxon_group_labels":TAXA,
        "quality_gate":"iNaturalist research grade, photos, not cultivated, unobscured, explicit 1-100m spatial accuracy, matched taxon",
        "profiles":profiles,
        "warnings":[
            "Source taxa predeclared from 2011 railway-weed survey; historical table counts cannot be quantitatively compared as a temporal trend with 2018-25 citizen-photo records.",
            "iNaturalist reported plants are opportunistically photographed and do not represent abundance, verified local absence, independent complete surveys, or chance of seed transport.",
            "Photographic rail proximity ≤250m does NOT establish growth on rail right-of-way, herbicide/spraying/cutting exposure or train-induced gust.",
            "Flower annotation absence is unknown, not biological absence. No flowering-onset dates or ROI measurements here.",
            "Genus Setaria observations cannot be interpreted as one biological species; accepted taxon names can change.",
            "Pueraria lobata source taxon aliases are restricted to the accepted Pueraria montana var. lobata; Pueraria montana species without variety is not automatically accepted.",
            "Taxon aliases are documented by Kew POWO: https://powo.science.kew.org/taxon/214449-2",
            "Source includes no real management records; environmental filtering is a hypothesis and historical descriptive comparison only.",
            "No third-party images, occurrence IDs, or precise observation coordinates are exported."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/japan_rail_herb_photo_coverage.json")
    a=p.parse_args()
    summary=run()
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "n_taxon_area_queries":len(summary["profiles"]),
        "quality_geo_total":sum(x.get("n_valid_geo",0) for x in summary["profiles"]),
        "rail_le_250_total":sum(sum(x.get("rail_bins",{}).get(k,0) for k in
            ("0-100m","100-250m")) for x in summary["profiles"]),
        "n_noncomplete_queries":sum(not x["complete_source_census"] for x in summary["profiles"]),
    },ensure_ascii=False))
    if not any(x.get("n_valid_geo",0)>0 for x in summary["profiles"]):
        raise SystemExit("No real Japanese georeferenced herbaceous photo data; no ecological analysis possible")


if __name__=="__main__":
    main()
