#!/usr/bin/env python3
"""NPN active-site/same-taxon audit near vs far from FRA railway.

Selection is independent of species and flowering status: choose a fixed
SHA256-ranked set of registered near-rail stations, then choose the nearest
unique far-rail station by projected coordinate within the same state. The
phenology service is queried ONLY AFTER these pairs are frozen.

Output is aggregate eligibility, never individual/site coordinates or IDs.
Not a railway effect estimator. Not a mowing/wind exposure measure.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import math
from pathlib import Path

import requests
from pyproj import Transformer

from audit_npn_registered_sites_rail import (
    STATES, STATIONS_URL, registered_sites, pick_field
)
from audit_npn_rail_proximity import (
    get_rail_features, UTM_EPSG, SESSION
)
from audit_usa_npn_status import ENDPOINT, parse_stream
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.npn_individuals import normalize_open_flowers, summarize_individuals


YEAR = 2024
MAX_PAIR_KM = 75.
SITE_MAX_ROWS = 24000
SITE_MAX_RESPONSE_BYTES = 16_000_000
CLIENT = requests.Session()
CLIENT.headers.update({
    "User-Agent": "tetsu-phenology-active-rail-site-check/0.7 (https://github.com/zuizui0223/tetsu)"
})


def validate_station(row):
    sid = pick_field(row, ("station_id", "site_id", "stationId", "id"))
    lat = pick_field(row, ("latitude", "lat", "station_latitude"))
    lon = pick_field(row, ("longitude", "long", "station_longitude"))
    try:
        sid = int(sid)
        lat, lon = float(lat), float(lon)
    except (ValueError, TypeError):
        return None
    if sid < 1 or not (math.isfinite(lat) and math.isfinite(lon)):
        return None
    if not (-125 <= lon <= -65 and 24 <= lat <= 50):
        return None
    return sid, lon, lat


def freeze_sites(state: str, n_pairs: int):
    """Freeze near 0-250m / far >=1000m coordinates before any outcome.

    Unique far control is picked by shortest in-state projected distance,
    breaking ties by integer ID, and must lie <=75 km away. A near site
    without an eligible matched far is logged as excluded, not replaced
    according to a floral outcome.
    """
    rows = registered_sites(state)
    closest, fra_info = get_rail_features(state)
    t = Transformer.from_crs(4326, UTM_EPSG[state], always_xy=True)
    valid, skipped = {}, collections.Counter()
    for row in rows:
        candidate = validate_station(row)
        if not candidate:
            skipped["bad_id_or_location"] += 1
            continue
        site_id, lon, lat = candidate
        if site_id in valid:
            skipped["duplicate_station_id"] += 1
            continue
        rail_dist = closest(lon, lat)
        xy = t.transform(lon, lat)
        valid[site_id] = (rail_dist, xy, lon, lat)
    near = [sid for sid, (distance, *_rest) in valid.items() if distance <= 250]
    far = [sid for sid, (distance, *_rest) in valid.items() if distance >= 1000]
    near.sort(key=lambda x: (hashlib.sha256(
        f"tetsu-v0.7-registered-site-sampling-{state}-{YEAR}-{x}".encode()
    ).hexdigest(), x))
    used_far = set()
    pairs = []
    for sid in near:
        if len(pairs) >= n_pairs:
            break
        a = valid[sid][1]
        candidate_far = [
            (math.hypot(a[0] - valid[j][1][0], a[1] - valid[j][1][1]), j)
            for j in far if j not in used_far
        ]
        if not candidate_far:
            skipped["no_far_registered_station"] += 1
            continue
        distance, matched = min(candidate_far)
        if distance > MAX_PAIR_KM * 1000:
            skipped["no_far_within_75km"] += 1
            continue
        used_far.add(matched)
        pairs.append((sid, matched, distance))
    summary = {
        "n_registered_api_rows": len(rows),
        "n_valid_station_ids": len(valid),
        "n_registered_within_250m": len(near),
        "n_registered_beyond_1000m": len(far),
        "n_pairs_selected_without_outcomes": len(pairs),
        "n_near_sites_selected": len({x[0] for x in pairs}),
        "n_far_sites_selected": len({x[1] for x in pairs}),
        "pair_distance_range_km": (
            [round(min(x[2] for x in pairs)/1000, 2),
             round(max(x[2] for x in pairs)/1000, 2)] if pairs else None
        ),
        "skips": dict(skipped),
        "rail_source_vintage": "FRA 2026",
        "rail_source_feature_count": fra_info["n_fra_features"],
    }
    return pairs, summary


def download_station_records(site_id, state, year=YEAR):
    """Official NPN station_id filter. Never silently accept an ignored filter."""
    if not isinstance(site_id, int) or site_id < 1 or year != YEAR:
        raise ValueError("invalid site/year")
    form = {
        "request_src": "tetsu_research_feasibility",
        "climate_data": "0",
        "station_id[1]": str(site_id),
        "state[1]": state,
        "start_date": f"{year}-01-01",
        "end_date": f"{year}-12-31",
    }
    with CLIENT.post(
        ENDPOINT, data=form, stream=True, timeout=(18, 140)
    ) as response:
        response.raise_for_status()
        records, sha256, bytes_read = parse_stream(response)
    if len(records) > SITE_MAX_ROWS or bytes_read > SITE_MAX_RESPONSE_BYTES:
        raise ValueError("site response exceeds per-site safety cap")
    if records:
        observed_ids = {str(r.get("site_id")) for r in records}
        if observed_ids != {str(site_id)}:
            raise ValueError(
                "NPN site filter failed: unexpected ID(s), refusing to count"
            )
        dates = {str(r.get("observation_date", ""))[:4] for r in records}
        if dates != {str(year)}:
            raise ValueError(
                "NPN year filter failed: unexpected dates, refusing to count"
            )
    return records, sha256, bytes_read


def reduce_station(rows):
    """Return site-level per-species active/eligible statistics.

    Organism level: individual_id repeated on distinct dates and at least
    one verified YES and NO for Open flowers. Missing phenophase is unknown.
    Site×species eligibility: >=2 distinct monitored individuals, >=4
    valid observed Open flower visits, both states observed.
    """
    open_rows = [r for r in rows if str(r.get("phenophase_description")) == "Open flowers"]
    by_taxon = collections.defaultdict(list)
    missing_species = 0
    for r in open_rows:
        species = r.get("species_id")
        try:
            key = int(species)
            if key < 1:
                raise ValueError
        except (TypeError, ValueError):
            missing_species += 1
            continue
        by_taxon[key].append(r)
    per_species = {}
    for species_id, items in by_taxon.items():
        visits, _ = normalize_open_flowers(items)
        individuals = summarize_individuals(visits)
        observed_yes = sum(i["n_yes"] for i in individuals)
        observed_no = sum(i["n_no"] for i in individuals)
        eligible_individuals = sum(
            i["n_dates"] >= 4 and i["n_yes"] >= 1 and i["n_no"] >= 1
            for i in individuals
        )
        eligible_site_species = (
            len(visits) >= 4 and
            observed_yes > 0 and observed_no > 0 and
            eligible_individuals >= 1
        )
        per_species[species_id] = {
            "n_valid_observations": len(visits),
            "n_individuals": len(individuals),
            "n_repeat_monitored_yes_and_no_individuals": eligible_individuals,
            "has_valid_yes": bool(observed_yes),
            "has_valid_no": bool(observed_no),
            "site_taxon_eligible": bool(eligible_site_species),
        }
    return {
        "n_records": len(rows),
        "n_openflower_raw": len(open_rows),
        "n_openflower_species_id_missing": missing_species,
        "species": per_species,
        "n_species_with_any_observed_openflower": len(by_taxon),
    }


def summarize_cohort(site_results, pairs, state):
    near_active = set()
    far_active = set()
    near_species = collections.Counter()
    far_species = collections.Counter()
    any_active_near, any_active_far = 0, 0
    paired_active = 0
    for near_id, far_id, _ in pairs:
        left = site_results.get(near_id)
        right = site_results.get(far_id)
        if left is None or right is None:
            continue
        if left["n_openflower_raw"]:
            any_active_near += 1
        if right["n_openflower_raw"]:
            any_active_far += 1
        if left["n_openflower_raw"] and right["n_openflower_raw"]:
            paired_active += 1
        for sid, result, counter, used in (
            (near_id, left, near_species, near_active),
            (far_id, right, far_species, far_active)
        ):
            for taxon, stats in result["species"].items():
                if stats["site_taxon_eligible"]:
                    counter[taxon] += 1
                    used.add(sid)
    common = near_species.keys() & far_species.keys()
    # Explicit rank by geographic replication, then taxon id: selection
    # is based on *eligibility*, not on a flowering effect/direction.
    common_rows = sorted([
        {
            "npn_species_id": k,
            "n_near_eligible_sites": near_species[k],
            "n_far_eligible_sites": far_species[k],
        }
        for k in common
    ], key=lambda x: (-min(x["n_near_eligible_sites"], x["n_far_eligible_sites"]),
                      -x["n_near_eligible_sites"]-x["n_far_eligible_sites"],
                      x["npn_species_id"]))
    return {
        "state": state,
        "n_complete_station_pairs_queried": sum(
            (left in site_results and right in site_results) for left, right, _ in pairs
        ),
        "n_near_sites_with_any_openflower_phenophase": any_active_near,
        "n_far_sites_with_any_openflower_phenophase": any_active_far,
        "n_pairs_with_both_sides_any_openflower": paired_active,
        "n_near_sites_with_eligible_species": len(near_active),
        "n_far_sites_with_eligible_species": len(far_active),
        "n_eligible_species_near": len(near_species),
        "n_eligible_species_far": len(far_species),
        "n_shared_eligible_species_near_far": len(common_rows),
        "shared_eligible_species_without_outcome_direction": common_rows[:25],
        "criteria": "≥4 distinct Open flower status dates on ≥1 individual with ≥1 YES and ≥1 NO; counted per site×species; requires both near and far sites to share species",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs-per-state", type=int, default=12)
    p.add_argument("--output", default="results/audit/npn_rail_active_sites_2024.json")
    args = p.parse_args()
    if not (1 <= args.pairs_per_state <= 50):
        p.error("pairs-per-state must be 1..50")
    result = {
        "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "design": "outcome-blind site-first geographic matching, then USA-NPN single-site 2024 records",
        "year": YEAR, "pairs_per_state_requested": args.pairs_per_state,
        "states": {}, "errors": {}, "source": {
            "rail": "FRA NARN (2026 geometry)",
            "NPN_station_list": STATIONS_URL,
            "NPN_status": ENDPOINT,
        },
    }
    for state in STATES:
        try:
            pairs, source = freeze_sites(state, args.pairs_per_state)
            site_results = {}
            errors = []
            for pair in pairs:
                for site_id in pair[:2]:
                    if site_id in site_results:
                        continue
                    try:
                        rec, _digest, _bytes = download_station_records(site_id, state)
                        site_results[site_id] = reduce_station(rec)
                    except Exception as exc:
                        errors.append(f"{type(exc).__name__}: {str(exc)[:155]}")
            report = summarize_cohort(site_results, pairs, state)
            report["site_frame"] = source
            report["n_station_requests_failed"] = len(errors)
            report["first_three_station_request_errors"] = errors[:3]
            result["states"][state] = report
            print(json.dumps({
                "state": state,
                "selected_pairs": len(pairs),
                "station_responses": len(site_results),
                "failed_requests": len(errors),
                "active_near": report["n_near_sites_with_any_openflower_phenophase"],
                "active_far": report["n_far_sites_with_any_openflower_phenophase"],
                "shared_eligible_species": report["n_shared_eligible_species_near_far"],
                "top_species": report["shared_eligible_species_without_outcome_direction"][:6],
            }, ensure_ascii=False), flush=True)
        except Exception as exc:
            result["errors"][state] = f"{type(exc).__name__}: {str(exc)[:320]}"
    result["limitations"] = [
        "Small predetermined station sample, not a census of all 170 registered rail-near sites.",
        "FRA 2026 railway geometry is not confirmed active line service in observation year 2024.",
        "Station registered proximity does not mean monitored individuals grow on railway land, experience train gusts or cutting.",
        "All flowering effects, absolute phenology shifts and causal inferences are deliberately omitted.",
        "Active site is based on actual downloaded 2024 records, not an NPN station registration or NPN species list.",
        "Matched geographic site pairs do not control urbanization, road verge, observer effort, land cover or climate.",
        "Selection based on independent site coordinates, not flowering direction. Never interpret null support as no biological railway effect.",
    ]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    if not result["states"]:
        raise SystemExit("No state source was fully audited. No inference.")
    # Leave nonzero when no complete pair has an API response on both sides.
    if not any(v["n_complete_station_pairs_queried"] > 0 for v in result["states"].values()):
        raise SystemExit("No complete site pair fetched; no evidence of active-site eligibility.")


if __name__ == "__main__":
    main()
