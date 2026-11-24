#!/usr/bin/env python3
"""Source-verified Slovenian railway-only BS/BT thermal indicator comparison.

Not a rail-vs-nonrail control, measured microclimate, or climate-offset test.
Input: unchanged publisher XLSX and v102 source audit JSON.
"""
import argparse
import json
import math
import random
import re
from collections import defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_SHA256, extract_original_xlsx, validate_original_keys
)

def value(v):
    if v is None or str(v).strip().casefold() in ("", "x", "na", "nan"):
        return None
    try:
        z = float(v)
    except (ValueError, TypeError):
        return None
    return z if math.isfinite(z) and 1 <= z <= 12 else None

def find_temperature_column(traits):
    columns = sorted(set().union(*(x.keys() for x in traits.values())))
    candidates = []
    for key in columns:
        vv = [t.get(key) for t in traits.values()]
        numeric = sum(value(x) is not None for x in vv)
        xs = sum(str(x).strip().casefold() == "x" for x in vv)
        missing = sum(x is None or str(x).strip() == "" for x in vv)
        if numeric == 142 and xs == 135 and missing == 35:
            candidates.append(key)
    if len(candidates) == 1:
        return candidates[0], "exact 142/135/35 source signature"
    named = [k for k in columns if re.search(r'temperat|(?:^|[_. ])eiv[_. ]?t$|^t$', k, re.I)
             and sum(value(t.get(k)) is not None for t in traits.values()) >= 30]
    if len(named) == 1:
        return named[0], "uniquely named temperature/EIV-T column"
    raise ValueError(f"Temperature column ambiguous: signature={candidates}; named={named}; headers={columns}")

def temperature_cwm(codes, traits, col):
    known = [value(traits[x].get(col)) for x in codes]
    valid = [v for v in known if v is not None]
    return (sum(valid) / len(valid) if valid else None,
            len(valid) / len(codes) if codes else 0, len(valid), len(codes))

def source_pairs(original):
    geo, zones, species, traits, _, qa = validate_original_keys(original)
    col, rule = find_temperature_column(traits)
    pairs = []
    for (site, zone_id), zone in zones.items():
        if zone["area"] != "BS":
            continue
        if zone_id not in ("LBS", "RBS"):
            raise ValueError("Unexpected BS zone id")
        toe_id = zone_id[0] + "BT"
        other = zones.get((site, toe_id))
        if other is None:
            continue
        if other["area"] != "BT" or zone.get("area_side") != other.get("area_side"):
            raise ValueError("Mismatched physical track side")
        a = temperature_cwm(species[(site, zone_id)], traits, col)
        b = temperature_cwm(species[(site, toe_id)], traits, col)
        pairs.append({"site": site, "side": zone_id[0],
                      "segment": str(geo[site]["track_section_number"]),
                      "shoulder": a[0], "toe": b[0],
                      "shoulder_coverage": a[1], "toe_coverage": b[1],
                      "shoulder_known": a[2], "toe_known": b[2],
                      "shoulder_total": a[3], "toe_total": b[3]})
    if len(pairs) != 173 or len({x["site"] for x in pairs}) != 89:
        raise ValueError("Original matched field-site numbers no longer agree")
    return pairs, col, rule, qa

def summarize(pairs, threshold, seed=318, repetitions=2500):
    bysite = defaultdict(list)
    for p in pairs:
        if p["shoulder"] is None or p["toe"] is None:
            continue
        if min(p["shoulder_coverage"], p["toe_coverage"]) < threshold:
            continue
        bysite[(p["segment"], p["site"])].append(p["shoulder"] - p["toe"])
    grouped = defaultdict(list)
    for (segment, site), diffs in bysite.items():
        grouped[segment].append(sum(diffs) / len(diffs))
    vals = [x for v in grouped.values() for x in v]
    out = {"min_coverage_each_habitat": threshold,
           "eligible_track_side_pairs": sum(len(v) for v in bysite.values()),
           "eligible_distinct_locations": len(bysite),
           "eligible_original_rail_segments": len(grouped),
           "mean_BS_minus_BT_EIV_T": sum(vals) / len(vals) if vals else None,
           "locations_BS_warmer": sum(v > 0 for v in vals),
           "locations_BT_warmer": sum(v < 0 for v in vals),
           "locations_equal": sum(v == 0 for v in vals)}
    if len(grouped) >= 3:
        rng = random.Random(seed)
        sections = sorted(grouped)
        draws = []
        for _ in range(repetitions):
            chosen = [rng.choice(sections) for _ in sections]
            samples = [x for segment in chosen for x in grouped[segment]]
            draws.append(sum(samples) / len(samples))
        draws.sort()
        out["rail_section_block_bootstrap_95pct_exploratory"] = [
            draws[int(repetitions * .025)], draws[int(repetitions * .975)]]
    else:
        out["rail_section_block_bootstrap_95pct_exploratory"] = None
    return out

def analyze(original):
    pairs, col, rule, qa = source_pairs(original)
    return {
      "version": "v318", "scope": "WITHIN ACTIVE RAILWAY BS-versus-BT ONLY",
      "cannot_estimate_cold_region_rail_nonrail_interaction": True,
      "source_doi": "10.1007/s44473-026-00155-z",
      "source_sha256": SOURCE_SHA256,
      "source_trait_column": col, "source_trait_identification_rule": rule,
      "source_original_qa": qa,
      "original_same_side_pairs": len(pairs),
      "original_distinct_field_locations": len({p["site"] for p in pairs}),
      "source_rail_sections": len({p["segment"] for p in pairs}),
      "zone_coverage_average": {
          "BS": sum(p["shoulder_coverage"] for p in pairs) / len(pairs),
          "BT": sum(p["toe_coverage"] for p in pairs) / len(pairs)},
      "coverage_sensitivity": [summarize(pairs, t) for t in (0, .5, .7, .9)],
      "limits": ["Ordinal ecological temperature preferences, not degrees Celsius",
          "Both habitats are railway microhabitats; BT is NOT nonrail",
          "Original field zones do not have verified equal areas",
          "Both sides and repeated locations are not independent rail line replicates",
          "Bootstrap groups source rail sections, not independent countries",
          "No measured background climate or exposure at each paired railway site"]}

def main():
    a = argparse.ArgumentParser()
    a.add_argument("--source-xlsx", required=True)
    a.add_argument("--verified-v102", required=True)
    a.add_argument("--output", required=True)
    p = a.parse_args()
    original = extract_original_xlsx(p.source_xlsx, p.verified_v102)
    result = analyze(original)
    dest = Path(p.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k != "source_original_qa"}, ensure_ascii=False), flush=True)

if __name__ == "__main__":
    main()
