#!/usr/bin/env python3
"""railstatus: world rail km per 100km/200km equal-area cell, split by
operating status and rail class, from the same 16 SHA-verified Overture
2026-11-18 rail-only shards used by v370.

Categories (each Overture rail segment goes to at most one):
  active_normal      standard/narrow/unknown gauge, no abandoned/disused/
                     under-construction flag. MUST reproduce v370 exactly.
  inactive_normal    standard/narrow/unknown gauge flagged is_abandoned or
                     is_disused, not is_under_construction.
  subway             class subway, no inactive flag.
  tram_lightrail     class tram or light_rail, no inactive flag.

Purpose: placebo / mechanism contrasts for the v371-v382 finding that rail
km predicts landscape-level iNaturalist plant richness. Output is coarse
grid km only; no GERS ids or coordinates are written.
"""
from __future__ import annotations
import argparse, json, time
from collections import Counter, defaultdict
from pathlib import Path

from build_global_Overture_physical_connector_components_v361 import (
    validate_files, SOURCE_RELEASE, NORMAL_CLASSES, EXCLUDED_INACTIVE)
from sample_entire_global_GERS_railkm_WGS84_v368b import (
    stream_source_rails, seg_m, EXPECTED_CANDIDATE_SEGMENTS, EXPECTED_CANDIDATE_KM)
from aggregate_entire_world_GERS_rail_100km_200km_grid_v370 import (
    clip_source_line_to_tiles)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import CELL_SIZE_M

CATEGORIES = ("active_normal", "inactive_normal", "subway", "tram_lightrail")
INACTIVE_STATUS = frozenset(("is_abandoned", "is_disused"))


def flags_of(row):
    out = set()
    for f in row.get("rail_flags") or ():
        if isinstance(f, dict):
            out.update(f.get("values") or ())
    return out


def category(row):
    cls = row.get("class")
    flags = flags_of(row)
    if cls in NORMAL_CLASSES:
        if not (flags & EXCLUDED_INACTIVE):
            return "active_normal"
        if (flags & INACTIVE_STATUS) and "is_under_construction" not in flags:
            return "inactive_normal"
        return None
    if flags & EXCLUDED_INACTIVE:
        return None
    if cls == "subway":
        return "subway"
    if cls in ("tram", "light_rail"):
        return "tram_lightrail"
    return None


def rollup(cells100, width):
    factor = width // CELL_SIZE_M
    g = defaultdict(float)
    for (ix, iy), m in cells100.items():
        g[ix // factor, iy // factor] += m
    return [{"source_equal_area_grid_ix": ix, "source_equal_area_grid_iy": iy,
             "km": round(m / 1000.0, 5)} for (ix, iy), m in sorted(g.items()) if m > 0]


def run(rail_dir, census_path):
    from pyproj import Transformer, Geod
    census = json.loads(Path(census_path).read_text(encoding="utf8"))
    paths = validate_files(rail_dir, census)
    fwd = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)
    inv = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)
    geod = Geod(ellps="WGS84")
    cells = {c: defaultdict(float) for c in CATEGORIES}
    totals = Counter(); segs = Counter(); unassigned = Counter(); qa = Counter()
    flag_by_class = Counter()
    t0 = time.monotonic()
    for row in stream_source_rails(paths):
        qa["all_rail_rows"] += 1
        fl = flags_of(row)
        for f in fl:
            flag_by_class[f"{row.get('class')}|{f}"] += 1
        cat = category(row)
        if cat is None:
            qa["not_in_any_category"] += 1
            continue
        geom, L = seg_m(row, geod)
        totals[cat] += L; segs[cat] += 1
        try:
            assign, omitted, status = clip_source_line_to_tiles(geom, L, fwd, inv, geod)
        except ValueError:
            if cat == "active_normal":
                raise
            unassigned[cat] += L; qa[f"{cat}_clip_error"] += 1
            continue
        unassigned[cat] += omitted
        qa[f"{cat}_{status}"] += 1
        for k, m in assign.items():
            cells[cat][k] += m
        if qa["all_rail_rows"] % 250000 == 0:
            print(json.dumps({"rows": qa["all_rail_rows"],
                              "km": {c: round(totals[c] / 1000, 1) for c in CATEGORIES},
                              "s": round(time.monotonic() - t0, 1)}), flush=True)
    if qa["all_rail_rows"] != 2054855:
        raise ValueError("Overture rail row count differs from v357 census")
    if segs["active_normal"] != EXPECTED_CANDIDATE_SEGMENTS:
        raise ValueError("active_normal does not reproduce v370 candidate segment count")
    if abs(totals["active_normal"] / 1000 - EXPECTED_CANDIDATE_KM) > .04:
        raise ValueError("active_normal km does not reproduce v370")
    grids = {c: {str(w): rollup(cells[c], w) for w in (100000, 200000)} for c in CATEGORIES}
    return {
        "version": "railstatus-v1-world-rail-km-by-status-and-class",
        "source_Overture_release": SOURCE_RELEASE,
        "definitions": {
            "active_normal": "standard/narrow/unknown gauge without is_abandoned/is_disused/is_under_construction (=v370)",
            "inactive_normal": "standard/narrow/unknown gauge flagged is_abandoned or is_disused, not under construction",
            "subway": "class subway without inactive flags",
            "tram_lightrail": "class tram or light_rail without inactive flags"},
        "total_km": {c: round(totals[c] / 1000, 3) for c in CATEGORIES},
        "segments": dict(segs),
        "unassigned_km": {c: round(unassigned[c] / 1000, 3) for c in CATEGORIES},
        "qa": dict(qa),
        "flag_by_class_counts": dict(flag_by_class),
        "grids_km_by_category_and_width": grids,
        "note": "Map status flags come from OSM-derived Overture tags; abandoned/disused coverage is incomplete and uneven by country."}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--rail-dir", required=True)
    p.add_argument("--census", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    obj = run(a.rail_dir, a.census)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf8")
    print(json.dumps({k: v for k, v in obj.items() if k != "grids_km_by_category_and_width"},
                     ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
