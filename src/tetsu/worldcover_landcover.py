"""ESA WorldCover 2021 v200 10m land-cover quality screens for rail candidates.

This map is processed Sentinel-1/2 imagery for 2021, NOT a field
inventory, railway vegetation map, at-grade railway source, or independent
replicate of 2025 Sentinel-2 SCL greenness. Tree cover may mask grasses
below; grassland class need not be on railway property.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import math
import statistics

SOURCE_BASE = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map"
SOURCE_ALT = "https://esa-worldcover.s3.amazonaws.com/v200/2021/map"
DOI = "10.5281/zenodo.7254221"
REF_YEAR = 2021
CLASS_LABELS = {
    10: "tree_cover",
    20: "shrubland",
    30: "grassland",
    40: "cropland",
    50: "built_up",
    60: "bare_sparse_vegetation",
    70: "snow_ice",
    80: "permanent_water",
    90: "herbaceous_wetland",
    95: "mangroves",
    100: "moss_lichen",
}
UNKNOWN = 0
POTENTIAL_GROUND_HERBS = {30, 90}
NON_TREE_POTENTIALLY_NONCROP = {20, 30, 60, 90}


def tile_name(lon: float, lat: float) -> str:
    """WorldCover 3x3-degree grid tile named by SOUTHWEST corner."""
    if not (-180 <= lon < 180 and -90 <= lat < 90):
        raise ValueError("Invalid WorldCover tile coordinate")
    southwest_lat = 3 * math.floor(lat / 3)
    southwest_lon = 3 * math.floor(lon / 3)
    ns = "N" if southwest_lat >= 0 else "S"
    ew = "E" if southwest_lon >= 0 else "W"
    return f"{ns}{abs(southwest_lat):02d}{ew}{abs(southwest_lon):03d}"


def tile_url(tile: str, base: str = SOURCE_BASE) -> str:
    if len(tile) != 7 or tile[0] not in "NS" or tile[3] not in "EW":
        raise ValueError("Unexpected WorldCover tile format")
    if not (tile[1:3].isdigit() and tile[4:].isdigit()):
        raise ValueError("Invalid WorldCover tile number")
    if base not in (SOURCE_BASE, SOURCE_ALT):
        raise ValueError("Untrusted WorldCover source")
    return f"{base}/ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"


def normalize_class(value) -> int | None:
    if value is None:
        return None
    try:
        number = int(value)
    except (TypeError,ValueError):
        raise ValueError(f"Invalid ESA WorldCover class {value!r}")
    if number == UNKNOWN:
        return None
    if number not in CLASS_LABELS:
        raise ValueError(f"Unknown ESA WorldCover class code {number}")
    return number


def summarize_neighborhood(pixel_classes) -> dict:
    """3x3 independent map pixel labels; 0 is missing, not bare ground."""
    observed = [normalize_class(v) for v in pixel_classes]
    valid = [v for v in observed if v is not None]
    c = Counter(valid)
    # A plurality requires a deterministic tie-break rather than
    # arbitrarily inferring forest or grassy railway at a mixed edge.
    top = sorted(c.items(), key=lambda kv: (-kv[1],kv[0]))
    tied = (len(top)>1 and top[0][1]==top[1][1])
    return {
        "n_window_pixels":len(observed),
        "n_classified":len(valid),
        "n_unclassified":len(observed)-len(valid),
        "tree_cover_fraction_among_classified":(
            round(c[10]/len(valid),5) if valid else None
        ),
        "grassland_fraction_among_classified":(
            round(c[30]/len(valid),5) if valid else None
        ),
        "plurality_class":(None if not top or tied else top[0][0]),
        "plurality_tied":tied,
    }


def summarize_route_source(site_rows: list[dict]) -> dict:
    """Site-level geometry-chosen rail/road WorldCover labels as aggregates."""
    near = Counter()
    road = Counter()
    both_coded = 0
    at_least_one_ground_herb = 0
    both_ground_herb = 0
    near_ground_herb = 0
    near_tree = 0
    near_plurals = Counter()
    tree_fractions=[]
    for row in site_rows:
        a = normalize_class(row.get("near_class"))
        b = normalize_class(row.get("control_class"))
        near[CLASS_LABELS[a] if a is not None else "unknown"] += 1
        road[CLASS_LABELS[b] if b is not None else "unknown"] += 1
        if a is not None and b is not None:
            both_coded += 1
        if a in POTENTIAL_GROUND_HERBS:
            near_ground_herb += 1
        if a == 10:
            near_tree += 1
        if a in POTENTIAL_GROUND_HERBS or b in POTENTIAL_GROUND_HERBS:
            at_least_one_ground_herb += 1
        if a in POTENTIAL_GROUND_HERBS and b in POTENTIAL_GROUND_HERBS:
            both_ground_herb += 1
        patch=row.get("near_window") or {}
        plurality=normalize_class(patch.get("plurality_class"))
        if plurality is not None:
            near_plurals[CLASS_LABELS[plurality]] += 1
        else:
            near_plurals["unknown_or_tied"] += 1
        frac=patch.get("tree_cover_fraction_among_classified")
        if isinstance(frac,(float,int)) and math.isfinite(frac):
            tree_fractions.append(frac)
    return {
        "n_rail_road_geometry_selected_pairs":len(site_rows),
        "near_2021_worldcover_class_counts":dict(near),
        "road_control_2021_worldcover_class_counts":dict(road),
        "n_both_2021_class_known":both_coded,
        "n_near_2021_worldcover_tree_cover":near_tree,
        "n_near_2021_worldcover_grassland_or_herbaceous_wetland":near_ground_herb,
        "n_near_or_control_potential_ground_herbaceous":at_least_one_ground_herb,
        "n_near_and_control_potential_ground_herbaceous":both_ground_herb,
        "near_3x3_plurality_class_counts":dict(near_plurals),
        "median_3x3_near_tree_fraction_if_known":(
            round(statistics.median(tree_fractions),5)
            if tree_fractions else None
        ),
        "2021_landcover_not_2025_or_actual_track_biota":True,
    }
