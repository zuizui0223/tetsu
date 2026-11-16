"""Pairwise similarity from preprocessed, independently sourced input tables.

No code here demonstrates physical seed movement by trains.
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


def read_incidence(path: str | Path) -> dict[str, dict[str, set[str]]]:
    """Return region -> alien/native taxon sets. Fail on bad labels or duplicates."""
    data = defaultdict(lambda: {"alien": set(), "native": set()})
    seen = set()
    taxon_status = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not {"region_id", "taxon", "status"}.issubset(reader.fieldnames or []):
            raise ValueError("incidence requires region_id,taxon,status")
        for row in reader:
            region, taxon, status = (row[k].strip() for k in ("region_id", "taxon", "status"))
            if not region or not taxon or status not in ("alien", "native"):
                raise ValueError("blank region/taxon or invalid taxon status")
            if (region, taxon) in seen:
                raise ValueError(f"duplicate region/taxon: {region}/{taxon}")
            if taxon in taxon_status and taxon_status[taxon] != status:
                raise ValueError(f"contradictory taxon status: {taxon}")
            seen.add((region, taxon))
            taxon_status[taxon] = status
            data[region][status].add(taxon)
    return dict(data)


def jaccard(left: set[str], right: set[str]) -> float | None:
    """Return Jaccard similarity; None if either endpoint lacks observation."""
    if not left or not right:
        return None
    return len(left & right) / len(left | right)


def nonnegative(row: dict[str, str], key: str, positive: bool = False) -> float:
    val = float(row[key])
    if not math.isfinite(val) or val < 0 or (positive and val == 0):
        raise ValueError(f"{key} must be finite and {'positive' if positive else 'nonnegative'}")
    return val


def prepare_pairwise(incidence: dict, pairs_csv: str | Path) -> list[dict]:
    outputs = []
    seen_pairs = set()
    columns = {"region_a", "region_b", "euclid_km", "road_km", "rail_km",
               "climate_distance", "urban_difference", "effort_a", "effort_b"}
    with open(pairs_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = columns - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"pairs missing columns: {sorted(missing)}")
        for row in reader:
            a, b = row["region_a"].strip(), row["region_b"].strip()
            if a == b or a not in incidence or b not in incidence:
                raise ValueError(f"invalid or unknown region pair: {a}/{b}")
            pair = tuple(sorted((a, b)))
            if pair in seen_pairs:
                raise ValueError(f"duplicate unordered pair: {a}/{b}")
            seen_pairs.add(pair)
            geo = nonnegative(row, "euclid_km", positive=True)
            road = nonnegative(row, "road_km", positive=True)
            climate = nonnegative(row, "climate_distance")
            urban = nonnegative(row, "urban_difference")
            if urban > 1:
                raise ValueError("urban_difference must be at most 1")
            ea = nonnegative(row, "effort_a", positive=True)
            eb = nonnegative(row, "effort_b", positive=True)
            rail = nonnegative(row, "rail_km", positive=True) if (row["rail_km"] or "").strip() else None
            aa, bb = incidence[a], incidence[b]
            aj = jaccard(aa["alien"], bb["alien"])
            nj = jaccard(aa["native"], bb["native"])
            outputs.append({
                "region_a": a, "region_b": b, "euclid_km": geo, "road_km": road,
                "rail_km": rail, "rail_connected": int(rail is not None),
                "rail_directness": (math.log1p(road) - math.log1p(rail)) if rail is not None else 0.,
                "climate_distance": climate, "urban_difference": urban,
                "log_effort_min": math.log1p(min(ea, eb)),
                "alien_jaccard": aj, "native_jaccard": nj,
                "alien_minus_native": (aj - nj) if aj is not None and nj is not None else None,
                "alien_richness_a": len(aa["alien"]), "alien_richness_b": len(bb["alien"]),
                "native_richness_a": len(aa["native"]), "native_richness_b": len(bb["native"])
            })
    return outputs
