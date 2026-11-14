#!/usr/bin/env python3
"""Japan real-data pilot: iNaturalist phenophase availability around rail corridors.

Retrieve only research-grade, photographed observations and explicit flower /
no-flowers-fruits annotations. Assess *observed* coordinates against the official
2025 MLIT passenger-rail geometry. Strictly a sample-coverage diagnostic:
not actual flowering onset, mowing effects, railway adaptation or causation.

No image downloads. Raw public coordinates/observation IDs are not exported.
See docs/inat_rail_phenophase_pilot.md.
"""
from __future__ import annotations

import argparse
import io
import json
import math
import time
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import requests

URL = "https://api.inaturalist.org/v1/observations"
N02 = "https://nlftp.mlit.go.jp/ksj/gml/data/N02/N02-25/N02-25_GML.zip"
STAGES = {"flowers": 13, "no_flowers_fruits": 21}
# Three bounding boxes deliberately fixed before reviewing results; not full
# administrative city boundaries. The same site definition is used each run.
AREAS = {
    "tokyo_bbox": (35.52, 139.56, 35.83, 139.92),
    "osaka_bbox": (34.53, 135.35, 34.83, 135.76),
    "sendai_bbox": (38.16, 140.72, 38.40, 141.04),
}
# First four are confirmed "accidental-only" in the Egawa & Koyama national
# checklist. Others are exploratory urban common plants; no pathway assigned.
TAXA = {
    "Briza minor": "accidental_only",
    "Lepidium virginicum": "accidental_only",
    "Rumex crispus": "accidental_only",
    "Sonchus asper": "accidental_only",
    "Trifolium repens": "exploratory_not_pathway_classified",
    "Plantago lanceolata": "exploratory_not_pathway_classified",
    "Solidago altissima": "exploratory_not_pathway_classified",
    "Erigeron annuus": "exploratory_not_pathway_classified",
}
ALLOWED_PHOTO_LICENSES = {"cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"}


def near_rail_bin(distance_m: float) -> str:
    if distance_m <= 250:
        return "0-250m"
    if distance_m <= 1000:
        return "250-1000m"
    return ">1000m"


def acceptable_observation(obs: dict, *, max_accuracy_m: float = 100.0) -> tuple[bool, str]:
    if obs.get("obscured") or obs.get("geoprivacy") in ("obscured", "private"):
        return False, "obscured_or_private"
    if obs.get("captive"):
        return False, "captive_or_cultivated"
    accuracy = obs.get("positional_accuracy")
    if isinstance(accuracy, bool) or not isinstance(accuracy, (float, int)) or not 0 < accuracy <= max_accuracy_m:
        return False, "accuracy_not_within_100m"
    coords = (obs.get("geojson") or {}).get("coordinates")
    if not isinstance(coords, list) or len(coords) != 2:
        return False, "no_valid_coordinates"
    lng, lat = coords
    if not isinstance(lng, (float, int)) or not isinstance(lat, (float, int)) or not (
        math.isfinite(lng) and math.isfinite(lat) and -180 <= lng <= 180 and -90 <= lat <= 90
    ):
        return False, "invalid_coordinates"
    if not obs.get("observed_on"):
        return False, "missing_observation_date"
    if not obs.get("photos"):
        return False, "no_photos"
    return True, "pass"


class RateLimitedClient:
    def __init__(self, delay_seconds=1.10):
        self.s = requests.Session()
        self.s.headers.update({"User-Agent": "tetsu-ecology-feasibility/0.4 (https://github.com/zuizui0223/tetsu)"})
        self.last = 0.0
        self.delay = delay_seconds

    def get(self, url, params=None, timeout=(15, 90)):
        for attempt in range(5):
            delay = self.delay - (time.monotonic() - self.last)
            if delay > 0:
                time.sleep(delay)
            self.last = time.monotonic()
            r = self.s.get(url, params=params, timeout=timeout)
            if r.status_code in (429, 502, 503, 504):
                time.sleep(min(20.0, 2.0 * (attempt + 1)))
                continue
            r.raise_for_status()
            return r
        raise RuntimeError(f"iNaturalist request did not succeed after retries: {url}")


def rail_index_for_bbox(raw_zip: bytes, bbox: tuple[float, float, float, float]):
    """Build local UTM spatial tree from official rail section geometry.

    The index is not a rail *service network*, and nearest-track proximity
    does not demonstrate arrival by train or train wind exposure.
    """
    import shapefile
    from pyproj import Transformer
    from shapely.geometry import LineString, Point, box
    from shapely.strtree import STRtree

    swlat, swlng, nelat, nelng = bbox
    zone = int(((swlng + nelng) / 2 + 180) // 6) + 1
    transformer = Transformer.from_crs(4326, 32600 + zone, always_xy=True)
    with zipfile.ZipFile(io.BytesIO(raw_zip)) as z:
        members = {n.lower(): n for n in z.namelist()}
        path = next(n for n in z.namelist() if n.endswith("UTF-8/N02-25_RailroadSection.shp"))
        stem = path[:-4]
        def stream(ext):
            return io.BytesIO(z.read(members[(stem + ext).lower()]))
        reader = shapefile.Reader(shp=stream(".shp"), shx=stream(".shx"), dbf=stream(".dbf"), encoding="utf-8")
        lines = []
        for shape in reader.iterShapes():
            # First cheap lat/lon bbox gate, with generous 0.04deg halo.
            a, b, c, d = shape.bbox
            if c < swlng - 0.05 or a > nelng + 0.05 or d < swlat - 0.05 or b > nelat + 0.05:
                continue
            parts = list(shape.parts) + [len(shape.points)]
            for p, q in zip(parts, parts[1:]):
                if q - p < 2:
                    continue
                pts = [transformer.transform(x, y) for x, y in shape.points[p:q]]
                ls = LineString(pts)
                if ls.length > 0:
                    lines.append(ls)
        reader.close()
    if not lines:
        raise RuntimeError(f"No MLIT rail geometry within {bbox}")
    tree = STRtree(lines)
    def distance_m(lng, lat):
        pt = Point(*transformer.transform(lng, lat))
        nearest = tree.nearest(pt)
        return float(pt.distance(lines[int(nearest)]))
    return distance_m, len(lines)


def run(*, species, areas, sample_limit: int, delay: float):
    client = RateLimitedClient(delay)
    raw_rail = client.get(N02, timeout=(15, 120)).content
    if not zipfile.is_zipfile(io.BytesIO(raw_rail)):
        raise ValueError("official 2025 N02 download was not a ZIP archive")
    n02_sha = __import__("hashlib").sha256(raw_rail).hexdigest()
    dist_maps = {}
    segment_counts = {}
    for place, bbox in areas.items():
        dist, n = rail_index_for_bbox(raw_rail, bbox)
        dist_maps[place] = dist
        segment_counts[place] = n
    out = []
    for place, bbox in areas.items():
        swlat, swlng, nelat, nelng = bbox
        for taxon, pathway in species.items():
            fixed = {
                "taxon_name": taxon, "quality_grade": "research", "photos": "true",
                "d1": "2018-01-01", "d2": "2025-12-31", "swlat": swlat,
                "swlng": swlng, "nelat": nelat, "nelng": nelng,
            }
            totals = client.get(URL, params={**fixed, "per_page": 1}).json()
            n_any = int(totals["total_results"])
            for stage, stage_id in STAGES.items():
                query = {
                    **fixed, "term_id": 12, "term_value_id": stage_id,
                    "per_page": sample_limit, "order_by": "id", "order": "asc",
                    "page": 1,
                }
                data = client.get(URL, params=query).json()
                n_phase = int(data["total_results"])
                counts = Counter()
                geo_bins = Counter()
                licensed_geo_bins = Counter()
                distinct_observers = set()
                for ob in data.get("results", []):
                    valid, why = acceptable_observation(ob)
                    counts[why] += 1
                    if valid:
                        lng, lat = ob["geojson"]["coordinates"]
                        d = dist_maps[place](lng, lat)
                        b = near_rail_bin(d)
                        geo_bins[b] += 1
                        u = ob.get("user") or {}
                        if u.get("id"):
                            distinct_observers.add(u["id"])
                        # Metadata only; no image fetched or mask inferred.
                        if any(str(ph.get("license_code") or "").lower() in ALLOWED_PHOTO_LICENSES for ph in ob["photos"]):
                            licensed_geo_bins[b] += 1
                out.append({
                    "bbox": place, "taxon": taxon, "pathway_group": pathway,
                    "phase": stage, "total_species_all_phenophases": n_any,
                    "annotated_phase_total": n_phase,
                    "returned_phase_records": len(data.get("results", [])),
                    "sample_truncated": bool(n_phase > sample_limit),
                    "valid_geo_sample_n": sum(geo_bins.values()),
                    "valid_geo_unique_observers": len(distinct_observers),
                    "geo_bins": dict(geo_bins), "licensed_photo_geo_bins": dict(licensed_geo_bins),
                    "rejection_reasons": dict(counts),
                })
            print(f"{place} / {taxon}: all={n_any} annotated={','.join(str(out[-j]['annotated_phase_total']) for j in (2,1))}", flush=True)
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "real source feasibility sampling; NOT a railway phenology effect",
        "inat_api": URL, "n02_sha256": n02_sha,
        "n02_segment_geometries": segment_counts, "bbox_definitions": areas,
        "taxa": species, "sample_cap_per_taxon_stage_area": sample_limit,
        "rows": out,
        "critical_caveats": [
            "iNaturalist annotations overlap and are not standardized complete-visit nonflowering denominators.",
            "No Flowers or Fruits refers to visibility in the photo, not biological nonflowering status.",
            "API per-stage result capped at sample_limit, ID-sorted and nonrandom; not a population estimate.",
            "No road-distance, microclimate, mowing or slipstream data have been matched.",
            "N02 is passenger rail; it excludes freight-only lines. A nearby line need not mean train-carried seed.",
            "Only observations with recorded <=100m coordinate accuracy, unobscured non-captive and visible photo are spatially classified.",
            "Photos are not downloaded; photo-mask ROI has not been computed.",
            "Sensitive coordinates and observation IDs are not exported; only aggregates."
        ]
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/inat/rail_phenophase_audit.json")
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--delay", type=float, default=1.10)
    args = parser.parse_args()
    if args.limit < 1 or args.limit > 200 or args.delay < 1.0:
        parser.error("limit must be 1..200 and delay >=1s to respect iNaturalist guidance")
    report = run(species=TAXA, areas=AREAS, sample_limit=args.limit, delay=args.delay)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"summary": {
        "n_strata": len(report["rows"]),
        "annotated_records_api_sum_with_possible_overlap": sum(r["annotated_phase_total"] for r in report["rows"]),
        "quality_filtered_sample": sum(r["valid_geo_sample_n"] for r in report["rows"]),
        "georeferenced_near_250m_sample": sum(r["geo_bins"].get("0-250m", 0) for r in report["rows"]),
        "truncated_strata": sum(r["sample_truncated"] for r in report["rows"])
    }}, indent=2))


if __name__ == "__main__":
    main()
