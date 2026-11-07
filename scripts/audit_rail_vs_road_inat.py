#!/usr/bin/env python3
"""First REAL rail-vs-road observation-contrast feasibility pilot for Tokyo.

Sources: official MLIT N02 (2025 passenger rail), N13 (2024 road
centrelines) and public iNaturalist taxon observations, 2023–2025.
Only an observational source/observer/process audit. It is NOT evidence of
mowing, train slipstream, heritable adaptation or train-carried seeds.
Raw record coordinates and IDs are never written to output artifacts.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from audit_inat_rail_phenophases import (
    AREAS, N02, URL, RateLimitedClient, acceptable_observation,
    rail_index_for_bbox
)
from audit_n13_road import N13_URL_PATTERN
import shapefile
from pyproj import Transformer
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.transport_overlap import (road_rail_status, weighted_case_control_table,
                                    monthly_mh_ratio)

METRO = "tokyo_bbox"
BBOX = AREAS[METRO]  # predeclared in v0.4, not shifted after seeing records
DATES = ("2023-01-01", "2025-12-31")
FOCUS = (
    "Trifolium repens", "Solidago altissima", "Plantago lanceolata",
    "Erigeron annuus"
)
PLANTS_INAT_TAXON_ID = 47126
N13_MESH = "5339"
MAX_API_WINDOW = 10000


def official_road_distance(raw, bbox):
    """Nearest 2024 MLIT road centreline; not route-based/traffic exposure.

    Reads only SHP + SHX to avoid loading giant DBF. Handles spatially
    intersecting line parts; geographic clipping includes a 0.03 degree
    halo so nearby lines outside bbox remain candidates.
    """
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("N13 downloaded file is not zip")
    swlat, swlng, nelat, nelng = bbox
    zone = int(((swlng + nelng) / 2 + 180) // 6) + 1
    transformer = Transformer.from_crs(4326, 32600 + zone, always_xy=True)
    lines = []
    n_shapes = 0
    in_bbox = (swlng-.03, swlat-.03, nelng+.03, nelat+.03)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names = z.namelist()
        shp = next((p for p in names if p.lower().endswith(".shp")), None)
        if shp is None:
            raise ValueError("No SHP in N13")
        stem = shp[:-4]
        member = {p.lower(): p for p in names}
        reader = shapefile.Reader(
            shp=io.BytesIO(z.read(member[(stem+".shp").lower()])),
            shx=io.BytesIO(z.read(member[(stem+".shx").lower()])),
        )
        for shape in reader.iterShapes(bbox=in_bbox):
            n_shapes += 1
            parts = list(shape.parts) + [len(shape.points)]
            for a, b in zip(parts, parts[1:]):
                if b-a < 2:
                    continue
                coords = shape.points[a:b]
                xx, yy = transformer.transform(
                    [p[0] for p in coords], [p[1] for p in coords]
                )
                line = LineString(zip(xx, yy))
                if line.length > 0:
                    lines.append(line)
        reader.close()
    if not lines:
        raise ValueError("No official road segments near Tokyo bounding box")
    tree = STRtree(lines)
    def nearest(lng, lat):
        q = Point(*transformer.transform(lng, lat))
        idx = int(tree.nearest(q))
        return float(q.distance(lines[idx]))
    return nearest, transformer, {"n_intersecting_shapefile_shapes": n_shapes,
                                   "n_indexed_line_parts": len(lines)}


def fixed_params():
    swlat, swlng, nelat, nelng = BBOX
    return {
        "quality_grade": "research", "photos": "true",
        "d1": DATES[0], "d2": DATES[1],
        "swlat": swlat, "swlng": swlng, "nelat": nelat, "nelng": nelng,
        "per_page": 200, "order_by": "id", "order": "asc"
    }


def pages_full(client, params, *, limit=MAX_API_WINDOW):
    params = {**params}
    first = client.get(URL, params={**params, "page": 1}).json()
    total = int(first.get("total_results") or 0)
    if total>limit:
        raise ValueError(f"{total} observations exceeds {limit} API window; narrow dates")
    rows = list(first.get("results") or [])
    for page in range(2, math.ceil(total/200)+1):
        rows.extend(client.get(URL, params={**params, "page": page}).json()["results"])
    if len(rows)!=total or len({r["id"] for r in rows})!=len(rows):
        raise ValueError(f"API paging inconsistency: {len(rows)} / {total}")
    return rows, total


def monthly_control_sample(client):
    """One hash-determined ID page per calendar month, NOT a probability sample."""
    records = []
    metadata = []
    base = fixed_params()
    for year in (2023, 2024, 2025):
        for month in range(1, 13):
            from calendar import monthrange
            lo = f"{year}-{month:02d}-01"
            hi = f"{year}-{month:02d}-{monthrange(year, month)[1]:02d}"
            args = {**base, "taxon_id": PLANTS_INAT_TAXON_ID, "d1": lo, "d2": hi}
            first = client.get(URL, params={**args, "page": 1}).json()
            total = int(first["total_results"])
            if total > MAX_API_WINDOW:
                raise ValueError(f"Reference month {lo} has {total}>10k; split by date")
            n_pages = max(1, math.ceil(total / 200))
            offset = int(hashlib.sha256(f"tetsu-v05-{lo}".encode()).hexdigest()[:8], 16)
            selected_page = 1 + offset % n_pages
            if selected_page == 1:
                chosen = first["results"]
            else:
                chosen = client.get(URL, params={**args, "page": selected_page}).json()["results"]
            if len({r["id"] for r in chosen}) != len(chosen):
                raise ValueError("duplicate control IDs within page")
            n = len(chosen)
            # For a randomly sampled page, Horvitz–Thompson weights would
            # require each page to have known inclusion prob. Here IDs are
            # sorted and page is selected via fixed hash; use approximation.
            approx_w = total / n if n else 1.0
            metadata.append({
                "year": year, "month": month,
                "total_inat_photo_plant_records": total,
                "sampled_page": selected_page, "available_pages": n_pages,
                "page_records_returned": n,
                "approx_inverse_page_sample_weight": round(approx_w, 4)
            })
            for ob in chosen:
                taxon = (ob.get("taxon") or {}).get("name") or ""
                if any(taxon==f or taxon.startswith(f+" ") for f in FOCUS):
                    continue
                records.append((ob, "reference", taxon, approx_w))
    return records, metadata


def aggregate():
    client = RateLimitedClient(delay_seconds=1.1)
    raw_rail = client.get(N02, timeout=(15, 120)).content
    raw_roads = client.get(N13_URL_PATTERN.format(mesh=N13_MESH), timeout=(25, 180)).content
    near_rail, n_rail = rail_index_for_bbox(raw_rail, BBOX)
    near_road, transformer, roadmeta = official_road_distance(raw_roads, BBOX)
    print(f"Loaded official Tokyo geo: rail_segments={n_rail}, road_indexed={roadmeta['n_indexed_line_parts']}",flush=True)

    observations = []
    counts = {}
    for taxon in FOCUS:
        data, n = pages_full(client, {**fixed_params(), "taxon_name":taxon})
        counts[taxon] = n
        observations.extend((x, "target", taxon, 1.0) for x in data)
        print(f"Full focal taxon {taxon} = {n}",flush=True)
    references, ref_sampling = monthly_control_sample(client)
    observations += references
    print(f"Reference pages yielded={len(references)} nonfocus plant observations",flush=True)

    quality = defaultdict(Counter)
    filtered = []
    observers = defaultdict(set)
    per_month = Counter()
    for ob, label, taxon, ref_weight in observations:
        okay, why = acceptable_observation(ob)
        quality[label][why] += 1
        if not okay:
            continue
        lon, lat = ob["geojson"]["coordinates"]
        date = (ob.get("observed_on") or "")[:10]
        try:
            year,month,_ = (int(x) for x in date.split("-"))
        except (ValueError, TypeError):
            quality[label]["date_parse_failed"] += 1
            continue
        if year not in (2023,2024,2025):
            quality[label]["outside_period"] += 1
            continue
        dist_rail = near_rail(lon,lat)
        dist_road = near_road(lon,lat)
        x,y = transformer.transform(lon,lat)
        grid2km = f"{int(math.floor(x/2000))}:{int(math.floor(y/2000))}"
        user = (ob.get("user") or {}).get("id")
        if user:
            observers[label].add(user)
        per_month[(label,month)] += 1
        filtered.append({
            "label":label,"focus_taxon":taxon if label=="target" else "other_plants",
            "rail_m":dist_rail,"road_m":dist_road,
            "grid2km":grid2km,"month":month,"year":year,
            "ref_weight":ref_weight
        })

    target = [r for r in filtered if r["label"]=="target"]
    reference = [r for r in filtered if r["label"]=="reference"]
    tables = {}
    for name, cases in [
        ("four_focal_species_pooled",target),
        *[(x,[r for r in target if r["focus_taxon"]==x]) for x in FOCUS]
    ]:
        paired = cases + reference
        primary = weighted_case_control_table(paired, rail_near=250., road_limit=100.)
        alternative = weighted_case_control_table(paired, rail_near=100., road_limit=100.)
        road150 = weighted_case_control_table(paired, rail_near=250., road_limit=150.)
        tables[name] = {
            "primary_rail_250m_road_100m":primary,
            "sensitivity_rail_100m_road_100m":alternative,
            "sensitivity_rail_250m_road_150m":road150,
            "month_stratified_primary":monthly_mh_ratio(paired),
        }
    # Check whether matched cells actually span different treatments.
    blockdata=defaultdict(lambda:{"target":set(),"reference":set()})
    for r in filtered:
        state=road_rail_status(r["rail_m"],r["road_m"])
        if state.endswith("_near"):
            blockdata[r["grid2km"]][r["label"]].add("near")
        elif state.endswith("_far"):
            blockdata[r["grid2km"]][r["label"]].add("far")
    both_rail_groups = sum(
        len(v["target"])>0 and len(v["reference"])>0 for v in blockdata.values()
    )
    return {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_type":"single-metro iNaturalist observation-distribution proxy; not causal, not railway adaptation",
        "metro_bbox": METRO,"bbox":BBOX, "observation_years": DATES,
        "rail": {"official_N02_year":2025,"sha256":hashlib.sha256(raw_rail).hexdigest(),
                 "rail_line_parts_indexed":n_rail},
        "road": {"official_N13_year":2024,"mesh":N13_MESH,
                 "sha256":hashlib.sha256(raw_roads).hexdigest(),**roadmeta},
        "full_focal_taxon_records":counts,
        "plant_reference_api_month_page_sampling":ref_sampling,
        "n_eligible_target":len(target), "n_eligible_reference":len(reference),
        "geo_quality_rejection":{k:dict(v) for k,v in quality.items()},
        "distinct_observers":{k:len(v) for k,v in observers.items()},
        "n_2km_blocks_with_both_labels_under_road_matching":both_rail_groups,
        "tables":tables,
        "limitations":[
            "Observations are observer-selected photographed plants; neither background nor non-observations is a plant absence.",
            "Plant reference is one hash-chosen ID-sorted iNaturalist page/month, not a probability sample. Month-specific size weights are exploratory only.",
            "Official 2024 road and 2025 passenger rail centrelines are geometry, NOT freight usage, train wind speeds or vegetation cutting exposure.",
            "Rail <=250m may be station gardens, parks, roadsides, or other non-trackside habitat; coordinate uncertainty up to 100m.",
            "Only a single Tokyo metro bbox and four focal herbs were inspected, no independent external-region holdout.",
            "A ratio above 1 would be a photographic-availability association, not an ecological seed dispersal, mowing, adaptation or phenological shift.",
            "No uncertainty intervals or p-values are issued until independent spatial/observer blocking and inclusion-bias models pass."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output", default="results/railroad/tokyo_rail_vs_road_photo_availability.json")
    args=p.parse_args()
    report=aggregate()
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "total_eligible_target":report["n_eligible_target"],
        "total_eligible_reference":report["n_eligible_reference"],
        "source_road_lines":report["road"]["n_indexed_line_parts"],
        "n_observer":report["distinct_observers"],
        "n_2km_blocks_with_rail_road_focal_reference":report["n_2km_blocks_with_both_labels_under_road_matching"],
        "comparisons":report["tables"],
    },ensure_ascii=False,indent=2), flush=True)


if __name__=="__main__":
    main()
