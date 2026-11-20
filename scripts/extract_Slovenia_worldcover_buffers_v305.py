"""Acquire 2021 ESA WorldCover landscape covariates; no botanical outcomes fitted.

Local-only cache keyed by original source site ID. Official class 50 includes
railways, so return both circles and 100m-excluded annuli. Pixel centers enter
by WGS84 geodesic distance, with ellipsoidal pixel-area weights. Tiles are
read separately across boundaries; no interpolation of categorical classes.
"""
from __future__ import annotations
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
import rasterio
from rasterio.windows import Window
from pyproj import Geod

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from tetsu.worldcover_landcover import tile_name, tile_url, CLASS_LABELS, DOI
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites

GEOD = Geod(ellps='WGS84')
RADII = (500, 1000)


def summarize_arrays(classes, distances, areas, radius, inner=0):
    mask = (distances <= radius) & (distances >= inner)
    valid = mask & (classes != 0)
    if np.any(valid & ~np.isin(classes, list(CLASS_LABELS))):
        raise ValueError('Unexpected WorldCover class')
    total, known = float(areas[mask].sum()), float(areas[valid].sum())
    if total <= 0 or known <= 0:
        raise ValueError('No mapped buffer support')
    result = {'radius_m': radius, 'inner_exclusion_m': inner,
              'n_pixels': int(mask.sum()), 'n_valid_pixels': int(valid.sum()),
              'sampled_area_m2': total, 'valid_area_fraction': known / total}
    for code, name in CLASS_LABELS.items():
        result[name + '_fraction'] = float(areas[valid & (classes == code)].sum()) / known
    return result


def read_site(site, datasets, stack):
    lon, lat = site['lon'], site['lat']
    west = GEOD.fwd(lon, lat, 270, 1020)[0]
    east = GEOD.fwd(lon, lat, 90, 1020)[0]
    south = GEOD.fwd(lon, lat, 180, 1020)[1]
    north = GEOD.fwd(lon, lat, 0, 1020)[1]
    tiles = sorted({tile_name(x, y) for x in (west, east) for y in (south, north)})
    arrays, distances, areas, proofs = [], [], [], []
    for tile in tiles:
        if tile not in datasets:
            datasets[tile] = stack.enter_context(rasterio.open(tile_url(tile)))
        src = datasets[tile]
        if src.crs.to_epsg() != 4326 or src.count != 1 or src.nodata not in (0, None):
            raise ValueError('Unexpected official raster schema')
        top, left = src.index(west, north)
        bottom, right = src.index(east, south)
        top, left = max(0, top), max(0, left)
        bottom, right = min(src.height, bottom + 1), min(src.width, right + 1)
        if bottom <= top or right <= left:
            continue
        window = Window(left, top, right-left, bottom-top)
        arr = src.read(1, window=window)
        tr = src.window_transform(window)
        xs = tr.c + (np.arange(arr.shape[1]) + .5) * tr.a
        ys = tr.f + (np.arange(arr.shape[0]) + .5) * tr.e
        xx, yy = np.meshgrid(xs, ys)
        _, _, dist = GEOD.inv(np.full(xx.shape, lon), np.full(yy.shape, lat), xx, yy)
        # Pixel area varies by latitude, constant across longitude within a row.
        row_areas = []
        for y in ys:
            a, _ = GEOD.polygon_area_perimeter(
                [xs[0]-tr.a/2, xs[0]+tr.a/2, xs[0]+tr.a/2, xs[0]-tr.a/2],
                [y-tr.e/2, y-tr.e/2, y+tr.e/2, y+tr.e/2])
            row_areas.append(abs(a))
        weights = np.broadcast_to(np.asarray(row_areas)[:, None], arr.shape)
        arrays.append(arr.ravel()); distances.append(dist.ravel()); areas.append(weights.ravel())
        proofs.append({'tile': tile, 'source_url': tile_url(tile),
                       'window_shape': list(arr.shape),
                       'pixel_bytes_sha256': hashlib.sha256(arr.tobytes()).hexdigest()})
    classes = np.concatenate(arrays)
    distance = np.concatenate(distances)
    weight = np.concatenate(areas)
    summaries = {}
    for radius in RADII:
        for inner in (0, 100):
            summaries[f'r{radius}_inner{inner}'] = summarize_arrays(classes, distance, weight, radius, inner)
    return {'source_id': site['source_id'], 'section': site['section'],
            'buffers': summaries, 'raster_windows': proofs}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='../si/source_original/SLO_2026_rail_vegetation_microhabitats_original.xlsx')
    p.add_argument('--audit', default='../si/results/audit/Slovenia_2026_RAIL_92sites_plant_supplement_source_v102.json')
    p.add_argument('--output', default='../v305/urbanization.json')
    args = p.parse_args()
    original = extract_original_xlsx(args.source, args.audit)
    grouped, qa = to_original_sites(original)
    sites = [s for group in grouped.values() for s in group]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {'source': 'ESA WorldCover 2021 v200', 'doi': DOI,
              'extraction_utc': datetime.now(timezone.utc).isoformat(),
              'source_plant_sha256': hashlib.sha256(Path(args.source).read_bytes()).hexdigest(),
              'source_cohort': qa['source_cohort'],
              'method': 'WGS84 geodesic pixel-center circles; ellipsoidal pixel-area weighting; no resampling',
              'warning': 'Built-up includes railroads; this is 2021 landscape context, not measured 2022/23 urbanization, soil, wind or management.',
              'sites': []}
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',
                      CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif', GDAL_HTTP_TIMEOUT='60'), ExitStack() as stack:
        datasets = {}
        for i, site in enumerate(sites):
            result['sites'].append(read_site(site, datasets, stack))
            output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
            print(f'WorldCover buffer {i+1}/{len(sites)}', flush=True)
    result['complete'] = len(result['sites']) == 89
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
