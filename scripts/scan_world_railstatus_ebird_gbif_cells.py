#!/usr/bin/env python3
"""Independent-network replication source: eBird (GBIF-mediated eBird
Observation Dataset, EOD) bird species, record counts and observer counts per
200 km EPSG:6933 cell, 2023-2025 pooled, for every cell that holds mapped
active rail (both ends of a candidate pair must have active rail).

Uses GBIF occurrence-search facets only (no record download): speciesKey and
recordedBy facets inside the exact lon/lat rectangle of each equal-area cell
(EPSG:6933 is cylindrical, so every cell is a lon/lat rectangle). Adjacent
pair Jaccard/Simpson components are computed inside CI; species lists and
observer ids are not written.
"""
from __future__ import annotations
import argparse, json, math, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

EOD = "4fa7b334-ce0d-4e88-aaae-2e0c138d049e"
AVES = 212
YEARS = "2023,2025"
W = 200_000
X0, Y0 = -17_400_000, -7_400_000
SP_LIMIT = 4000
OBS_LIMIT = 200000
API = "https://api.gbif.org/v1/"
GATE = (50, 5, 20)  # observations, observers, species: same as v369 pair gate


def get(path, params, tries=6):
    url = API + path + ("?" + urllib.parse.urlencode(params, doseq=True) if params else "")
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "tetsu-railstatus/1"}), timeout=300) as r:
                return json.loads(r.read())
        except Exception as e:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(5 * 2 ** i)


def cell_bounds(ix, iy, inv):
    x0 = X0 + ix * W; x1 = x0 + W; y0 = Y0 + iy * W; y1 = y0 + W
    lon0, lat0 = inv.transform(x0, y0); lon1, lat1 = inv.transform(x1, y1)
    lat0 = -90.0 if not math.isfinite(lat0) else max(-90.0, lat0)
    lat1 = 90.0 if not math.isfinite(lat1) else min(90.0, lat1)
    return lon0, lon1, lat0, lat1


def query(key, inv):
    lon0, lon1, lat0, lat1 = cell_bounds(*key, inv)
    p = {"datasetKey": EOD, "taxonKey": AVES, "year": YEARS, "occurrenceStatus": "PRESENT",
         "decimalLatitude": f"{lat0:.6f},{lat1:.6f}", "decimalLongitude": f"{lon0:.6f},{lon1:.6f}",
         "limit": 0, "facet": ["speciesKey", "recordedBy"], "facetMincount": 1,
         "speciesKey.facetLimit": SP_LIMIT, "recordedBy.facetLimit": OBS_LIMIT}
    js = get("occurrence/search", p)
    f = {x["field"]: x["counts"] for x in js.get("facets", [])}
    sp = [int(c["name"]) for c in f.get("SPECIES_KEY", [])]
    ob = f.get("RECORDED_BY", [])
    return key, {"records": int(js["count"]), "species": sp, "observers": len(ob),
                 "species_truncated": len(sp) >= SP_LIMIT, "observers_truncated": len(ob) >= OBS_LIMIT,
                 "bounds": [round(v, 6) for v in (lon0, lon1, lat0, lat1)]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    from pyproj import Transformer
    inv = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)
    ds = get("dataset/" + EOD, {})
    if "ebird" not in ds.get("title", "").lower():
        raise ValueError("GBIF dataset is not eBird: " + ds.get("title", ""))
    per_year = {y: get("occurrence/search", {"datasetKey": EOD, "taxonKey": AVES, "year": y, "limit": 0})["count"]
                for y in ("2023", "2024", "2025")}
    print(json.dumps({"dataset": ds["title"], "records_by_year": per_year}), flush=True)
    if min(per_year.values()) < 1_000_000:
        raise ValueError("EOD year coverage incomplete")
    status = json.loads(Path(a.status).read_text())
    keys = sorted({(z["source_equal_area_grid_ix"], z["source_equal_area_grid_iy"])
                   for z in status["grids_km_by_category_and_width"]["active_normal"][str(W)]})
    print(json.dumps({"cells_with_active_rail": len(keys)}), flush=True)
    cells = {}; t0 = time.monotonic(); fails = []
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(query, k, inv): k for k in keys}
        for i, fu in enumerate(as_completed(futs), 1):
            try:
                k, v = fu.result(); cells[k] = v
            except Exception as e:  # noqa: BLE001
                fails.append([list(futs[fu]), str(e)[:200]])
            if i % 100 == 0:
                print(json.dumps({"done": i, "of": len(keys), "fail": len(fails), "s": round(time.monotonic() - t0)}), flush=True)
    if fails:
        raise ValueError(f"{len(fails)} cell queries failed, e.g. {fails[:3]}")
    trunc = sum(v["species_truncated"] or v["observers_truncated"] for v in cells.values())
    grid = [{"source_equal_area_grid_ix": k[0], "source_equal_area_grid_iy": k[1], "observations": v["records"],
             "species": len(v["species"]), "observers": v["observers"], "truncated": v["species_truncated"] or v["observers_truncated"]}
            for k, v in sorted(cells.items()) if v["records"] > 0]
    ok = {k: v for k, v in cells.items()
          if v["records"] >= GATE[0] and v["observers"] >= GATE[1] and len(v["species"]) >= GATE[2]}
    pairs = []
    for (i, j), v in sorted(ok.items()):
        for dx, dy in ((1, 0), (0, 1)):
            o = ok.get((i + dx, j + dy))
            if o is None:
                continue
            A, B = set(v["species"]), set(o["species"]); sh = len(A & B); un = len(A | B)
            pairs.append({"a": [i, j], "b": [i + dx, j + dy], "full_a": [v["records"], v["observers"], len(A)],
                          "full_b": [o["records"], o["observers"], len(B)],
                          "species_full": {"shared": sh, "union": un, "jaccard": sh / un}})
    out = {"version": "railstatus-eBird-GBIF-EOD-200km-2023-2025-v1", "dataset_key": EOD, "dataset_title": ds["title"],
           "records_by_year": per_year, "years": YEARS, "width_m": W, "gate_obs_observers_species": list(GATE),
           "n_cells_queried": len(keys), "n_cells_with_records": len(grid), "n_cells_truncated": trunc,
           "n_gated_cells": len(ok), "grid": grid, "pairs": pairs,
           "note": "GBIF-mediated EOD; recordedBy facet counts distinct eBird observer ids; 2023-2025 pooled"}
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out))
    print(json.dumps({k: v for k, v in out.items() if k not in ("grid", "pairs")} | {"n_pairs": len(pairs)}), flush=True)


if __name__ == "__main__":
    main()
