#!/usr/bin/env python3
"""railstatus analysis: does the v371-v382 'rail km predicts landscape plant
record richness' signal come from operating rail, from the ballast corridor
itself (abandoned/disused lines), or from urban access (subway, tram)?

Same frozen v382 panel (v369 plants, v370 rail, GHSL built-up, EarthEnv tree,
GRIP4 roads, WorldClim, Natural Earth), same continent-heldout ridge with
1000 km block weights. Adds four rail-status predictors from the railstatus
aggregate, checks that active_normal reproduces v370 cell by cell, and runs a
climate/latitude-band shuffle placebo per category.

Interpretation key (pre-stated):
  inactive adds after active  -> corridor habitat / legacy, not trains
  active adds after inactive  -> something tied to operating lines
  subway or tram predicts     -> urban/access signal, not ballast habitat
All responses are whole-landscape iNaturalist plant records, not trackside flora.
"""
import argparse, hashlib, json, math, random
from collections import Counter
from pathlib import Path
import numpy as np

from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    verify_original_sources, source_climate_lookup, original_cell_rows,
    continent_holdout, GRID_WIDTHS, CONTINENTS, MIN_CONTINENT_RICH_CELLS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import (
    read_source_builtup, inject_built, ORIGINAL_NE_SHA)
from generate_global_equal_area_climate_land_grid_v351 import (
    WORLDCLIM_SHA, source_continents, worldclim_tif, climate_scale, sha_file)
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
from analyze_global_rail_road_tree_richness_v382 import (
    read, make_lookup, enrich_cells, within_continent_direction, BASE, LAND_THRESHOLDS)

CATS = ("active_normal", "inactive_normal", "subway", "tram_lightrail")
LOG = {c: "log_" + c for c in CATS}
N_PLACEBO = 39
SEED = "tetsu-railstatus-20261125-band-shuffle"
ROAD = ("log_road_density",)


def inject_status(rows, status, width):
    out = []
    lookups = {c: {(x["source_equal_area_grid_ix"], x["source_equal_area_grid_iy"]): float(x["km"])
                   for x in status["grids_km_by_category_and_width"][c][str(width)]} for c in CATS}
    worst = 0.0
    for r in rows:
        r = dict(r)
        for c in CATS:
            km = lookups[c].get(r["key"], 0.0)
            r[c + "_km"] = km
            r[LOG[c]] = math.log1p(km)
        worst = max(worst, abs(r["active_normal_km"] - r["rail_km"]))
        out.append(r)
    if worst > 0.01:
        raise ValueError(f"active_normal does not reproduce v370 cell km (max diff {worst} km)")
    return out, worst


def band(r):
    t = r["bio1"]
    tb = 0 if t < 8 else (1 if t < 20 else 2)
    return (r["continent"], tb, int(r["abs_lat"] // 10))


def shuffled(rows, field, rng):
    groups = {}
    for i, r in enumerate(rows):
        groups.setdefault(band(r), []).append(i)
    vals = [r[field] for r in rows]
    new = list(vals)
    for idx in groups.values():
        v = [vals[i] for i in idx]
        rng.shuffle(v)
        for i, x in zip(idx, v):
            new[i] = x
    return [dict(r, **{field: x}) for r, x in zip(rows, new)]


def ho(rows, base, ext, label):
    res = continent_holdout(rows, base, ext, "log_species", False,
                            MIN_CONTINENT_RICH_CELLS, MIN_CONTINENT_RICH_CELLS * 2, label)
    return {"mean_dMSE": res["mean_holdout_continent_delta_MSE"],
            "n_positive": res["n_holdout_continents_with_positive_delta_MSE"],
            "n_folds": res["n_original_continents_eligible_as_heldout_folds"],
            "per_continent": {x["heldout_continent"]: x["delta_MSE_positive_means_rail_added_prediction_helped"]
                              for x in res["source_independent_continent_holdouts"]},
            "base_mean_MSE": (sum(x["baseline_test_blockweighted_MSE"] for x in res["source_independent_continent_holdouts"])
                              / max(1, len(res["source_independent_continent_holdouts"])))}


def describe(rows):
    d = {}
    for c in CATS:
        pos = [r for r in rows if r[c + "_km"] > 0]
        d[c] = {"cells_positive": len(pos),
                "by_continent": dict(Counter(r["continent"] for r in pos)),
                "median_km_if_positive": float(np.median([r[c + "_km"] for r in pos])) if pos else None}
    a = np.array([r[LOG["active_normal"]] for r in rows])
    for c in CATS[1:]:
        b = np.array([r[LOG[c]] for r in rows])
        d[c]["pearson_log_with_active"] = float(np.corrcoef(a, b)[0, 1]) if b.std() > 0 else None
    return d


def main():
    p = argparse.ArgumentParser()
    for n in ("plant", "rail", "status", "built", "tree", "road", "ne", "climate", "out"):
        p.add_argument("--" + n, type=Path, required=True)
    a = p.parse_args()
    plant, rail, status, tree, road = map(read, (a.plant, a.rail, a.status, a.tree, a.road))
    checks = verify_original_sources(plant, rail)
    if tree["source_file_SHA256"] != EXPECTED_SHA:
        raise ValueError("Tree source not v380 verified")
    if road["source_zip_MD5"] != MD5:
        raise ValueError("Road source not GRIP4 verified")
    if sha_file(a.ne) != ORIGINAL_NE_SHA or sha_file(a.climate) != WORLDCLIM_SHA:
        raise ValueError("NE / WorldClim mismatch")
    if not status["version"].startswith("railstatus-v1"):
        raise ValueError("Not railstatus aggregate")
    built, _ = read_source_builtup(a.built)
    land, _ = source_continents(a.ne)
    raster, mem, _ = worldclim_tif(a.climate)
    rng = random.Random(hashlib.sha256(SEED.encode()).hexdigest())
    panels = []
    try:
        scale, _ = climate_scale(raster)
        for width in GRID_WIDTHS:
            g = plant["source_grids"][str(width)]
            original = g["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            keys = {(x["source_equal_area_grid_ix"], x["source_equal_area_grid_iy"]) for x in original}
            clim, _ = source_climate_lookup(land, raster, scale, keys, width)
            rows, _, _ = original_cell_rows(original, rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)], clim)
            rows, _ = inject_built(rows, built[str(width)])
            rows, worst = inject_status(rows, status, width)
            tr = make_lookup(tree["source_global_100km_200km_tree_consensus"][str(width)])
            rd = make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(width)])
            for minland in LAND_THRESHOLDS:
                m, qa = enrich_cells(rows, rd, tr, width, minland)
                RR = BASE + ROAD
                L = LOG
                C = {}
                for c in CATS:
                    C["noroad_" + c] = ho(m, BASE, BASE + (L[c],), c)
                    C["road_" + c] = ho(m, RR, RR + (L[c],), c)
                C["road_active_then_inactive"] = ho(m, RR + (L["active_normal"],), RR + (L["active_normal"], L["inactive_normal"]), "a->i")
                C["road_inactive_then_active"] = ho(m, RR + (L["inactive_normal"],), RR + (L["inactive_normal"], L["active_normal"]), "i->a")
                C["road_active_then_subway"] = ho(m, RR + (L["active_normal"],), RR + (L["active_normal"], L["subway"]), "a->s")
                C["road_active_then_tram"] = ho(m, RR + (L["active_normal"],), RR + (L["active_normal"], L["tram_lightrail"]), "a->t")
                C["road_subway_tram_then_active"] = ho(m, RR + (L["subway"], L["tram_lightrail"]),
                                                       RR + (L["subway"], L["tram_lightrail"], L["active_normal"]), "st->a")
                full = RR + tuple(L[c] for c in CATS)
                direction = {c: within_continent_direction(m, "log_species", full, L[c]) for c in CATS}
                placebo = {}
                for c in CATS:
                    real = C["road_" + c]["mean_dMSE"]
                    null = []
                    for _ in range(N_PLACEBO):
                        sh = shuffled(m, L[c], rng)
                        null.append(ho(sh, RR, RR + (L[c],), "placebo")["mean_dMSE"])
                    placebo[c] = {"real_mean_dMSE": real, "null_mean": float(np.mean(null)),
                                  "null_min": float(np.min(null)), "null_max": float(np.max(null)),
                                  "n_null_ge_real": int(sum(x >= real for x in null)), "n_null": N_PLACEBO}
                cond = {}
                for c in ("inactive_normal",):
                    real = C["road_active_then_inactive"]["mean_dMSE"]
                    null = []
                    for _ in range(N_PLACEBO):
                        sh = shuffled(m, L[c], rng)
                        null.append(ho(sh, RR + (L["active_normal"],), RR + (L["active_normal"], L[c]), "placebo")["mean_dMSE"])
                    cond["inactive_after_active"] = {"real_mean_dMSE": real, "null_mean": float(np.mean(null)),
                                                     "null_max": float(np.max(null)),
                                                     "n_null_ge_real": int(sum(x >= real for x in null)), "n_null": N_PLACEBO}
                panel = {"width_m": width, "min_land_fraction": minland, "n_cells": len(m),
                         "active_vs_v370_max_abs_km_diff": worst, "excluded": qa,
                         "describe": describe(m), "contrasts": C, "within_continent_direction_full_model": direction,
                         "band_shuffle_placebo_road_plus_category": placebo,
                         "band_shuffle_placebo_conditional": cond}
                panels.append(panel)
                print(json.dumps({"width": width, "minland": minland, "n": len(m),
                                  "dMSE": {k: (round(v["mean_dMSE"], 6) if v["mean_dMSE"] is not None else None, v["n_positive"], v["n_folds"]) for k, v in C.items()},
                                  "placebo": {k: (v["n_null_ge_real"], v["n_null"]) for k, v in placebo.items()},
                                  "cond": {k: (v["n_null_ge_real"], v["n_null"]) for k, v in cond.items()}}), flush=True)
    finally:
        raster.close(); mem.close()
    out = {"version": "railstatus-v1-analysis", "source_replay": checks,
           "status_totals_km": status["total_km"], "panels": panels,
           "scope": "whole-landscape 2025 iNaturalist plant record richness; NOT trackside flora; not causal",
           "limitations": ["OSM abandoned/disused tagging is incomplete and uneven by country",
                           "subway km includes underground track; tram/light rail mostly street-embedded",
                           "same panel and model family as v382; continent folds share training data",
                           "placebo shuffles values within continent x temperature band x 10-degree latitude band"]}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, allow_nan=False))


if __name__ == "__main__":
    main()
