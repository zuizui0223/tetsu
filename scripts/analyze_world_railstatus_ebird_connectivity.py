#!/usr/bin/env python3
"""Independent observation-network replication of the iNaturalist bird
crossing-similarity association, using eBird (GBIF EOD) 2023-2025 pooled.

Pre-stated primary (fixed before seeing eBird outcomes): 200 km, 50% land,
endpoint-FE Jaccard effect of an active-rail boundary crossing, same pair
design, covariates, gates and block-clustered intervals as the iNaturalist
bird analysis (direct-boundary / observer-split pipelines). Secondary: block
FE, Simpson, inactive crossing, active-minus-inactive, 80% land, continent-
heldout prediction. Replication = same sign with interval excluding zero.
"""
import argparse, json
from pathlib import Path
from analyze_world_railstatus_direct_boundary_beta import (
    build_pairs, fit, PRED, DIFF, CROSS, ACTIVE, CATS, join_population, inject_status, LAND_THRESHOLDS,
    read, make_lookup, enrich_cells, read_source_builtup, ORIGINAL_NE_SHA, WORLDCLIM_SHA, EXPECTED_SHA, MD5,
    sha_file, source_continents, worldclim_tif, climate_scale, source_climate_lookup, original_cell_rows,
    inject_built, continent_holdout, POP_SHA)
from analyze_world_railstatus_observer_split_connectivity import adapt_pair
from analyze_world_railstatus_cross_taxa import convert

W = 200_000
ENDPOINT = DIFF + ("horizontal",) + CROSS


def summarize(f):
    if f.get("status") != "estimated":
        return {"status": f.get("status"), "n_pairs": f.get("n_pairs"), "n_blocks": f.get("n_blocks"),
                "n_blocks_residual_active": f.get("n_blocks_residual_active")}
    return {"status": "estimated", "n_pairs": f["n_pairs"], "n_blocks": f["n_blocks"],
            "n_blocks_residual_active": f["n_blocks_residual_active"], "n_cross_active": f["n_cross_active"],
            "active": f["effects"][ACTIVE], "inactive": f["effects"]["cross_inactive_normal"],
            "active_minus_inactive": f["active_minus_inactive"]}


def main():
    p = argparse.ArgumentParser()
    for n in ("ebird", "rail", "status", "edges", "built", "tree", "road", "population", "ne", "climate", "out"):
        p.add_argument("--" + n, type=Path, required=True)
    a = p.parse_args()
    eb, rail, status, edge, tree, road, pop = map(read, (a.ebird, a.rail, a.status, a.edges, a.tree, a.road, a.population))
    if eb["version"] != "railstatus-eBird-GBIF-EOD-200km-2023-2025-v1":
        raise ValueError("wrong eBird source")
    if edge["version"] != "railstatus-direct-grid-boundary-v1" or not edge["source_verified"]:
        raise ValueError("edge source")
    if any(abs(edge["km"][c] - status["total_km"][c]) > .01 for c in CATS):
        raise ValueError("status km mismatch")
    if tree["source_file_SHA256"] != EXPECTED_SHA or road["source_zip_MD5"] != MD5 or pop["source_archive_SHA256"] != POP_SHA:
        raise ValueError("frozen covariate mismatch")
    if sha_file(a.ne) != ORIGINAL_NE_SHA or sha_file(a.climate) != WORLDCLIM_SHA:
        raise ValueError("NE/WorldClim mismatch")
    built, _ = read_source_builtup(a.built); land, _ = source_continents(a.ne); raster, mem, _ = worldclim_tif(a.climate)
    panels = []
    try:
        scale, _ = climate_scale(raster)
        keys = {(z["source_equal_area_grid_ix"], z["source_equal_area_grid_iy"]) for z in eb["grid"]}
        clim, _ = source_climate_lookup(land, raster, scale, keys, W)
        tr = make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
        rd = make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
        pp = make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
        el = {c: {(tuple(z["a"]), tuple(z["b"])): z["crossing_source_segments"] for z in edge["edges"][str(W)][c]} for c in CATS}
        rows, _, _ = original_cell_rows([convert(z) for z in eb["grid"]], rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)], clim)
        rows, _ = inject_built(rows, built[str(W)]); rows, _ = inject_status(rows, status, W)
        for th in LAND_THRESHOLDS:
            m, _ = enrich_cells(rows, rd, tr, W, th); m, _ = join_population(m, pp, rd)
            pr = build_pairs([adapt_pair(z) for z in eb["pairs"]], m, el)
            res = {"width_m": W, "min_land_fraction": th, "n_pairs": len(pr),
                   "n_pairs_cross_active": sum(r[ACTIVE] > 0 for r in pr)}
            for t in ("jaccard", "simpson"):
                res[t] = {"endpoint_FE": summarize(fit(pr, t, True, ENDPOINT)), "block_FE": summarize(fit(pr, t, False))}
            base = tuple(k for k in PRED if k != ACTIVE)
            cv = continent_holdout(pr, base, base + (ACTIVE,), "jaccard", True, 50, 300, "eBird Jaccard") if pr else None
            res["jaccard_continent_holdout"] = None if cv is None else {
                "mean_dMSE": cv["mean_holdout_continent_delta_MSE"],
                "n_positive": cv["n_holdout_continents_with_positive_delta_MSE"],
                "n_folds": cv["n_original_continents_eligible_as_heldout_folds"]}
            panels.append(res)
            print(json.dumps({"ebird_panel": res}), flush=True)
    finally:
        raster.close(); mem.close()
    out = {"version": "railstatus-eBird-connectivity-v1", "source": {k: eb[k] for k in ("dataset_title", "records_by_year", "years", "n_cells_queried", "n_cells_truncated", "n_gated_cells")},
           "n_source_pairs": len(eb["pairs"]), "panels": panels,
           "primary": "200km land50 endpoint-FE Jaccard active crossing effect",
           "limits": ["eBird and iNaturalist observers overlap partly", "pooled 2023-2025, not annual",
                      "GBIF EOD excludes sensitive species and some records", "mapped crossing is not train traffic"]}
    a.out.parent.mkdir(parents=True, exist_ok=True); a.out.write_text(json.dumps(out, allow_nan=False))


if __name__ == "__main__":
    main()
