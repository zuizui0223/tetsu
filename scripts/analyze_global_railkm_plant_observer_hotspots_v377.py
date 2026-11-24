#!/usr/bin/env python3
"""v377 globally original rail-km macro vs observed plant richness WITHOUT observer hotspots.

This is a SOURCE-GATED, post-hoc robustness control for v371/v375 after
seeing original worldwide real rail map prediction gain. Before inspecting
the sensitivity, drop only the original 2025 high-observer-count top 5%
and top 10% of eligible source cells PER CONTINENT; keep all other recorded
species and the true original rail geography unchanged. Never select cells
by species outcome or rail density. Full 13.7GB worldwide botanical source
is already QA-authenticated by v369; original 2.237M km railway by v370.

These are whole landscape photo-recorded taxa, NOT railroad-field flora.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
from generate_global_equal_area_climate_land_grid_v351 import (
    sha_file,WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale
)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    load_json,verify_original_sources,source_climate_lookup,original_cell_rows,
    pair_rows,continent_holdout,CONTINENTS,GRID_WIDTHS,
    RICH_BASE,RICH_RAIL,BETA_BASE,BETA_RAIL,
    MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_BETA_PAIRS,
    MIN_TRAIN_BETA_PAIRS,MIN_CONTINENTS
)

SOURCE_NE_SHA="718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b"
FRACTIONS=(0.0,0.05,0.10)
CORRECTED_V371={
    100000:{"richness":0.004463745492835951,"jaccard":0.0001476673508208617},
    200000:{"richness":0.0066142438317320406,"jaccard":-0.00006120015031440586}
}

def select_unhot_cells(eligible,fraction):
    if fraction not in FRACTIONS:raise ValueError("Not a frozen source observer hotspot percentile")
    by_continent=defaultdict(list)
    for r in eligible:
        if r["continent"] not in CONTINENTS:raise ValueError("Original source continent unexpected")
        by_continent[r["continent"]].append(r)
    kept=[];qa={}
    for c in CONTINENTS:
        items=by_continent[c]
        # Sorted using ORIGINAL 2025 distinct observer count, not source
        # plant richness, real rail km, or predictor residual.
        sorted_obs=sorted(items,key=lambda r:(-int(r["observers"]),int(r["key"][0]),int(r["key"][1])))
        exclude=math.ceil(len(sorted_obs)*fraction)
        kept+=sorted_obs[exclude:]
        qa[c]={"original_effort_eligible_plant_cells":len(items),
               "source_observer_top_fraction_removed":fraction,
               "n_excluded_source_public_plant_observer_hotspot_cells":exclude,
               "n_retained_botanical_landcells":len(sorted_obs)-exclude,
               "source_plant_observer_count_minimum_retained":(
                    min(x["observers"] for x in sorted_obs[exclude:])
                    if len(sorted_obs)>exclude else None),
               "railkm_was_not_used_for_source_observer_hotspot_selection":True}
    if len(kept)!=len({x["key"] for x in kept}):
        raise AssertionError("Original source observer high-effort drop duplicated grid cells")
    return kept,qa

def source_models(cells,pair_source,original_rail):
    bycell={x["key"]:x for x in cells}
    pairs,pairqa=pair_rows(pair_source,bycell,original_rail)
    if len(bycell)!=len(cells):raise ValueError("Original source eligible cell IDs duplicated")
    r=continent_holdout(cells,RICH_BASE,RICH_RAIL,"log_species",False,
        MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2,
        "original_global_2025_landscape_recorded_plants_logrichness")
    b=continent_holdout(pairs,BETA_BASE,BETA_RAIL,"jaccard",True,
        MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,
        "source_2025_adjacent_landscape_Jaccard")
    return r,b,pairqa,len(pairs)

def run(plantpath,railpath,geojson,worldclim):
    plants=load_json(plantpath)
    rails=load_json(railpath)
    source=verify_original_sources(plants,rails)
    if sha_file(geojson)!=SOURCE_NE_SHA or sha_file(worldclim)!=WORLDCLIM_SHA:
        raise ValueError("Original source NaturalEarth or WorldClim source checksum changed")
    land,landqa=source_continents(geojson)
    raster,mem,rasqa=worldclim_tif(worldclim)
    try:
        scale,phys=climate_scale(raster)
        results=[]
        for width in GRID_WIDTHS:
            k=str(width)
            orig=plants["source_grids"][k]
            plantcells=orig["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            sourcekeys={(row["source_equal_area_grid_ix"],row["source_equal_area_grid_iy"])
                        for row in plantcells}
            climate,climqa=source_climate_lookup(land,raster,scale,sourcekeys,width)
            rawrail=rails["source_global_gridded_rail_exposure_km_by_cell_width"][k]
            eligible,sourceqa,rail_bykey=original_cell_rows(plantcells,rawrail,climate)
            original_j=orig["source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"]
            item={"original_global_grid_width_m":width,
               "original_source_valid_plant_landcells_before_hotspot_drop":len(eligible),
               "original_source_plant_species_and_2026_real_GERS_allworld_checksum_valid":True,
               "source_climate_QA":climqa,"source_botany_QA":sourceqa,
               "source_observer_hotspot_exclusions":[]}
            for fraction in FRACTIONS:
                current,remove_qa=select_unhot_cells(eligible,fraction)
                richness,beta,pairqa,npairs=source_models(current,original_j,rail_bykey)
                if fraction==0:
                    corrected=CORRECTED_V371[width]
                    for label,model in (("richness",richness),("jaccard",beta)):
                        delta=model["mean_holdout_continent_delta_MSE"]
                        if delta is None or abs(delta-corrected[label])>1e-9:
                            raise AssertionError("Original corrected v371 exact normalized worldwide source baseline changed")
                item["source_observer_hotspot_exclusions"].append({
                     "fraction_highest_original_source_distinct_observers_removed_per_continent":fraction,
                     "source_removed_original_observer_hotspot_count_by_continent":remove_qa,
                     "n_source_original_plant_cells_retained":len(current),
                     "n_source_adjacent_plant_Jaccard_pairs_retained":npairs,
                     "source_pair_eligibility":pairqa,
                     "source_2025_plant_recorded_richness_railkm_continent_holdout":richness,
                     "source_2025_adjacent_plant_Jaccard_railkm_continent_holdout":beta})
            results.append(item)
    finally:raster.close();mem.close()
    return {"version":"v377-original-fullworld-railkm-vs-photo-plant-recording-effort-hotspot-exclusion",
       "independent_source_authentication":source,
       "source_original_WorldClim":rasqa,"source_original_NaturalEarth":landqa,
       "original_SOURCE_freeze_fractions_no_outcome_selection":list(FRACTIONS),
       "source_total_Global_photo_linked_plant_2025_record_count":source[
         "original_2025_full_world_plant_detections_in_grid"],
       "global_100km_200km_original_source_results":results,
       "caveats":[
         "Post-hoc source observer-hotspot exclusion after seeing v371/v375 globally positive regional photo-richness predictive lift; sensitivity not confirmatory new ecological field finding",
         "Original iNaturalist public plant observer numbers are a detection effort proxy; species counts are whole landscapes, not true sampled railway plots",
         "No exact population/road/built-up adjustment here; source v376c separate global urbanization-adjusted test",
         "No source rail GIS line route operation, service frequency, vegetation microsite, direct seed vector or global thermal traits inferred",
         "Origin input true original per-cell species counts and neighboring Jaccard do not publish individual observer/taxon GPS or field sample matrices",
         "Repeated all-world 100km and 200km are not independent experiments"]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-v369-global-plants",required=True)
    p.add_argument("--original-v370-global-rail",required=True)
    p.add_argument("--natural-earth",required=True)
    p.add_argument("--worldclim-zip",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=run(a.original_v369_global_plants,a.original_v370_global_rail,
            a.natural_earth,a.worldclim_zip)
    path=Path(a.output);path.parent.mkdir(exist_ok=True,parents=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf8")
    print(json.dumps({"version":out["version"],"original_source_results":[
      {"width":r["original_global_grid_width_m"],
       "hotspot_drop":[
         {"fraction":q["fraction_highest_original_source_distinct_observers_removed_per_continent"],
          "cells":q["n_source_original_plant_cells_retained"],
          "pairs":q["n_source_adjacent_plant_Jaccard_pairs_retained"],
          "rich_delta":q["source_2025_plant_recorded_richness_railkm_continent_holdout"][
             "mean_holdout_continent_delta_MSE"],
          "rich_positive_continent_folds":q["source_2025_plant_recorded_richness_railkm_continent_holdout"][
             "n_holdout_continents_with_positive_delta_MSE"],
          "rich_n_continent_folds":q["source_2025_plant_recorded_richness_railkm_continent_holdout"][
             "n_original_continents_eligible_as_heldout_folds"],
          "jacc_delta":q["source_2025_adjacent_plant_Jaccard_railkm_continent_holdout"][
             "mean_holdout_continent_delta_MSE"],
          "jacc_positive_continent_folds":q["source_2025_adjacent_plant_Jaccard_railkm_continent_holdout"][
             "n_holdout_continents_with_positive_delta_MSE"],
          "jacc_n_continent_folds":q["source_2025_adjacent_plant_Jaccard_railkm_continent_holdout"][
             "n_original_continents_eligible_as_heldout_folds"]} for q in r["source_observer_hotspot_exclusions"]]}
       for r in out["global_100km_200km_original_source_results"]]},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
