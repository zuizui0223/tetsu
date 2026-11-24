#!/usr/bin/env python3
"""v375 GLOBAL observed plant taxonomic richness and landscape Jaccard:
conditional railway GIS placement negative control.

Use SOURCE-AUTHENTIC actual v369 full worldwide 2025 public plant source
and v370 actual worldwide original 2.237million GERS railway km.
Reproduce exact v371 six-continent heldout improvement before testing.

Keep original species richness, adjacent whole-landscape Jaccard, effort,
geo/BIO1, latitudes and blocked heldout folds unchanged. Reassign the original
railway km INCLUDING true zeros among eligible plant grid cells ONLY inside
continent x annual macrotemperature x abs latitude 10deg strata. This is
NOT railway geospatial process randomization or direct railway vegetation.
Do not redistribute any cell-level original plant species or observer data.
"""
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    load_json,verify_original_sources,
    source_climate_lookup,original_cell_rows,pair_rows,continent_holdout,
    MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_BETA_PAIRS,
    MIN_TRAIN_BETA_PAIRS,MIN_CONTINENTS,
    RICH_BASE,RICH_RAIL,BETA_BASE,BETA_RAIL,GRID_WIDTHS
)
from generate_global_equal_area_climate_land_grid_v351 import (
    source_continents,worldclim_tif,climate_scale,sha_file,WORLDCLIM_SHA,temp_band
)

SEED="tetsu-v375-20261124-conditional-rail-map-placebo"
DRAWS=39
SOURCE_NE_SHA="718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b"
EXPECTED={
  100000:{"richness":0.0034565261140192885,
          "jaccard":0.00003170576215421156},
  200000:{"richness":0.0036532861030028514,
          "jaccard":-0.00010563827123567526}
}
BANDS=("cool","temperate","warm")

def shuffle_stratum(row):
    if row["continent"] not in (
         "Africa","Asia","Europe","North America","South America","Oceania"):
        raise ValueError("Original unknown continent in matched landscape flora")
    t=temp_band(row["bio1"])
    if t not in BANDS:raise ValueError("Source BIO1 invalid")
    lat_band=int(row["abs_lat"]//10)
    return (row["continent"],t,lat_band)

def permute_railkm_in_eligible_observed_cells(rows,rng):
    """Keep every source plant observation and actual original species fixed.

    Condition on whole-continent, BIO1 macroclimate, |latitude| 10degree.
    The source botanical sampling set is fixed; permutations do not
    infer railway absence in low-effort or missing botanical cells.
    """
    before=[float(r["rail_km"]) for r in rows]
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        if r["rail_km"]<0 or not math.isfinite(r["rail_km"]):
            raise ValueError("Original world GERS railway source negative or missing")
        groups[shuffle_stratum(r)].append(i)
    candidate_groups=[idx for idx in groups.values() if len(idx)>=3
                      and len({before[i] for i in idx})>=2]
    after=before[:]
    for indices in candidate_groups:
        vals=[before[i] for i in indices]
        rng.shuffle(vals)
        for i,x in zip(indices,vals):
            after[i]=x
    # all-world and per-stratum actual source exposure and zeros preserved.
    if not math.isclose(sum(before),sum(after),rel_tol=0,abs_tol=1e-6):
        raise AssertionError("Conditional rail-map placebo lost eligible source railway km")
    for indices in groups.values():
        if sorted(before[i] for i in indices)!=sorted(after[i] for i in indices):
            raise AssertionError("Original climate-region rail km source multiset changed")
    for r,km in zip(rows,after):
        r["rail_km"]=km
        r["log_rail_km"]=math.log1p(km)
    return {"original_source_climate_latitude_strata":len(groups),
        "strata_with_n_at_least_three_and_nontrivial_rail_distribution":len(candidate_groups),
        "source_valid_botanical_cells_railkm_movable":sum(len(s) for s in candidate_groups),
        "actual_source_botanical_cells_with_different_assigned_rail_km":
            sum(a!=b for a,b in zip(before,after)),
        "original_total_GERS_km_in_observation_eligible_source_cells":sum(before),
        "source_railkm_is_reassigned_predictor_only_not_original_GERS_line_geometry":True}

def compute_holdouts(cells,pair_source,rail_by_key):
    bycell={z["key"]:z for z in cells}
    pairs,qa=pair_rows(pair_source,bycell,rail_by_key)
    rich=continent_holdout(cells,RICH_BASE,RICH_RAIL,"log_species",False,
             MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2,"source_world_richness")
    beta=continent_holdout(pairs,BETA_BASE,BETA_RAIL,"jaccard",True,
             MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,"source_world_jaccard")
    if (rich["n_original_continents_eligible_as_heldout_folds"]<MIN_CONTINENTS or
        beta["n_original_continents_eligible_as_heldout_folds"]<MIN_CONTINENTS):
        raise ValueError("Original full world cannot support valid external continent holdouts")
    return {"richness":rich,"jaccard":beta,"n_beta_original_source":len(pairs),
            "pair_QA":qa}

def effect(model):
    return model["mean_holdout_continent_delta_MSE"]

def tech_summary(observed,values):
    if not values:return None
    a=sorted(values)
    n=len(a)
    return {"observed_real_rail_map_delta_MSE":observed,
            "n_fake_rail_map_source_placebo_draws":n,
            "fake_map_null_mean_delta_MSE":mean(a),
            "fake_map_null_median_delta_MSE":median(a),
            "fake_map_null_low_high":[a[0],a[-1]],
            "fake_map_null_technical_95percentile":[a[int(.025*(n-1))],a[int(.975*(n-1))]],
            "source_placement_fraction_fake_at_least_as_predictive_as_real_plus_one":
                (sum(x>=observed-1e-12 for x in values)+1)/(n+1),
            "fake_rail_maps_with_positive_holdout_prediction_lift":sum(x>0 for x in values),
            "exploratory_technical_negative_control_not_CI_or_causal_randomization":True}

def run(plantpath,railpath,nepath,wcpath):
    source_plants=load_json(plantpath)
    source_rails=load_json(railpath)
    integrity=verify_original_sources(source_plants,source_rails)
    if sha_file(nepath)!=SOURCE_NE_SHA or sha_file(wcpath)!=WORLDCLIM_SHA:
        raise ValueError("Source exact original NE/WorldClim source hashes changed")
    land,landqa=source_continents(nepath)
    raster,mem,rasqa=worldclim_tif(wcpath)
    panels=[]
    try:
        scale,scaleqa=climate_scale(raster)
        for width in GRID_WIDTHS:
            original=source_plants["source_grids"][str(width)]
            plantgrid=original["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            source_pairs=original["source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"]
            keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for r in plantgrid}
            climate,climateqa=source_climate_lookup(land,raster,scale,keys,width)
            raw_rail=source_rails["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)]
            cells,cellqa,all_rail=original_cell_rows(plantgrid,raw_rail,climate)
            existing=compute_holdouts(cells,source_pairs,all_rail)
            baseline={name:effect(existing[name]) for name in ("richness","jaccard")}
            for typ,measured in baseline.items():
                if abs(measured-EXPECTED[width][typ])>1e-9:
                    raise AssertionError(
                       f"Original v371 genuine 2025 source M0/M1 continent-heldout result changed {width} {typ}: {measured}")
            original_rail=[r["rail_km"] for r in cells]
            rng=random.Random(hashlib.sha256((SEED+"|"+str(width)).encode()).hexdigest())
            null_values={name:[] for name in ("richness","jaccard")}
            placebos=[]
            for i in range(DRAWS):
                changed=permute_railkm_in_eligible_observed_cells(cells,rng)
                fitted=compute_holdouts(cells,source_pairs,all_rail)
                for name in null_values:
                    null_values[name].append(effect(fitted[name]))
                placebos.append(changed)
                # Reset exact source rail exposure for each random draw,
                # no sequential accumulation of shuffles from prior draw.
                for z,km in zip(cells,original_rail):
                    z["rail_km"]=km
                    z["log_rail_km"]=math.log1p(km)
            if any(z["rail_km"]!=orig for z,orig in zip(cells,original_rail)):
                raise AssertionError("Original observed GERS source not restored after placebo draws")
            qa={
                "source_climate_region_latitude_conditioning":placebos[0],
                "n_placebo_maps":DRAWS,
                "min_changed_source_observation_cells":min(q["actual_source_botanical_cells_with_different_assigned_rail_km"] for q in placebos),
                "max_changed_source_observation_cells":max(q["actual_source_botanical_cells_with_different_assigned_rail_km"] for q in placebos),
                "source_observed_botanical_cells":len(cells),
                "original_observed_adjoining_plant_Jaccard_pairs":existing["n_beta_original_source"]
            }
            panels.append({
                "source_grid_width_m":width,
                "n_eligible_original_2025_public_botanical_cells":len(cells),
                "n_eligible_original_global_plant_neighbor_Jaccard_pairs":existing["n_beta_original_source"],
                "original_2025_source_actual_railkm_holdout_M0minusM1":baseline,
                "n_original_heldout_continents":{"richness":existing["richness"]["n_original_continents_eligible_as_heldout_folds"],
                                                "jaccard":existing["jaccard"]["n_original_continents_eligible_as_heldout_folds"]},
                "original_climate_and_source_observation_QA":climateqa,
                "source_rail_shuffling_QA":qa,
                "source_rail_map_spatial_placebo_diagnostics":{
                    name:tech_summary(baseline[name],null_values[name])
                    for name in null_values
                }})
    finally:
        raster.close();mem.close()
    return {"version":"v375-genuine-world-railkm-conditional-spatial-placebo-for-2025-public-plant-richer-landscapes",
        "world_plant_source_obs_gzip_SHA256":integrity["original_iNat_obs_gzip_SHA256"],
        "world_plant_source_taxon_gzip_SHA256":integrity["original_iNat_taxa_gzip_SHA256"],
        "original_GERS_population_all_128_source_2026_candidate_railkm":2237043.481,
        "source_world_plant_2025_public_observations_with_geo_and_research_quality":integrity["original_2025_full_world_plant_detections_in_grid"],
        "source_original_worldclim_sha256":WORLDCLIM_SHA,
        "source_original_NaturalEarth_sha256":SOURCE_NE_SHA,
        "placebo_frozen_seed_sha256":hashlib.sha256(SEED.encode()).hexdigest(),
        "source_technical_conditional_placebo_draws_each_grid":DRAWS,
        "source_results":panels,
        "interpretation_guardrails":[
            "The 39 placement shuffles are posthoc exploratory negative controls after seeing v371 observed global result, not pre-registered randomized intervention p values",
            "Placebo moves source rail km between OTHER QUALIFIED BOTANICAL CELLS of same continent/temperature band/abs latitude decile; original global railway geometry is NOT physically shuffled",
            "This controls broad bioclimatic and continental rail distributions, but does NOT control unmeasured urbanization, roads, population or iNaturalist observation accessibility",
            "All source observed plant records are full LANDSCAPE counts and neighbor Jaccard, never near-track or true ballast vegetation plots",
            "With known geoprivacy fields absent in the open-data plant metadata, no source botanical location is assigned to within 80 m of track",
            "Actual rail annual service, freight movement, direct seed carriage, climate range limits and establishment success not measured",
            "Source 100km/200km are two resolution views of the SAME original source population, not independent global ecological replications",
            "No original individual observer IDs, plant taxa or source rare-species coordinates/railway segment original IDs exported"]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-v369-global-plants",required=True)
    p.add_argument("--original-v370-global-rail",required=True)
    p.add_argument("--natural-earth",required=True)
    p.add_argument("--worldclim-zip",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    obj=run(a.original_v369_global_plants,a.original_v370_global_rail,
            a.natural_earth,a.worldclim_zip)
    out=Path(a.output);out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding="utf8")
    print(json.dumps({"version":obj["version"],"source_results":obj["source_results"]},
                     ensure_ascii=False),flush=True)
if __name__=="__main__":main()
