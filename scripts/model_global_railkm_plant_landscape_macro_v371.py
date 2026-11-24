#!/usr/bin/env python3
"""v371: GLOBAL source-complete macro prediction from original 2026 railway GIS,
2025 photo-linked plant community 100km/200km regions and source WorldClim.

Observed outcomes are ENTIRE regional photo-linked flora, not railway
ballast presence. No source original observer IDs, rare plant coordinates,
individual taxon names or original Overture geometries are redistributed.
Worldwide one-continent-out spatial validation; dependent source cell-pairs
never treated as independent field studies or randomized railway plots.
"""
from __future__ import annotations
import argparse,hashlib,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    OBS_ETAG,TAXA_ETAG,X_START,Y_START
)
from generate_global_equal_area_climate_land_grid_v351 import (
    WORLDCLIM_SHA,CONTINENTS,source_continents,worldclim_tif,
    climate_scale,read_temp
)

GRID_SIZES=("100000","200000")
RIDGE_LAMBDA=.10
MIN_BIO_OBS=50
MIN_BIO_OBSERVERS=5
MIN_BIO_SPECIES=20
MIN_CONTINENT_TEST_CELLS=100
MIN_CONTINENT_TEST_PAIRS=50
MIN_TRAIN_REGIONAL_CELLS=500
MIN_TRAIN_PAIRS=300
MIN_INDEPENDENT_CONTINENT_FOLDS=4
BLOCK_M=1_000_000
SOURCE_RAIL_KM=2237043.481
SOURCE_GRIDDING_TOLERANCE_KM=.08

RICH_BASE=("log_observations","log_observers","bio1_C","bio1_C2",
           "abs_latitude","abs_latitude_squared")
RICH_RAIL=("log_rail_km",)
BETA_BASE=("log_geometric_observations","abs_log_observation_contrast",
           "log_geometric_observers","abs_log_observer_contrast",
           "log_geometric_species","abs_log_species_contrast",
           "bio1_mean_C","abs_bio1_difference_C","mean_abs_latitude",
           "north_south_neighbor")
BETA_RAIL=("log_mean_rail_km","abs_log_rail_contrast")

def source_load(path):
    return json.loads(Path(path).read_text(encoding="utf8"))

def original_sources(plant_file,rail_file,worldclim,nefile):
    from shapely.geometry import Point
    from pyproj import Transformer
    flora=source_load(plant_file)
    rail=source_load(rail_file)
    pmeta=flora.get("source_summary") or {}
    if (pmeta.get("original_source_iNat_obs_S3_ETag")!=OBS_ETAG or
        pmeta.get("original_source_iNat_taxa_S3_ETag")!=TAXA_ETAG or
        not pmeta.get("original_coarse_grid_100km_and_200km_are_OBSERVATION_OPPORTUNITY_not_causal_rail_species_presence")
        or not pmeta.get("source_iNat_original_2026_10_27_observations_gzip_SHA256_full")):
        raise ValueError("Original full 13.7GB photo-linked source did not completely authenticate")
    if (rail.get("source_fullworld_candidate_normal_rail_km")!=SOURCE_RAIL_KM or
        not rail.get("source_original_fullworld_rail_128_files_SHA_verified") or
        "source_global_gridded_rail_exposure_km_by_cell_width" not in rail):
        raise ValueError("Original all-128 rail source grid incomplete")
    for width in GRID_SIZES:
        total=rail["source_100km_200km_grid_rail_km_total_within_mapped_cells"].get(width)
        if not isinstance(total,(float,int)) or abs(total-SOURCE_RAIL_KM)>SOURCE_GRIDDING_TOLERANCE_KM:
            raise ValueError("Original worldwide true railkm clipped grid source total invalid")
    checksum=hashlib.sha256(Path(worldclim).read_bytes()).hexdigest()
    if checksum!=WORLDCLIM_SHA:
        raise ValueError("Original WorldClim macroclimate source ZIP SHA changed")
    lands,neqa=source_continents(nefile)
    rast,mem,meta=worldclim_tif(worldclim)
    try:
        scale,phys=climate_scale(rast)
        transformer=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
        coords={}
        result={}
        qa=Counter()
        for width in GRID_SIZES:
            cell_m=int(width)
            source_cells=flora["source_grids"][width][
                "source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            original_pairs=flora["source_grids"][width][
                "source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"]
            exposed=rail["source_global_gridded_rail_exposure_km_by_cell_width"][width]
            rail_by_key={}
            for val in exposed:
                key=(val["source_equal_area_grid_ix"],val["source_equal_area_grid_iy"])
                if key in rail_by_key:raise ValueError("Duplicate global rail km grid cell")
                km=val["source_real_original_Overture_normal_rail_km_clipped_into_cell"]
                if not math.isfinite(km) or km<=0:
                    raise ValueError("Original gridded mapped railway km must be >0 in nonzero list")
                rail_by_key[key]=km
            plants={}
            for entry in source_cells:
                key=(entry["source_equal_area_grid_ix"],entry["source_equal_area_grid_iy"])
                if key in plants:raise ValueError("Duplicate original 2025 worldwide plant grid cell")
                nobs=entry["source_year_2025_photo_linked_plants_observation_count"]
                nob=entry["source_2025_distinct_observers_with_photo_linked_plant_records"]
                ntax=entry["source_2025_plant_species_richness"]
                if not all(isinstance(z,int) and z>=0 for z in (nobs,nob,ntax)):
                    raise ValueError("Invalid independent worldwide botanical source per-cell effort")
                if nobs<MIN_BIO_OBS or nob<MIN_BIO_OBSERVERS or ntax<MIN_BIO_SPECIES:
                    qa["source_plant_cells_with_insufficient_observation_effort_"+width]+=1
                    continue
                x=X_START+(key[0]+.5)*cell_m
                y=Y_START+(key[1]+.5)*cell_m
                lon,lat=transformer.transform(x,y)
                if not (math.isfinite(lon) and math.isfinite(lat)):
                    qa["source_botanical_cell_center_invalid_projection_"+width]+=1
                    continue
                matches=[c for c in CONTINENTS if lands[c].covers(Point(lon,lat))]
                if len(matches)!=1:
                    qa["source_botanical_cell_center_not_one_continent_land_"+width]+=1
                    continue
                temp=read_temp(rast,lat,lon)
                if temp is None:
                    qa["source_botanical_cell_without_original_WorldClim_BIO1_"+width]+=1
                    continue
                temp=temp/scale
                if not -70<=temp<=50:
                    raise ValueError("World source site macrotemperature outside plausible °C")
                railway_km=rail_by_key.get(key,0.0)
                if railway_km<0:raise ValueError("Complete GERS source no negative rail")
                record={"cell":key,"continent":matches[0],"source_grid_width_m":cell_m,
                    "n_photo_plant_observations":nobs,
                    "n_public_distinct_plant_observers":nob,
                    "n_observed_plant_species":ntax,
                    "source_original_bio1_C":temp,
                    "abs_lat":abs(lat),
                    "rail_km":railway_km,
                    "original_equal_area_1000km_block":(int(x//BLOCK_M),int(y//BLOCK_M))}
                plants[key]=record
            betas=[]
            for pair in original_pairs:
                a=tuple(pair["source_cell_A_original_grid_ix_iy"])
                b=tuple(pair["source_cell_B_original_grid_ix_iy"])
                if a not in plants or b not in plants:
                    qa["original_neighbor_beta_pair_removed_due_to_source_cell_effort_climate_or_land_"+width]+=1
                    continue
                A=plants[a];B=plants[b]
                if A["continent"]!=B["continent"]:
                    qa["original_neighbor_beta_pair_crosses_continent_boundary_"+width]+=1
                    continue
                jaccard=pair["source_biological_Jaccard_similarity_between_ENTIRE_landscape_cells"]
                if not isinstance(jaccard,(float,int)) or not 0<=jaccard<=1:
                    raise ValueError("Original whole-world source Jaccard is not in 0..1")
                if abs(a[0]-b[0])+abs(a[1]-b[1])!=1:
                    raise ValueError("Original source neighbor Jaccard not adjacent or duplicate grid")
                if pair["source_grid_A_original_plant_species_n"]!=A["n_observed_plant_species"]:
                    raise ValueError("Original 2025 plant species richness mismatch across source grids and pairs")
                if pair["source_grid_B_original_plant_species_n"]!=B["n_observed_plant_species"]:
                    raise ValueError("Original 2025 plant species richness mismatch across source grids and pairs")
                x=(A["cell"][0]+B["cell"][0]+1)*cell_m/2+X_START
                y=(A["cell"][1]+B["cell"][1]+1)*cell_m/2+Y_START
                betas.append({"A":A,"B":B,"y":float(jaccard),
                    "continent":A["continent"],
                    "block":(int(x//BLOCK_M),int(y//BLOCK_M)),
                    "direction":int(a[0]==b[0])})
            result[width]={"cells":list(plants.values()),"pairs":betas,
                "eligible_rail_geographic_grid_original_source_cell_n":len(rail_by_key),
                "source_botanical_input_cells_n":len(source_cells),
                "source_original_neighbor_jaccard_input_pair_n":len(original_pairs)}
        return result,dict(qa),{
            "original_WorldClim_BIO1_SHA":WORLDCLIM_SHA,
            "original_WorldClim_BIO1_units_divisor":scale,
            "source_original_NaturalEarth_world_admin0_source_countries_QA":neqa,
            "source_original_botanical_full_obs_hash":pmeta[
                "source_iNat_original_2026_10_27_observations_gzip_SHA256_full"]}
    finally:
        rast.close();mem.close()

def richness_features(r):
    t=r["source_original_bio1_C"]
    lat=r["abs_lat"]
    return {
      "log_observations":math.log1p(r["n_photo_plant_observations"]),
      "log_observers":math.log1p(r["n_public_distinct_plant_observers"]),
      "bio1_C":t,"bio1_C2":t*t,
      "abs_latitude":lat,"abs_latitude_squared":lat*lat,
      "log_rail_km":math.log1p(r["rail_km"])}

def beta_features(pair):
    a,b=pair["A"],pair["B"]
    def both(name):
        x=math.log1p(a[name]);y=math.log1p(b[name])
        return .5*(x+y),abs(x-y)
    ob=both("n_photo_plant_observations")
    obs=both("n_public_distinct_plant_observers")
    rich=both("n_observed_plant_species")
    rails=both("rail_km")
    t1,t2=a["source_original_bio1_C"],b["source_original_bio1_C"]
    return {
       "log_geometric_observations":ob[0],
       "abs_log_observation_contrast":ob[1],
       "log_geometric_observers":obs[0],
       "abs_log_observer_contrast":obs[1],
       "log_geometric_species":rich[0],
       "abs_log_species_contrast":rich[1],
       "bio1_mean_C":.5*(t1+t2),
       "abs_bio1_difference_C":abs(t1-t2),
       "mean_abs_latitude":.5*(a["abs_lat"]+b["abs_lat"]),
       "north_south_neighbor":float(pair["direction"]),
       "log_mean_rail_km":math.log1p(.5*(a["rail_km"]+b["rail_km"])),
       "abs_log_rail_contrast":rails[1]}

def original_source_block_key(entry):
    # True original entire-world source richness CELL records preserve the
    # field name original_equal_area_1000km_block, while adjacent beta
    # PAIRS use block. They are the same 1000km equal-area source block
    # semantics and must not be treated as two independent studies.
    if "block" in entry:
        return tuple(entry["block"])
    if "original_equal_area_1000km_block" in entry:
        return tuple(entry["original_equal_area_1000km_block"])
    raise ValueError("Genuine global plant source cell lacks geographic block")

def block_weights(items):
    count=Counter(original_source_block_key(x) for x in items)
    if not count:raise ValueError("No original global source spatial blocks")
    weights=np.array([1.0/count[original_source_block_key(x)] for x in items],dtype=float)
    return weights/weights.sum()

def ridge_estimate(train,test,base,extra,target,kind):
    fields=tuple(base)+tuple(extra)
    f=richness_features if kind=="richness" else beta_features
    x=np.asarray([[f(z)[k] for k in fields] for z in train],dtype=np.float64)
    xt=np.asarray([[f(z)[k] for k in fields] for z in test],dtype=np.float64)
    y=np.asarray([target(z) for z in train],dtype=float)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(xt)) or not np.all(np.isfinite(y)):
        raise ValueError("Original global 100km macro model contains NaN")
    w=block_weights(train)
    mu=(w[:,None]*x).sum(axis=0)
    sd=np.sqrt((w[:,None]*(x-mu)**2).sum(axis=0))
    sd=np.where(sd<1e-12,1.0,sd)
    X=np.column_stack([np.ones(len(train)),(x-mu)/sd])
    Xt=np.column_stack([np.ones(len(test)),(xt-mu)/sd])
    penal=np.eye(X.shape[1]);penal[0,0]=0
    coef=np.linalg.solve(X.T@(w[:,None]*X)+RIDGE_LAMBDA*penal,X.T@(w*y))
    p=Xt@coef
    if kind=="beta":p=np.clip(p,0,1)
    return p

def cv_country_holdout(items,kind):
    if kind not in ("richness","beta"):raise ValueError("Unknown global macroecology outcome")
    base=RICH_BASE if kind=="richness" else BETA_BASE
    extra=RICH_RAIL if kind=="richness" else BETA_RAIL
    def y(z):
        return math.log1p(z["n_observed_plant_species"]) if kind=="richness" else z["y"]
    required=MIN_CONTINENT_TEST_CELLS if kind=="richness" else MIN_CONTINENT_TEST_PAIRS
    trainminimum=MIN_TRAIN_REGIONAL_CELLS if kind=="richness" else MIN_TRAIN_PAIRS
    folds=[];support={}
    for continent in CONTINENTS:
        train=[z for z in items if z["continent"]!=continent]
        test=[z for z in items if z["continent"]==continent]
        support[continent]={"n_original_test_source_cells_or_neighbor_pairs":len(test),
            "n_actual_source_training_cells_or_neighbor_pairs":len(train)}
        if len(test)<required or len(train)<trainminimum:
            support[continent]["not_identified_low_original_source_sample"]=True
            continue
        expected=np.asarray([y(z) for z in test],dtype=float)
        p0=ridge_estimate(train,test,base,(),y,kind)
        p1=ridge_estimate(train,test,base,extra,y,kind)
        weights=block_weights(test)
        mse0=float(np.sum(weights*(expected-p0)**2))
        mse1=float(np.sum(weights*(expected-p1)**2))
        mae0=float(np.sum(weights*np.abs(expected-p0)))
        mae1=float(np.sum(weights*np.abs(expected-p1)))
        folds.append({"heldout_continent":continent,
          "n_test_original_macro_source_units_not_independent_field_studies":len(test),
          "n_training_original_macro_source_units":len(train),
          "n_independent_geographic_1000km_test_blocks":len({original_source_block_key(q) for q in test}),
          "M0_geo_climate_observer_effort_test_MSE":mse0,
          "M1_plus_source_GERS_rail_km_test_MSE":mse1,
          "M0_minus_M1_test_MSE_delta_positive_means_better":mse0-mse1,
          "M0_test_MAE":mae0,"M1_test_MAE":mae1})
    if len(folds)<MIN_INDEPENDENT_CONTINENT_FOLDS:
        return {"identified":False,
            "reason":"<4 continents have >=source sample support for leave-entire-continent-out",
            "original_source_continent_support":support,
            "n_qualified_original_continent_folds":len(folds)}
    aggregate={
       "source_continent_equal_heldout_M0_test_MSE":statistics.mean(
           x["M0_geo_climate_observer_effort_test_MSE"] for x in folds),
       "source_continent_equal_heldout_M1_with_rail_km_test_MSE":statistics.mean(
           x["M1_plus_source_GERS_rail_km_test_MSE"] for x in folds),
       "source_continent_equal_heldout_incremental_railkm_M0_minus_M1_MSE":statistics.mean(
           x["M0_minus_M1_test_MSE_delta_positive_means_better"] for x in folds),
       "source_heldout_continents_where_rail_predictor_IMPROVES_test_MSE":sum(
           x["M0_minus_M1_test_MSE_delta_positive_means_better"]>0 for x in folds),
       "source_heldout_continents_where_rail_predictor_WORSENS_test_MSE":sum(
           x["M0_minus_M1_test_MSE_delta_positive_means_better"]<0 for x in folds),
       "n_original_independent_macroregion_continent_folds":len(folds),
       "test_results_not_railway_specific_PLOT_vegetation_or_causal_train_effect":True}
    if aggregate["source_continent_equal_heldout_M0_test_MSE"]>0:
        aggregate["relative_mse_improvement_vs_baseline_pct"]=100*aggregate[
         "source_continent_equal_heldout_incremental_railkm_M0_minus_M1_MSE"]/aggregate[
         "source_continent_equal_heldout_M0_test_MSE"]
    return {"identified":True,
        "source_frozen_ridge_penalty":RIDGE_LAMBDA,
        "original_global_climate_geo_observer_control_columns":list(base),
        "new_GERS_railway_exposure_columns":list(extra),
        "source_1000km_geographic_training_and_test_block_equal_weights":True,
        "original_source_continent_support":support,
        "source_per_continent_heldout_reproducible_losses":folds,
        "source_global_continent_equal_heldout_summary":aggregate,
        "pseudo_replicated_neighbor_source_pairs_not_independent_studies":True}

def run(plant,rail,worldclim,naturalearth):
    x,qa,meta=original_sources(plant,rail,worldclim,naturalearth)
    results={}
    for width,panel in x.items():
        richness=panel["cells"]
        beta=panel["pairs"]
        results[width]={
           "original_source_equalarea_grid_width_m":int(width),
           "original_nonzero_mapped_rail_source_grid_cells":panel[
              "eligible_rail_geographic_grid_original_source_cell_n"],
           "original_botanical_source_global_2025_publicly_observed_cells_raw":panel[
              "source_botanical_input_cells_n"],
           "original_botanical_source_global_2025_public_neighbor_jaccard_pairs_raw":panel[
              "source_original_neighbor_jaccard_input_pair_n"],
           "original_botanical_source_effort_qualified_climate_valid_cells":len(richness),
           "original_botanical_source_effort_qualified_actual_neighbor_pairs":len(beta),
           "source_qualified_2025_public_plant_cells_with_nonzero_mapped_rail":sum(
              a["rail_km"]>0 for a in richness),
           "source_qualified_2025_public_plant_cells_with_zero_mapped_rail":sum(
              a["rail_km"]==0 for a in richness),
           "source_2025_original_macro_richness_predictive_continent_holdout":
              cv_country_holdout(richness,"richness"),
           "source_2025_original_macro_neighbor_Jaccard_predictive_continent_holdout":
              cv_country_holdout(beta,"beta")
        }
    return {"version":"v371-authentic-world-rail-km-2026-vs-coarse-original-iNat-plants-2025-global-macro",
      "source_two_original_global_100km_200km_study_scales":results,
      "original_source_provenance":meta,
      "source_original_genus_floristic_sources_proof_not_raw_IDs":True,
      "source_data_quality_filters_and_missingness":qa,
      "source_2025_iNat_observations_without_original_unobscured_GPS_flags":True,
      "whole_landscape_flora_only_NOT_railway-trackbed_botanical_diversity":True,
      "rail_km_is_2026_geographic_infrastructure_NOT_train_traffic":True,
      "GHSL_urbanization_and_road_density_not_yet_covariate_source":True,
      "source_WorldClim_BIO1_is_climate_covariate_NOT_species_thermal_preference":True,
      "independent_global_functional_traits_not_yet_reconciled":True,
      "source_no_world_railway_plant_adaptation_seed_disperal_causal_claim":True,
      "no_original_botanical_species_ID_observer_ID_plant_GPS_or_rail_GERS_geometry_republished":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--world-plant-source",required=True)
    p.add_argument("--world-rail-grid-source",required=True)
    p.add_argument("--worldclim-zip",required=True)
    p.add_argument("--natural-earth",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=run(a.world_plant_source,a.world_rail_grid_source,
               a.worldclim_zip,a.natural_earth)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
