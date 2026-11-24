#!/usr/bin/env python3
"""v371: TRUE worldwide regional plant observations x original rail-km macro model.

Direct joins of source-authenticated FULL-WORLD 2025 iNaturalist coarse
observation and botanical Jaccard data (v369) to source-authenticated
2.237 million physical railway km 100/200km equal-area exposures (v370).

Does NOT measure railway-specific botanical field plots, real ballast flora,
operated train service, measured railway microclimate, individual seed
carriage or plant species thermal preference; BIO1 is regional climate only.

Frozen before full 2025 source result: see docs/tetsu_v371_true_global_
rail_plant_turnover_macro_source_frozen_ja.md.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

from generate_global_equal_area_climate_land_grid_v351 import (
    source_continents,worldclim_tif,read_temp,climate_scale,sha_file,
    WORLDCLIM_SHA,NATURAL_EARTH_COMMIT,CONTINENTS
)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    X_START,Y_START,OBS_ETAG,TAXA_ETAG,MIN_PAIR_CELL_OBS,
    MIN_PAIR_CELL_OBSERVERS,MIN_PAIR_CELL_PLANT_SPECIES
)

RAIL_SOURCE_KM=2237043.481
RAIL_RELEASE="2026-11-18.0"
SOURCE_YEAR=2025
GRID_WIDTHS=(100_000,200_000)
MIN_CONTINENT_RICH_CELLS=100
MIN_CONTINENT_BETA_PAIRS=50
MIN_TRAIN_BETA_PAIRS=300
MIN_CONTINENTS=4
RIDGE=0.1
SOURCE_GEOMETRY_COORDINATE_EPSG="EPSG:6933"
BLOCK_WIDTH_M=1_000_000

RICH_BASE=("log_observations","log_observers","bio1","bio1_sq","abs_lat","abs_lat_sq")
RICH_RAIL=RICH_BASE+("log_rail_km",)
BETA_BASE=("mean_log_observations","absdiff_log_observations",
           "mean_log_observers","absdiff_log_observers",
           "mean_log_species","absdiff_log_species",
           "mean_bio1","absdiff_bio1","mean_abs_lat","orientation_horizontal")
BETA_RAIL=BETA_BASE+("log_mean_rail_km","absdiff_log_rail_km")

def verify_original_sources(plants,rails):
    if plants.get("source_summary",{}).get("original_source_iNat_obs_S3_ETag")!=OBS_ETAG:
        raise ValueError("Not frozen official full global original iNaturalist source 2026-10-27")
    if plants["source_summary"].get("original_source_iNat_taxa_S3_ETag")!=TAXA_ETAG:
        raise ValueError("Original global iNaturalist taxonomy snapshot changed")
    qa=plants["source_summary"].get("source_global_photo_observations_QA",{})
    if qa.get("original_global_iNat_photo_linked_rows_read",0)<50_000_000:
        raise ValueError("Original 13.7GB worldwide observations stream not complete")
    if set(plants.get("source_grids",{}))!={"100000","200000"}:
        raise ValueError("Original global plant grids incomplete")
    if not rails.get("source_original_fullworld_rail_128_files_SHA_verified"):
        raise ValueError("Not full-world all 128 Overture rail source files")
    if rails.get("source_Oveture_release")!=RAIL_RELEASE:
        raise ValueError("Original worldwide railway release differs from frozen source")
    if abs(rails.get("source_fullworld_candidate_normal_rail_km",0)-RAIL_SOURCE_KM)>.06:
        raise ValueError("Not source genuine 2.237 million km active-candidate railway population")
    if set(rails.get("source_global_gridded_rail_exposure_km_by_cell_width",{}))!={"100000","200000"}:
        raise ValueError("Original 100km and 200km rail exposure census missing")
    sums=[]
    for width in ("100000","200000"):
        entries=rails["source_global_gridded_rail_exposure_km_by_cell_width"][width]
        val=sum(float(x["source_real_original_Overture_normal_rail_km_clipped_into_cell"])
                for x in entries)
        if abs(val-RAIL_SOURCE_KM)>0.15:
            raise ValueError("Original global gridded railway km do not sum to true source world total")
        sums.append(val)
    if abs(sums[0]-sums[1])>.05:
        raise ValueError("Source railway 100km-200km grid totals disagree")
    plant_counts=[]
    for width in ("100000","200000"):
        x=plants["source_grids"][width]["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
        n=sum(int(t["source_year_2025_photo_linked_plants_observation_count"]) for t in x)
        plant_counts.append(n)
    if plant_counts[0]!=plant_counts[1]:
        raise ValueError("Original same-world iNaturalist 2025 100km/200km counts do not match")
    return {"original_2025_full_world_plant_detections_in_grid":plant_counts[0],
        "source_2026_normal_rail_100km_200km_length_totals_km":sums,
        "original_iNat_obs_gzip_SHA256":plants["source_summary"].get(
            "source_iNat_original_2026_10_27_observations_gzip_SHA256_full"),
        "original_iNat_taxa_gzip_SHA256":plants["source_summary"].get(
            "source_iNat_original_2026_10_27_taxa_gzip_SHA256_full")}

def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf8"))

def valid_source_count(obs,observers,species):
    return (obs>=MIN_PAIR_CELL_OBS and observers>=MIN_PAIR_CELL_OBSERVERS
            and species>=MIN_PAIR_CELL_PLANT_SPECIES)

def cell_center(ix,iy,width):
    if not 0<=ix<=(34_800_000//width) or not 0<=iy<=(14_800_000//width):
        raise ValueError("Source global EPSG6933 cell index outside source frame")
    return (X_START+(ix+0.5)*width,Y_START+(iy+0.5)*width)

def block_id(ix,iy,width):
    x,y=cell_center(ix,iy,width)
    return (math.floor(x/BLOCK_WIDTH_M),math.floor(y/BLOCK_WIDTH_M))

def source_climate_lookup(land,raster,scale,points,width):
    from shapely.geometry import Point
    from pyproj import Transformer
    t=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    out={};qa=Counter()
    for ix,iy in points:
        x,y=cell_center(ix,iy,width)
        lon,lat=t.transform(x,y)
        if not all(math.isfinite(k) for k in (lat,lon)) or not(-90<=lat<=90):
            qa["invalid_proj"]+=1
            continue
        pp=Point(lon,lat)
        found=[cont for cont in CONTINENTS if land[cont].covers(pp)]
        if len(found)!=1:
            qa["cell_center_on_ocean_or_ambiguous"]+=1
            continue
        bio=read_temp(raster,lat,lon)
        if bio is None:
            qa["original_worldclim_temperature_missing"]+=1
            continue
        tc=bio/scale
        if not -70<=tc<=50:
            raise ValueError("Original WorldClim BIO1 physical unit implausible")
        out[ix,iy]={"continent":found[0],"bio1":tc,
                     "abs_lat":abs(lat),"block":block_id(ix,iy,width)}
        qa["source_land_center_and_climate_valid"]+=1
    return out,dict(qa)

def original_cell_rows(original_plant_grid,original_rail_grid,clim):
    rail={}
    for z in original_rail_grid:
        key=z["source_equal_area_grid_ix"],z["source_equal_area_grid_iy"]
        if key in rail:raise ValueError("Original global rail grid cell repeated")
        km=float(z["source_real_original_Overture_normal_rail_km_clipped_into_cell"])
        if not math.isfinite(km) or km<=0:
            raise ValueError("Original positive rail census has invalid km")
        rail[key]=km
    count=Counter()
    rows=[]
    ids=set()
    for source in original_plant_grid:
        key=(source["source_equal_area_grid_ix"],source["source_equal_area_grid_iy"])
        if key in ids:raise ValueError("Original global 2025 plant cell repeated")
        ids.add(key)
        count["original_source_2025_plant_occupied_grid_cells"]+=1
        o=int(source["source_year_2025_photo_linked_plants_observation_count"])
        b=int(source["source_2025_distinct_observers_with_photo_linked_plant_records"])
        s=int(source["source_2025_plant_species_richness"])
        if any(x<0 for x in (o,b,s)) or s>o or b>o:
            raise ValueError("Original plant source richness/effort inconsistent")
        if not valid_source_count(o,b,s):
            count["source_low_observer_or_detection_effort_cells"]+=1
            continue
        geo=clim.get(key)
        if geo is None:
            count["source_plant_cell_no_land_center_or_worldclim"]+=1
            continue
        rk=rail.get(key,0.)
        log_o,log_b,log_s=map(math.log1p,(o,b,s))
        row={"key":key,"continent":geo["continent"],"block":geo["block"],
             "observations":o,"observers":b,"species":s,"rail_km":rk,
             "log_species":log_s,
             "log_observations":log_o,"log_observers":log_b,
             "bio1":geo["bio1"],"bio1_sq":geo["bio1"]**2,
             "abs_lat":geo["abs_lat"],"abs_lat_sq":geo["abs_lat"]**2,
             "log_rail_km":math.log1p(rk)}
        rows.append(row)
        count["eligible_original_worldland_plant_cells"]+=1
        if rk>0:count["eligible_cells_with_original_rail_km"]+=1
    if len(rows)!=len({v["key"] for v in rows}):
        raise AssertionError("Duplicate eligible plant source cells")
    return rows,dict(count),rail

def pair_rows(pair_source,bycell,rail):
    rows=[];qa=Counter()
    for p in pair_source:
        qa["original_world_landscape_neighbor_pairs"]+=1
        a=tuple(p["source_cell_A_original_grid_ix_iy"])
        b=tuple(p["source_cell_B_original_grid_ix_iy"])
        if len(a)!=2 or len(b)!=2:
            raise ValueError("Invalid original two-dimensional source grid key")
        dx,dy=b[0]-a[0],b[1]-a[1]
        if (dx,dy) not in ((0,1),(1,0)):
            raise ValueError("Original source beta pair not a unique horizontally/vertically adjacent land-grid pair")
        x=bycell.get(a);y=bycell.get(b)
        if not x or not y:
            qa["neighbor_endpoint_fails_plant_or_continent_or_climate_gate"]+=1
            continue
        if x["continent"]!=y["continent"]:
            qa["source_pair_crosses_continent_boundary"]+=1
            continue
        j=float(p["source_biological_Jaccard_similarity_between_ENTIRE_landscape_cells"])
        if not 0<=j<=1:
            raise ValueError("Original source botanical beta Jaccard outside [0,1]")
        if int(p["source_grid_A_original_plant_species_n"])!=x["species"] or int(
           p["source_grid_B_original_plant_species_n"])!=y["species"]:
            raise ValueError("Original source Jaccard record flora richness differs from original cell census")
        avg=lambda key:(x[key]+y[key])/2
        diff=lambda key:abs(x[key]-y[key])
        ra=math.log1p(x["rail_km"]);rb=math.log1p(y["rail_km"])
        row={"continent":x["continent"],
             "block":(math.floor((x["block"][0]+y["block"][0])/2),
                      math.floor((x["block"][1]+y["block"][1])/2)),
             "jaccard":j,
             "mean_log_observations":avg("log_observations"),
             "absdiff_log_observations":diff("log_observations"),
             "mean_log_observers":avg("log_observers"),
             "absdiff_log_observers":diff("log_observers"),
             "mean_log_species":avg("log_species"),
             "absdiff_log_species":diff("log_species"),
             "mean_bio1":avg("bio1"),
             "absdiff_bio1":diff("bio1"),
             "mean_abs_lat":avg("abs_lat"),
             "orientation_horizontal":1. if dx else 0.,
             "log_mean_rail_km":math.log1p((x["rail_km"]+y["rail_km"])/2.),
             "absdiff_log_rail_km":abs(ra-rb)}
        rows.append(row)
        qa["valid_landscape_beta_pairs"]+=1
    return rows,dict(qa)

def block_weights(rows):
    count=Counter((r["continent"],r["block"]) for r in rows)
    w=np.array([1/count[(r["continent"],r["block"])] for r in rows],dtype=float)
    if not np.all(np.isfinite(w)) or np.any(w<=0):
        raise AssertionError("Nonpositive original block weights")
    return w

def fit_predict(train,test,predictors,target,clip_predictions):
    if len(train)<len(predictors)*4:
        raise ValueError("Insufficient source training support for frozen ridge model")
    weights=block_weights(train)
    X=np.array([[float(r[k]) for k in predictors] for r in train],dtype=float)
    Z=np.array([[float(r[k]) for k in predictors] for r in test],dtype=float)
    y=np.array([float(r[target]) for r in train],dtype=float)
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(Z)) and np.all(np.isfinite(y))):
        raise ValueError("Original global predictor has nonfinite values")
    mu=np.average(X,axis=0,weights=weights)
    var=np.average((X-mu)**2,axis=0,weights=weights)
    sd=np.sqrt(var);sd[sd<1e-12]=1.
    A=np.column_stack((np.ones(len(train)),(X-mu)/sd))
    B=np.column_stack((np.ones(len(test)),(Z-mu)/sd))
    # Each source 1000km block has total weight exactly one.
    weighted=A*weights[:,None]
    gram=A.T@weighted
    penalty=np.eye(gram.shape[0])*RIDGE
    penalty[0,0]=0.
    beta=np.linalg.solve(gram+penalty,A.T@(weights*y))
    pred=B@beta
    return np.clip(pred,0,1) if clip_predictions else pred

def metrics(rows,yp,actual_key):
    y=np.array([float(r[actual_key]) for r in rows],dtype=float)
    w=block_weights(rows)
    residual=y-yp
    mse=np.average(residual**2,weights=w)
    mae=np.average(np.abs(residual),weights=w)
    return float(mse),float(mae)

def continent_holdout(rows,base,extended,target,clip_predictions,min_holdout,
                       min_training,metric_name):
    bycontinent=Counter(r["continent"] for r in rows)
    preds=[]
    for continent in CONTINENTS:
        test=[r for r in rows if r["continent"]==continent]
        train=[r for r in rows if r["continent"]!=continent]
        if len(test)<min_holdout or len(train)<min_training:
            continue
        p0=fit_predict(train,test,base,target,clip_predictions)
        p1=fit_predict(train,test,extended,target,clip_predictions)
        mse0,mae0=metrics(test,p0,target)
        mse1,mae1=metrics(test,p1,target)
        preds.append({"heldout_continent":continent,
            "n_original_heldout_source_records":len(test),
            "n_training_source_records":len(train),
            "n_independent_heldout_1000km_grid_blocks":len(set(r["block"] for r in test)),
            "n_independent_training_1000km_blocks":len(set((r["continent"],r["block"]) for r in train)),
            "baseline_test_blockweighted_MSE":mse0,
            "railkm_extended_test_blockweighted_MSE":mse1,
            "delta_MSE_positive_means_rail_added_prediction_helped":mse0-mse1,
            "baseline_test_blockweighted_MAE":mae0,
            "railkm_extended_test_blockweighted_MAE":mae1,
            "delta_MAE_positive_means_rail_added_prediction_helped":mae0-mae1})
    return {"outcome":metric_name,
      "training_frozen_baseline_predictors":list(base),
      "rail_extended_predictors":list(extended),
      "original_source_continent_counts_before_holdout":dict(bycontinent),
      "minimum_heldout_original_records_per_continent":min_holdout,
      "minimum_training_original_records":min_training,
      "ridge_penalty_lambda":RIDGE,
      "n_original_continents_eligible_as_heldout_folds":len(preds),
      "frozen_minimum_continental_holdout_folds_for_interpretation":MIN_CONTINENTS,
      "independent_cross_continent_test_supported":len(preds)>=MIN_CONTINENTS,
      "n_holdout_continents_with_positive_delta_MSE":sum(x[
          "delta_MSE_positive_means_rail_added_prediction_helped"]>0 for x in preds),
      "mean_holdout_continent_delta_MSE":(sum(x[
          "delta_MSE_positive_means_rail_added_prediction_helped"] for x in preds)/len(preds) if preds else None),
      "mean_holdout_continent_delta_MAE":(sum(x[
          "delta_MAE_positive_means_rail_added_prediction_helped"] for x in preds)/len(preds) if preds else None),
      "source_independent_continent_holdouts":preds,
      "all_source_model_comparisons_are_REGIONAL_landscape_species_observations_NOT_railway_track_flora":True}

def analyze(plantpath,railpath,geojson,worldclim):
    plant=load_json(plantpath);rail=load_json(railpath)
    source=verify_original_sources(plant,rail)
    if sha_file(geojson)!="718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b":
        raise ValueError("Original pinned Natural Earth global land source changed")
    if sha_file(worldclim)!=WORLDCLIM_SHA:
        raise ValueError("Original pinned WorldClim BIO ZIP changed")
    land,landqa=source_continents(geojson)
    raster,mem,rasqa=worldclim_tif(worldclim)
    try:
        scale,climqa=climate_scale(raster)
        results=[]
        for width in GRID_WIDTHS:
            k=str(width)
            plants=plant["source_grids"][k]
            plantcell=plants["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            sourcekeys={(p["source_equal_area_grid_ix"],p["source_equal_area_grid_iy"])
                         for p in plantcell}
            clim,geoqa=source_climate_lookup(land,raster,scale,sourcekeys,width)
            rails=rail["source_global_gridded_rail_exposure_km_by_cell_width"][k]
            eligible,cqa,rail_bykey=original_cell_rows(plantcell,rails,clim)
            bycell={r["key"]:r for r in eligible}
            beta,bqa=pair_rows(
                plants["source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"],
                bycell,rail_bykey)
            rich=continent_holdout(eligible,RICH_BASE,RICH_RAIL,"log_species",False,
                                   MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2,
                                   "log1p_global_2025_PHOTO_RECORDED_plant_species_richness")
            turn=continent_holdout(beta,BETA_BASE,BETA_RAIL,"jaccard",True,
                                   MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,
                                   "adjacent_landscape_plant_species_Jaccard_similarity")
            results.append({"original_source_grid_m":width,
               "original_source_land_climate_cell_QA":geoqa,
               "original_plant_source_effort_QA":cqa,
               "original_source_adjoining_landscape_flora_QA":bqa,
               "original_full_world_railkm_2026_positive_cells_n":len(rail_bykey),
               "original_quality_support_eligible_2025_photo_plant_cells":len(eligible),
               "original_quality_support_valid_adjacent_2025_photo_plant_pairs":len(beta),
               "richness_model":rich,"landscape_botanical_Jaccard_model":turn})
    finally:
        raster.close();mem.close()
    return {"version":"v371-original-fullworld-railkm-landscape-plant-richness-Jaccard-cross-continent-holdout",
      "source_original_rail_exposure_DOI":"Overture Maps 2026-11-18.0 original 128/128 Parquet",
      "source_original_global_iNaturalist":"open-data observations.csv.gz 2026-10-27 pinned S3 ETag, 2025 original 2025 active research-grade wild photo-linked Plantae species",
      "source_original_climate":"WorldClim 2.1 1970-2000 BIO1 physical annual temperature",
      "source_original_land":"Natural Earth 110m original pinned map-units 2026 source",
      "source_integrity":source,"original_source_land_continent_QA":landqa,
      "original_source_WorldClim_QA":rasqa,
      "original_climate_scale_QA":climqa,
      "grids":results,
      "conclusions_precommitment":[
        "Only regional whole-landscape original recorded botanical richness and adjacent-cell Jaccard, NEVER actual near-rail ballast botanical communities",
        "The 2025 iNaturalist public original location data lack source coordinate obscuration flags; this analysis does NOT use proximity to individual plants",
        "Source 2026 physical railway centreline km may include unoperated lines and does not measure train traffic or seed transport",
        "Neither positive nor negative continent heldout predictive performance identifies causal railway vegetation mechanisms",
        "Research-grade public photo detection effort is geographically nonrandom and richness is not species completeness or botanical abundance",
        "Original 1970-2000 WorldClim site BIO1 is not species-specific temperature affinity or current rail microclimate",
        "Continent holdouts 100km vs 200km analyses share same global underlying observations and do NOT constitute independent biological studies",
        "Source 1000km source spatial block weighting reduces but does not remove spatial phylogenetic and sample-effort confounding",
        "No global seed dispersal, plant functional convergence, climate-range expansion, or within-train adaptation inferred without new independent source data"]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-v369-global-plants",required=True)
    p.add_argument("--original-v370-global-rail",required=True)
    p.add_argument("--natural-earth",required=True)
    p.add_argument("--worldclim-zip",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    d=analyze(a.original_v369_global_plants,a.original_v370_global_rail,
              a.natural_earth,a.worldclim_zip)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({"version":d["version"],"source_integrity":d["source_integrity"],
        "grids":[{"source_grid_m":x["original_source_grid_m"],
             "n_eligible_source_plant_cells":x["original_quality_support_eligible_2025_photo_plant_cells"],
             "n_eligible_source_pairs":x["original_quality_support_valid_adjacent_2025_photo_plant_pairs"],
             "richness_m0_m1":{k:x["richness_model"][k] for k in (
                 "n_original_continents_eligible_as_heldout_folds",
                 "independent_cross_continent_test_supported",
                 "n_holdout_continents_with_positive_delta_MSE",
                 "mean_holdout_continent_delta_MSE",
                 "source_independent_continent_holdouts")},
             "Jaccard_m0_m1":{k:x["landscape_botanical_Jaccard_model"][k] for k in (
                 "n_original_continents_eligible_as_heldout_folds",
                 "independent_cross_continent_test_supported",
                 "n_holdout_continents_with_positive_delta_MSE",
                 "mean_holdout_continent_delta_MSE",
                 "source_independent_continent_holdouts")}}
             for x in d["grids"]]},ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
