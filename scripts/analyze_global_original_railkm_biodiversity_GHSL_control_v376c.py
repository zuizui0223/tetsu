#!/usr/bin/env python3
"""v376c: FULL GLOBAL GHSL 2020 built-up area adjustment for GERS rail–plants.

True original worldwide 2025 iNat plant photo observations v369,
full-source Overture rail km GIS 2026 v370 and independent JRC GHSL2020
full-world 1km built-up m2 aggregated to exact 100/200km grid v376b.

Three models with identical source photo/land-climate + genuine GHSL
availability: M0 baseline effort+latitude+BIO1, MU with settlement/built
fraction, MUR with settlement plus physical GERS railway km. Each uses
source v371 normalized 1000km-block ridge lambda 0.10 and six entire
continent heldouts, no botanical near-track or true field plot presence.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path
from statistics import mean

from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    load_json,verify_original_sources,source_climate_lookup,
    original_cell_rows,pair_rows,continent_holdout,
    RICH_BASE,RICH_RAIL,BETA_BASE,BETA_RAIL,GRID_WIDTHS,
    MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_BETA_PAIRS,
    MIN_TRAIN_BETA_PAIRS,MIN_CONTINENTS
)
from generate_global_equal_area_climate_land_grid_v351 import (
    WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
)

ORIGINAL_NE_SHA="718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b"
BUILT_DOI="10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA"
RICH_BUILT="log1p_1000_times_original_2020_builtup_fraction"
BETA_BUILT=("mean_log1p_builtup","absdiff_log1p_builtup")
ACTUAL_V371_CORRECTED={
 100000:{"richness":.004463745492835951,"jaccard":.0001476673508208617},
 200000:{"richness":.0066142438317320406,"jaccard":-.00006120015031440586}
}

def read_source_builtup(path):
    original=load_json(path)
    if not original.get("version","").startswith("v376b-publisher-original-GHSL2020"):
        raise ValueError("Not original verified source JRC GHSL 2020 global built up 1km-to-100km")
    if original.get("source_European_Commission_JRC_data_DOI")!=BUILT_DOI:
        raise ValueError("GHSL original source DOI changed")
    if len(original.get("source_JRC_GHSL2020_original_fullZIP_SHA256",""))!=64:
        raise ValueError("Original JRC full 2020 ZIP SHA256 missing")
    if original.get("global_GHSL_2020_publisher_full_ZIP_bytes")!=152527564:
        raise ValueError("Original GHSL source full ZIP extent differs")
    grids=original.get("original_global_GHSL2020_built_m2_by_100km_200km_cell_width",{})
    if set(grids)!={"100000","200000"}:
        raise ValueError("Original global 100km and 200km JRC built-up fractions missing")
    bygrid={}
    for width,sources in grids.items():
        m={}
        for z in sources:
            key=(z["source_equal_area_grid_ix"],z["source_equal_area_grid_iy"])
            if key in m:raise ValueError("Duplicate original GHSL equal-area cell")
            frac=z["source_GHSL2020_builtup_fraction_in_original_cell"]
            area=z["source_GHSL2020_total_builtup_m2_in_cell"]
            if (not isinstance(frac,(int,float)) or not math.isfinite(frac)
                 or not 0<=frac<=1.005):
                raise ValueError("Original builtup fraction invalid")
            if not math.isclose(frac*int(width)**2,area,abs_tol=1.5,rel_tol=1e-8):
                raise ValueError("Source GHSL observed builtup m2 and 100km fraction inconsistent")
            m[key]=frac
        bygrid[width]=m
    return bygrid,{
        "source_original_JRC_GHSL_2020_fullZIP_SHA256":original[
            "source_JRC_GHSL2020_original_fullZIP_SHA256"],
        "source_100km_builtup_complete_grid_cells":len(bygrid["100000"]),
        "source_200km_builtup_complete_grid_cells":len(bygrid["200000"]),
        "source_200km_is_2x2_100km_m2_rollup_verified":original[
            "aggregate_source_QA"]["total_m2_conserved_in_grid_rollup"]}

def inject_built(cells,frac):
    qualified=[]
    qa=Counter()
    for a in cells:
        value=frac.get(a["key"])
        if value is None:
            qa["source_quality_ok_plant_cell_which_lacks_original_GHSL2020_data"]+=1
            continue
        if not 0<=value<=1.005:
            raise ValueError("Impossible source settlement fraction")
        b=dict(a)
        b[RICH_BUILT]=math.log1p(1000*value)
        qualified.append(b)
        qa["qualifying_source_plant_cells_with_GHSL2020"]=1+qa[
           "qualifying_source_plant_cells_with_GHSL2020"]
        if value==0:qa["source_qualified_true_zero_builtup"]+=1
    return qualified,dict(qa)

def enrich_neighbor_pairs(original_pairs,full_cells):
    bycell={tuple(x["key"]):x for x in full_cells}
    pair_records,_=pair_rows(original_pairs,bycell,{})
    # pair_rows uses its rail predictor only from the rows; the
    # original source values remain the actual rail map, not a fake map.
    expected_keys=[]
    for source in original_pairs:
        a=tuple(source["source_cell_A_original_grid_ix_iy"])
        b=tuple(source["source_cell_B_original_grid_ix_iy"])
        if a in bycell and b in bycell and bycell[a]["continent"]==bycell[b]["continent"]:
            expected_keys.append((a,b))
    if len(expected_keys)!=len(pair_records):
        raise AssertionError("Original fixed global Jaccard source pair match lost")
    for row,(a,b) in zip(pair_records,expected_keys):
        x,y=bycell[a][RICH_BUILT],bycell[b][RICH_BUILT]
        row[BETA_BUILT[0]]=(x+y)/2
        row[BETA_BUILT[1]]=abs(x-y)
    return pair_records

def fit_same_study_continent_cv(rich,beta):
    r0=continent_holdout(rich,RICH_BASE,RICH_RAIL,"log_species",False,
             MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2,
             "original_global_recorded_logplant_species_richness")
    r1=continent_holdout(rich,RICH_BASE+(RICH_BUILT,),
             RICH_RAIL+(RICH_BUILT,),"log_species",False,
             MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2,
             "original_global_richness_add_original_GHSL_builtup")
    j0=continent_holdout(beta,BETA_BASE,BETA_RAIL,"jaccard",True,
             MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,
             "original_global_adjacent_landscape_Jaccard")
    j1=continent_holdout(beta,BETA_BASE+BETA_BUILT,
             BETA_RAIL+BETA_BUILT,"jaccard",True,
             MIN_CONTINENT_BETA_PAIRS,MIN_TRAIN_BETA_PAIRS,
             "original_global_Jaccard_add_GHSL_builtup")
    for v in (r0,r1,j0,j1):
        if v["n_original_continents_eligible_as_heldout_folds"]<MIN_CONTINENTS:
            raise ValueError("Too few source continent holdouts with 2020 independent builtup")
    return {"richness_no_GHSL_original":r0,"richness_with_GHSL":r1,
            "landscape_Jaccard_no_GHSL_original":j0,"landscape_Jaccard_with_GHSL":j1}

def model_compare(model_b,model_city):
    a=model_b["source_independent_continent_holdouts"]
    b=model_city["source_independent_continent_holdouts"]
    mapping={x["heldout_continent"]:x for x in a}
    per=[]
    for x in b:
        c=x["heldout_continent"];o=mapping[c]
        m0=o["baseline_test_blockweighted_MSE"]
        m1=o["railkm_extended_test_blockweighted_MSE"]
        mu=x["baseline_test_blockweighted_MSE"]
        mur=x["railkm_extended_test_blockweighted_MSE"]
        per.append({"heldout_continent":c,
          "observed_original_source_records":x["n_original_heldout_source_records"],
          "original_base_MSE":m0,
          "original_base_plus_rail_only_MSE":m1,
          "base_plus_JRC_2020_builtup_MSE":mu,
          "base_plus_JRC2020_and_real_GERS_rail_MSE":mur,
          "original_GERS_rail_lift_without_settlement":m0-m1,
          "JRC2020_settlement_lift_over_original_base":m0-mu,
          "GERS_rail_lift_after_JRC_settlement":mu-mur})
    n=len(per)
    means={k:mean(x[k] for x in per) for k in (
       "original_base_MSE","original_base_plus_rail_only_MSE",
       "base_plus_JRC_2020_builtup_MSE","base_plus_JRC2020_and_real_GERS_rail_MSE",
       "original_GERS_rail_lift_without_settlement",
       "JRC2020_settlement_lift_over_original_base",
       "GERS_rail_lift_after_JRC_settlement")}
    means["relative_railkm_prediction_lift_after_settlement_pct"]=(
       100*means["GERS_rail_lift_after_JRC_settlement"]/
       means["base_plus_JRC_2020_builtup_MSE"])
    return {"n_entire_continents_heldout":n,
            "mean_six_continent_equal_source_MSE":means,
            "n_continents_where_railkm_still_helps_after_builtup":
               sum(x["GERS_rail_lift_after_JRC_settlement"]>0 for x in per),
            "n_continents_where_urbanization_helps_original_base":
               sum(x["JRC2020_settlement_lift_over_original_base"]>0 for x in per),
            "source_heldout_continent_details":per,
            "no_causal_rail_or_plant_source_claim":True}

def analyze(plants_path,rail_path,ghsl_path,worldclim_path,ne_path):
    rawplant=load_json(plants_path)
    rawrail=load_json(rail_path)
    checks=verify_original_sources(rawplant,rawrail)
    ghsl,urbanq=read_source_builtup(ghsl_path)
    if sha_file(worldclim_path)!=WORLDCLIM_SHA or sha_file(ne_path)!=ORIGINAL_NE_SHA:
        raise ValueError("World original climate or land source hashes differ")
    land,landqa=source_continents(ne_path)
    raster,mem,rasterqa=worldclim_tif(worldclim_path)
    panels=[]
    try:
        scale,scaleqa=climate_scale(raster)
        for width in GRID_WIDTHS:
            source=rawplant["source_grids"][str(width)]
            plantgrids=source[
                "source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
            neighbors=source[
                "source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs"]
            sourcekeys={(x["source_equal_area_grid_ix"],x["source_equal_area_grid_iy"])
                        for x in plantgrids}
            climate,climqa=source_climate_lookup(land,raster,scale,sourcekeys,width)
            rails=rawrail[
              "source_global_gridded_rail_exposure_km_by_cell_width"][str(width)]
            valid,qa,rail_bycell=original_cell_rows(plantgrids,rails,climate)
            # Baseline and settlement models have EXACT same botanical
            # source units after GHSL availability gate.
            subset,builtqa=inject_built(valid,ghsl[str(width)])
            b=pair_rows(neighbors,{x["key"]:x for x in subset},rail_bycell)[0]
            # No difference in acceptable pair source ordering/metadata
            # when built-up source is added.
            enriched=enrich_neighbor_pairs(neighbors,subset)
            if len(b)!=len(enriched):
                raise AssertionError("GHSL pair method changes botanical Jaccard source")
            scored=fit_same_study_continent_cv(subset,enriched)
            comparison_rich=model_compare(
                scored["richness_no_GHSL_original"],scored["richness_with_GHSL"])
            comparison_beta=model_compare(
                scored["landscape_Jaccard_no_GHSL_original"],scored["landscape_Jaccard_with_GHSL"])
            if len(subset)==len(valid):
                corrected=ACTUAL_V371_CORRECTED[width]
                for typeout,observed in (
                    (scored["richness_no_GHSL_original"],corrected["richness"]),
                    (scored["landscape_Jaccard_no_GHSL_original"],corrected["jaccard"])):
                    if abs(typeout["mean_holdout_continent_delta_MSE"]-observed)>1e-9:
                        raise AssertionError("Original v371 normalized source result not reproducible in GHSL analysis")
            panels.append({"source_grid_width_m":width,
              "original_botanical_cells_which_qualified_without_builtup":len(valid),
              "source_GHSL2020_join_QA":builtqa,
              "n_botanical_source_cells_with_full_covariates":len(subset),
              "n_original_adjacent_source_landscape_Jaccard_pairs_with_GHSL":len(enriched),
              "source_richness_global_2020_settlement_and_true_railkm":comparison_rich,
              "source_Jaccard_global_2020_settlement_and_true_railkm":comparison_beta})
    finally:
        raster.close();mem.close()
    return {"version":"v376c-GHSL-2020-world-builtup-adjusted-true-global-railkm-photo-plant-biodiversity-macro",
        "source_GHSL_2020_EC_JRC_DOI":BUILT_DOI,
        "source_GHSL_2020_orig_fingerprint":urbanq,
        "source_original_global_2025_photo_plant_SHA256":checks.get("original_iNat_obs_gzip_SHA256"),
        "source_original_fullworld_2026_conventional_railkm":2237043.481,
        "original_NaturalEarth_source_QA":landqa,
        "source_WorldClim_BIO1_QA":rasterqa,
        "original_continents_count":6,
        "source_macro_by_width":panels,
        "scientific_limits":[
           "GHSL2020 built surface fraction is a partial independent urbanization proxy, NOT observed roads, population, local train service or actual railway ballast plant habitat",
           "2025 iNaturalist public photo-linked plant records are 100km/200km whole-landscape detections, not standardized trackbed flora; original open metadata have no precise unobscured-privacy audit",
           "All worldwide source 100km/200km records share same global observations and are not independent field ecological replicates",
           "All 6 entire continent heldouts/1000km block-weighted ridge cannot identify causality or full local human-development/observability confounding",
           "Urban built-up 2020 predates rail2026 and plant2025 spatial layers, physical geo statistics have temporal mismatch",
           "WorldClim original BIO1 is site/landscape climate not independent plant species thermal affinity",
           "No source rail train-mediated dispersal, seed establishment, species adaptation or functional convergence concluded",
           "No individual original GPS, species/observer IDs, private source railway linework or original 1km GHSL raster pixels published"]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-v369-global-plants",required=True)
    p.add_argument("--original-v370-global-rail",required=True)
    p.add_argument("--original-JRC-GHSL2020-global-builtup",required=True)
    p.add_argument("--natural-earth",required=True)
    p.add_argument("--worldclim-zip",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=analyze(a.original_v369_global_plants,a.original_v370_global_rail,
        a.original_JRC_GHSL2020_global_builtup,a.worldclim_zip,a.natural_earth)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf8")
    print(json.dumps({"version":out["version"],
        "source_macro_by_width":out["source_macro_by_width"],
        "source_GHSL_2020_orig_fingerprint":out["source_GHSL_2020_orig_fingerprint"]},
        ensure_ascii=False),flush=True)
if __name__=="__main__":main()
