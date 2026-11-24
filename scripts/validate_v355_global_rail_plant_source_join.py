#!/usr/bin/env python3
"""v355: Pure standard-library audit of source-authenticated global 32-cell join.

All original data already sourced, validated and source SHA recorded:
- v351 official Natural Earth×WorldClim source success run 45330702677
- v354 original six-continent fixed-release Overture source success 45331600841
- v353 original official iNaturalist 2025 source complete 45331986000.
No network calls or plant species/GPS/user data. Check no 32-cell omissions,
original sample hashes, original stratum N and exact inclusion probabilities,
and independent published source totals. Estimate rail presence only for the
defined finite GLOBAL LAND GRID CENTER population, not rail-network length.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict
from pathlib import Path

EXPECTED_SAMPLE_SHA="b65a13e32d65e8fa818bbefe2d7fd1f0e15c6c5cb97ff6f6212aea23cf138bf2"
EXPECTED_WORLDCLIM_SHA="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
EXPECTED_NE_SHA="718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b"
EXPECTED_SOURCE_RUNS=(45331600841,45331986000)
EXPECTED_N=13437
EXPECTED_SAMPLE_N=32
EXPECTED_RAIL_PRESENT=13
EXPECTED_RAIL_KM=1245.517
EXPECTED_INAT_RECTANGLE_SUM=43025
EXPECTED_ZERO_INAT_RECTANGLES=11
EXPECTED_GE100_INAT_RECTANGLES=15
SUM_TOL=5e-7

def validate(frame,result):
    if frame.get("sampled_cell_unique_hash")!=EXPECTED_SAMPLE_SHA:
        raise ValueError("Underlying source-drawn global probability sample SHA drift")
    qa=frame.get("source_input_QA",{})
    if qa.get("WorldClim_BIO1",{}).get("source_zip_sha256")!=EXPECTED_WORLDCLIM_SHA:
        raise ValueError("Original world climatology source checksum drift")
    if qa.get("Natural_Earth",{}).get("source_file_sha256")!=EXPECTED_NE_SHA:
        raise ValueError("Original Natural Earth source checksum drift")
    N=frame["source_full_cell_grid_QA"]["source_climate_supported_land_cells"]
    if N!=EXPECTED_N:
        raise ValueError("Original finite land cell population drift")
    prov=result.get("provenance",{})
    if prov.get("world_cell_sample_sha256")!=EXPECTED_SAMPLE_SHA:
        raise ValueError("v355 read wrong source sampling frame")
    for name,runid in (("rail_32_original_source_actions_url",45331600841),
                       ("public_plant_32_original_source_actions_url",45331986000)):
        if not prov.get(name,"").endswith("/"+str(runid)):
            raise ValueError("Original source Actions provenance changed: "+name)
    actual=result.get("source_cells_aggregate_only",[])
    if len(actual)!=EXPECTED_SAMPLE_N:
        raise ValueError("Source 32 independent original cells missing")
    by_id={c["grid_cell_ix_iy"]:c for c in actual}
    if len(by_id)!=EXPECTED_SAMPLE_N:
        raise ValueError("Source original grid cell duplicated")

    strata={}
    for s in frame["strata"]:
        key=(s["continent"],s["source_macrotemperature_band_C"])
        cells=s["sampled_cells"]
        pi=s["source_finite_population_inclusion_probability_per_cell"]
        n=s["n_equal_area_source_land_climate_cells"]
        strata[key]=(n,pi,cells)
        for x in cells:
            id=x["source_original_fixed_grid_index"]
            if id not in by_id:
                raise ValueError("Original globally sampled cell absent: "+id)
            y=by_id[id]
            if y["source_cell_sha256"]!=x["source_grid_cell_sha256"]:
                raise ValueError("Original source cell hash drift")
            if y["continent"]!=s["continent"] or y["climate_band"]!=key[1]:
                raise ValueError("Original continent/band mislabelled")
            if y["source_landcell_stratum_N"]!=n:
                raise ValueError("Finite original stratum population changed")
            if y["source_design_inclusion_probability"]!=pi:
                raise ValueError("Original sampling inclusion pi changed")
            if abs(y["source_BIO1_C"]-x["original_macroclimate_BIO1_C"])>0.0005001:
                raise ValueError("Original macroclimate WorldClim temperature drift")
            rail=y["original_mapped_rail_km_in_exact_100km_cell"]
            plants=y["original_iNaturalist_2025_PLANTAE_public_annual_record_count_in_RECTANGULAR_bbox"]
            if (not isinstance(rail,(int,float)) or not math.isfinite(rail) or rail<0
                    or bool(y["source_mapped_rail_present"])!=(rail>0)):
                raise ValueError("Original mapped rail source invalid")
            if not isinstance(plants,int) or isinstance(plants,bool) or plants<0:
                raise ValueError("Original iNat source public count invalid")
    if len(strata)!=18 or sum(n for n,p,cells in strata.values())!=N:
        raise ValueError("Original 18 source climate×continent strata or population drift")
    if set(by_id)!={x["source_original_fixed_grid_index"] for n,p,cells in strata.values()
                   for x in cells}:
        raise ValueError("v355 source rows outside predeclared 32-cell frame")
    good=list(by_id.values())
    nrail=sum(x["source_mapped_rail_present"] for x in good)
    length=math.fsum(x["original_mapped_rail_km_in_exact_100km_cell"] for x in good)
    inat=sum(x["original_iNaturalist_2025_PLANTAE_public_annual_record_count_in_RECTANGULAR_bbox"] for x in good)
    zeros=sum(x["original_iNaturalist_2025_PLANTAE_public_annual_record_count_in_RECTANGULAR_bbox"]==0 for x in good)
    ge100=sum(x["original_iNaturalist_2025_PLANTAE_public_annual_record_count_in_RECTANGULAR_bbox"]>=100 for x in good)
    if (nrail!=EXPECTED_RAIL_PRESENT or abs(length-EXPECTED_RAIL_KM)>SUM_TOL or
            inat!=EXPECTED_INAT_RECTANGLE_SUM or zeros!=EXPECTED_ZERO_INAT_RECTANGLES or
            ge100!=EXPECTED_GE100_INAT_RECTANGLES):
        raise ValueError("Genuine source v353 or v354 independent full-run summary mismatched")
    value=0.0
    for (continent,band),(Nh,pi,fields) in strata.items():
        if not fields:continue
        yh=sum(by_id[x["source_original_fixed_grid_index"]]["source_mapped_rail_present"]
               for x in fields)/len(fields)
        value+=Nh/N*yh
    saved=result["global_sample"]
    if not math.isclose(value,saved["stratified_SRSWOR_design_weighted_estimate_mapped_rail_presence_per_SOURCE_LAND_GRID_CELL"],abs_tol=1e-12):
        raise ValueError("Calculated global finite-grid HT mapped rail descriptor mismatch")
    if not saved["source_citizen_counts_are_rectangular_BBOX_not_exact_100km_landcell_or_near_track"]:
        raise ValueError("Cannot relabel source plant query bounding boxes as actual rail-near flora")
    if not saved["no_global_railway_plant_warmth_or_train_adaptation_inferred"]:
        raise ValueError("Source presence/observation opportunity cannot claim biological thermal effect")
    return {"version":"v355-original-source-32-global-grid-join-revalidated",
      "original_source_sample_cells":len(good),
      "original_source_land_grid_population_cells":N,
      "mapped_rail_present_source_cells":nrail,
      "original_sampled_source_clipped_rail_km":length,
      "original_2025_iNat_rectangular_bbox_total_source_observations_not_unique":inat,
      "source_2025_no_plant_observation_bboxes":zeros,
      "source_2025_ge100_plant_observation_bboxes":ge100,
      "source_rail_presence_finite_land_GRID_weighted_sample_descriptor":value,
      "no_source_species_temperature_trait_or_rail_biological_effect":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--frame",required=True)
    p.add_argument("--joined",required=True)
    args=p.parse_args()
    f=json.loads(Path(args.frame).read_text(encoding="utf8"))
    j=json.loads(Path(args.joined).read_text(encoding="utf8"))
    print(json.dumps(validate(f,j),ensure_ascii=False,indent=2))

if __name__=="__main__":main()
