#!/usr/bin/env python3
"""v354: All 32 originally selected world land-grid cells vs true Overture rail.

No cherry-picking of rail-present sampled locations: every v351 world source
candidate checked or labelled source unavailable. Original 2026-11-18 Overture
fixed release and railway filter identical to independently successful v352.
"""
from __future__ import annotations
import argparse,json,statistics,time
from collections import Counter,defaultdict
from pathlib import Path

from audit_global_landcell_Overture_rail_v352 import (
    source_overture,RAIL_RELEASE,SOURCE_FRAME_SHA
)
from audit_original_iNat_world_probability_cells_v353 import selected_cells

CONT_CODES={"af":"Africa","as":"Asia","eu":"Europe",
            "na":"North America","sa":"South America","oc":"Oceania"}
WAVE1_REFERENCE=(
  "data/derived/v352_original_global_temperate_landcell_Overture_rail_first6_aggregate.json"
)

def source_partition(frame,continent):
    if continent not in CONT_CODES.values():
        raise ValueError("Unknown frozen 6-continent source stratum")
    cohort=[r for r in selected_cells(frame) if r["continent"]==continent]
    if len(cohort) not in (4,6):
        raise ValueError("Unexpected frozen world source cell count for continent")
    return sorted(cohort,key=lambda r:r["cell"]["source_original_fixed_grid_index"])

def analyze_partition(frame_path,continent,output):
    frame=json.loads(Path(frame_path).read_text(encoding="utf8"))
    cells=source_partition(frame,continent)
    out={"version":"v354-true-Overture-2026-global-land-grid-original-source-partition",
         "source_frame_SHA":SOURCE_FRAME_SHA,"source_Overture_release":RAIL_RELEASE,
         "source_continent":continent,
         "n_original_cells_in_continent":len(cells),
         "srswor_original_landcell_sample_not_global_rail_length_weighted":True,
         "source_real_rail_and_plant_2025_original_geoposition_NOT_reexported":True,
         "source_original_cells":[]}
    target=Path(output);target.parent.mkdir(parents=True,exist_ok=True)
    for i,cell in enumerate(cells):
        if i:time.sleep(2.0)
        source=source_overture(cell["cell"])
        row={"continent":continent,
             "source_climate_band":cell["macrotemp_stratum"],
             "source_original_grid_ix_iy":cell["cell"]["source_original_fixed_grid_index"],
             "source_cell_sha256":cell["cell"]["source_grid_cell_sha256"],
             "source_original_BIO1_C":cell["cell"]["original_macroclimate_BIO1_C"],
             "source_population_landcell_N":cell["population_n"],
             "source_original_design_inclusion_probability":cell["source_pi"],
             "source_Oveture_rail":source}
        out["source_original_cells"].append(row)
        # Only genuine source results examined so far, never substitute unobserved
        # sections with rail-positive ones after looking at biology.
        target.write_text(json.dumps({**out,"status":"PARTIAL-not-full-32-cell-source-result"},
                                     ensure_ascii=False,indent=2),encoding="utf8")
        print(json.dumps({"continent":continent,"n_cells_processed":i+1,
             "source_cell":row["source_original_grid_ix_iy"],
             "source_rail_status":source["status"],
             "source_rail_km":source["source_rail_geodesic_km_in_original_cell"]}),flush=True)
    out["status"]="COMPLETE-continent-source-partition"
    out["source_n_rail_mapped_present"]=sum(
        x["source_Oveture_rail"]["source_rail_candidate_in_original_land_cell"] is True for x in out["source_original_cells"])
    out["source_n_rail_mapped_absent"]=sum(
        x["source_Oveture_rail"]["source_rail_candidate_in_original_land_cell"] is False for x in out["source_original_cells"])
    out["source_n_rail_unavailable"]=sum(
        not x["source_Oveture_rail"]["source_complete"] for x in out["source_original_cells"])
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    return out

def consolidate(frame_path,directory,reference):
    j=json.loads(Path(frame_path).read_text(encoding="utf8"))
    selected=selected_cells(j)
    baseline=json.loads(Path(reference).read_text(encoding="utf8"))
    if baseline["original_land_grid_frame_sha256"]!=SOURCE_FRAME_SHA:
        raise ValueError("v352 genuine first-wave source baseline SHA changed")
    first={r["source_grid"]:r["source_mapped_rail"] for r in baseline["source_rail_candidate_100km_cells"]}
    if len(first)!=6:raise ValueError("Not 6 genuine v352 original rail source baseline cells")
    d=Path(directory)
    source_rows=[]
    for code,continent in CONT_CODES.items():
        files=list(d.glob(f"v354_{code}.json"))
        if len(files)!=1:
            raise FileNotFoundError(f"Original source artifact for {code} missing or duplicate")
        group=json.loads(files[0].read_text(encoding="utf8"))
        if (group.get("source_frame_SHA")!=SOURCE_FRAME_SHA
           or group.get("status")!="COMPLETE-continent-source-partition"
           or group.get("source_continent")!=continent
           or group.get("source_Overture_release")!=RAIL_RELEASE):
            raise ValueError("Original continent Overture source job incomplete or incorrect")
        source_rows+=group["source_original_cells"]
    ids={x["source_original_grid_ix_iy"] for x in source_rows}
    source_ids={r["cell"]["source_original_fixed_grid_index"] for r in selected}
    if len(source_rows)!=32 or len(ids)!=32 or ids!=source_ids:
        raise AssertionError("World source original 32 fixed sample cells changed/duplicated/dropped")
    for row in source_rows:
        x=row["source_original_grid_ix_iy"]
        if x in first:
            now=row["source_Oveture_rail"]
            old=first[x]
            if now["source_rail_candidate_in_original_land_cell"]!=old["source_rail_candidate_in_original_land_cell"]:
                raise AssertionError("Same fixed Overture original source v352 rail present indicator drift")
            if not now["source_complete"]:
                raise AssertionError("Previous genuine original v352 source cell failed v354 verification")
            if abs(now["source_rail_geodesic_km_in_original_cell"]-old["source_rail_geodesic_km_in_original_cell"])>0.001:
                raise AssertionError("Same fixed source release v352 rail geodesic length baseline drift")
    summary=[]
    for cont in CONT_CODES.values():
        for band in ("cool","temperate","warm"):
            rows=[r for r in source_rows if r["continent"]==cont and r["source_climate_band"]==band]
            eligible=[r["source_Oveture_rail"] for r in rows if r["source_Oveture_rail"]["source_complete"]]
            have=[r for r in eligible if r["source_rail_candidate_in_original_land_cell"]]
            summary.append({"continent":cont,"source_temp_band":band,
              "n_original_sample_cells":len(rows),"source_rail_complete_n":len(eligible),
              "source_rail_unavailable_n":len(rows)-len(eligible),
              "source_rail_present_n":len(have),
              "source_rail_absent_n":sum(x["source_rail_candidate_in_original_land_cell"] is False for x in eligible),
              "source_rail_total_clipped_km_of_sampled_cells":round(
                  sum(x["source_rail_geodesic_km_in_original_cell"] for x in eligible),3)})
    return {"version":"v354-real-global-probability-landcell-all32-Overture-network-source",
       "original_v351_global_frame_sha256":SOURCE_FRAME_SHA,
       "original_v352_six_source_baseline_reproduced":True,
       "Overture_source_release":RAIL_RELEASE,
       "original_32_sample_land_cells_included_without_outcome_selection":True,
       "total_original_source_landcells_n":32,
       "original_source_complete_cell_n":sum(r["source_Oveture_rail"]["source_complete"] for r in source_rows),
       "original_source_unavailable_cell_n":sum(not r["source_Oveture_rail"]["source_complete"] for r in source_rows),
       "source_rail_present_n":sum(r["source_Oveture_rail"]["source_rail_candidate_in_original_land_cell"] is True for r in source_rows),
       "source_rail_absent_n":sum(r["source_Oveture_rail"]["source_rail_candidate_in_original_land_cell"] is False for r in source_rows),
       "original_rail_km_only_within_selected_landcells":round(sum(
           (r["source_Oveture_rail"]["source_rail_geodesic_km_in_original_cell"] or 0) for r in source_rows),3),
       "continent_macrotemperature_strata":summary,
       "source_rows":source_rows,
       "not_global_network_length_weighted":True,
       "no_train_activity_or_railway_plant_warming_inferred":True,
       "no_original_source_rail_geometries_or_plant_geo_IDs_republished":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--frame",required=True)
    p.add_argument("--continent",choices=list(CONT_CODES.values()))
    p.add_argument("--output",required=True)
    p.add_argument("--aggregate-input-dir")
    p.add_argument("--v352-baseline",default=WAVE1_REFERENCE)
    a=p.parse_args()
    if a.aggregate_input_dir:
        result=consolidate(a.frame,a.aggregate_input_dir,a.v352_baseline)
        target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
        print(json.dumps({k:result[k] for k in (
          "total_original_source_landcells_n","original_source_complete_cell_n",
          "original_source_unavailable_cell_n","source_rail_present_n",
          "source_rail_absent_n","original_rail_km_only_within_selected_landcells"
        )},ensure_ascii=False),flush=True)
    else:
        if not a.continent:p.error("--continent is required for original source partition")
        analyze_partition(a.frame,a.continent,a.output)

if __name__=="__main__":main()
