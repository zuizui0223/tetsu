#!/usr/bin/env python3
"""v353 actual global 100km finite land-cell source iNat plant observation counts.

Count the 2025 public research-grade plant records reported by official
iNaturalist API in WGS bbox enclosing EACH sampled projected 100km cell.
No source individual observation IDs, names or exact GPS returned or stored.
This is observation OPPORTUNITY, not actual botanical plot vegetation.
"""
from __future__ import annotations
import argparse,hashlib,json,statistics,time,urllib.parse,urllib.request
from pathlib import Path
from collections import defaultdict,Counter
from urllib.error import HTTPError,URLError
from audit_global_landcell_Overture_rail_v352 import (
    SOURCE_FRAME_SHA,projected_cell_square
)
from generate_global_equal_area_climate_land_grid_v351 import CONTINENTS,TEMP_BANDS

INAT="https://api.inaturalist.org/v1/observations"
MIN_DELAY_S=1.1
SOURCE_YEAR=2025
LIMIT_SOURCE_RECORDS_PER_PAGE=1

def selected_cells(frame):
    if frame.get("sampled_cell_unique_hash")!=SOURCE_FRAME_SHA:
        raise ValueError("Original source world probability sample not pinned")
    rows=[]
    for stratum in frame["strata"]:
        for cell in stratum["sampled_cells"]:
            rows.append({"continent":stratum["continent"],
                        "macrotemp_stratum":stratum["source_macrotemperature_band_C"],
                        "population_n":stratum["n_equal_area_source_land_climate_cells"],
                        "source_pi":stratum["source_finite_population_inclusion_probability_per_cell"],
                        "cell":cell})
    if len(rows)!=32 or len({row["cell"]["source_original_fixed_grid_index"] for row in rows})!=32:
        raise ValueError("World source SRSWOR 32 original selected land cells drifted")
    return rows

def count_source_cell(row,client=None):
    cell=row["cell"]
    _,bbox=projected_cell_square(cell)
    w,s,e,n=bbox
    params={"taxon_id":47126,"quality_grade":"research","captive":"false",
        "geo":"true","d1":"2025-01-01","d2":"2025-12-31",
        "swlat":f"{s:.5f}","swlng":f"{w:.5f}",
        "nelat":f"{n:.5f}","nelng":f"{e:.5f}",
        "per_page":LIMIT_SOURCE_RECORDS_PER_PAGE}
    url=INAT+"?"+urllib.parse.urlencode(params)
    errors=[]
    for attempt in (1,2):
        try:
            if client is None:
                req=urllib.request.Request(url,headers={"User-Agent":"tetsu-v353-sampled-global-land-grid-plants/1.0",
                    "Accept":"application/json"})
                with urllib.request.urlopen(req,timeout=35) as response:
                    data=response.read(2_000_001)
                    if len(data)>2_000_000:
                        raise ValueError("Original one-record iNat JSON >2MB source bound")
                    status=response.status
            else:
                data,status=client(url)
            if status!=200:raise ValueError("Original public iNaturalist API HTTP!=200")
            obj=json.loads(data)
            count=obj.get("total_results")
            if not isinstance(count,int) or count<0 or not isinstance(obj.get("results"),list):
                raise ValueError("Original iNaturalist source API count schema unexpected")
            return {"source_api_2025_public_plants_bbox_count":count,
                    "source_HTTP":status,
                    "source_response_SHA256":hashlib.sha256(data).hexdigest(),
                    "source_original_results_per_page":LIMIT_SOURCE_RECORDS_PER_PAGE,
                    "source_rectangular_bbox_not_exact_landcell_geometry":True,
                    "source_original_observations_QUALITY_filtered_with_per_record_GPS":False,
                    "source_plant_records_not_survey_plot_absences":True,
                    "source_API_attempts":attempt}
        except (HTTPError,URLError,TimeoutError,ValueError,OSError,json.JSONDecodeError) as err:
            errors.append({"attempt":attempt,"type":type(err).__name__,
                          "source_http":getattr(err,"code",None),
                          "short_reason":str(err)[:145]})
            if attempt==1:time.sleep(10.0)
    return {"source_api_2025_public_plants_bbox_count":None,
            "source_status":"unavailable_not_zero_biology",
            "source_errors":errors}

def summarize(rows):
    out=[]
    for cont in CONTINENTS:
        for band in TEMP_BANDS:
            selected=[r for r in rows if r["continent"]==cont and r["source_climate_stratum"]==band]
            good=[r["API"]["source_api_2025_public_plants_bbox_count"] for r in selected
                  if r["API"]["source_api_2025_public_plants_bbox_count"] is not None]
            out.append({"continent":cont,"macrotemp_band":band,
               "source_probability_draws":len(selected),
               "original_source_available_n":len(good),
               "source_unavailable_n":len(selected)-len(good),
               "source_2025_plant_observations_bbox_sum":sum(good) if good else None,
               "source_2025_plant_observations_bbox_min":min(good) if good else None,
               "source_2025_plant_observations_bbox_max":max(good) if good else None,
               "source_2025_plant_observations_bbox_median":statistics.median(good) if good else None,
               "source_cells_zero_observed_plants_n":sum(z==0 for z in good),
               "source_cells_ge100_observed_plants_n":sum(z>=100 for z in good),
               "source_cells_ge1000_observed_plants_n":sum(z>=1000 for z in good),
               "source_cells_gt10000_public_API_paging_cap_n":sum(z>10000 for z in good)})
    return out

def run(frame_path,checkpoint=None):
    frame=json.loads(Path(frame_path).read_text(encoding="utf8"))
    all_selected=selected_cells(frame)
    rows=[]
    for i,row in enumerate(all_selected):
        if i>0:time.sleep(MIN_DELAY_S)
        api=count_source_cell(row)
        cell=row["cell"]
        rows.append({"continent":row["continent"],
          "source_climate_stratum":row["macrotemp_stratum"],
          "source_grid_ID":cell["source_original_fixed_grid_index"],
          "source_original_grid_SHA256":cell["source_grid_cell_sha256"],
          "source_grid_WorldClim_BIO1_C":cell["original_macroclimate_BIO1_C"],
          "source_original_cell_pi":row["source_pi"],
          "API":api})
        print(json.dumps({"source_progress_n":i+1,"source_continent":row["continent"],
             "source_climate_stratum":row["macrotemp_stratum"],
             "source_API_2025_plant_observations_bbox":api["source_api_2025_public_plants_bbox_count"]}),
             flush=True)
        if checkpoint:
            dst=Path(checkpoint);dst.parent.mkdir(parents=True,exist_ok=True)
            dst.write_text(json.dumps({"version":"v353-PARTIAL-not-full-world-census",
              "source_successful_original_cells_processed":len(rows),"rows":rows,
              "no_global_inference_if_partial":True},ensure_ascii=False,indent=2),encoding="utf8")
    summary=summarize(rows)
    return {"version":"v353-original-iNat-source-2025-plant-opportunity-in-random-global-landcells",
      "source_WorldClim_NaturalEarth_sample_frame_SHA":SOURCE_FRAME_SHA,
      "source_NaturalEarth_world_temperatures_continent_strata_n":18,
      "source_100km_equal_area_landcells_sampled_n":32,
      "source_year":SOURCE_YEAR,
      "source_known_public_plant_observations_sum_only_in_rectangular_bboxes":
        sum(r["API"]["source_api_2025_public_plants_bbox_count"] or 0 for r in rows),
      "source_count_available_n":sum(r["API"]["source_api_2025_public_plants_bbox_count"] is not None for r in rows),
      "source_cells_zero_iNat_public_observations_n":sum(r["API"]["source_api_2025_public_plants_bbox_count"]==0 for r in rows),
      "source_cells_with_over_100_public_observations_n":sum(
        (r["API"]["source_api_2025_public_plants_bbox_count"] or 0)>=100 for r in rows),
      "original_Plant_API_is_observation_census_NOT_standardized_field_plots":True,
      "original_source_rectangular_query_encloses_100km_projected_cell_not_equal_area_count":True,
      "railway_source_not_used_for_sampling_or_plant_access_inclusion":True,
      "source_BIO1_macrotemperature_is_NOT_global_plant_species_thermal_niche":True,
      "estimated_global_railway_warm_species_effect":None,
      "no_original_observer_plant_GPS_or_taxonomic_identity_republished":True,
      "strata_source_counts":summary,"source_original_cells":rows}

def main():
    p=argparse.ArgumentParser();p.add_argument("--frame",required=True);p.add_argument("--output",required=True)
    args=p.parse_args();obj=run(args.frame,checkpoint=args.output)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({"version":obj["version"],
        "total_original_public_source_count":obj["source_known_public_plant_observations_sum_only_in_rectangular_bboxes"],
        "original_sampling_cells":obj["source_100km_equal_area_landcells_sampled_n"],
        "source_count_available_n":obj["source_count_available_n"],
        "source_zeros_n":obj["source_cells_zero_iNat_public_observations_n"],
        "source_ge100_n":obj["source_cells_with_over_100_public_observations_n"],
        "strata":obj["strata_source_counts"]},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
