#!/usr/bin/env python3
"""v352: Real railway presence in SIX probability-selected global land-grid cells.

Selected original cells are a deterministic temperate-first wave from v351
BEFORE any Overture railway or iNaturalist observation results were queried.
Original source: Overture 2026-11-18 transportation segments, filtered rail,
clipped to EXACT original EPSG:6933 100 km square. Geo/plant raw never shared.
"""
from __future__ import annotations
import argparse,hashlib,importlib.metadata,json,math,time,urllib.parse,urllib.request
from collections import Counter
from pathlib import Path

from generate_global_equal_area_climate_land_grid_v351 import (
    CELL_M,X_START,Y_START,CONTINENTS,SEED,WORLDCLIM_SHA
)
from audit_Overture_global_rail_segments_v349 import eligible_segment

SOURCE_FRAME_SHA="b65a13e32d65e8fa818bbefe2d7fd1f0e15c6c5cb97ff6f6212aea23cf138bf2"
NATURAL_EARTH_SHA="718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b"
RAIL_RELEASE="2026-11-18.0"
SOURCE_MAX_TRANSPORT_ROWS=3_000_000
SOURCE_STAC_ERROR_BACKOFF=(8,30)
SOURCE_I_NAT="https://api.inaturalist.org/v1/observations"

def wave_a(frame):
    if frame["sampled_cell_unique_hash"]!=SOURCE_FRAME_SHA:
        raise ValueError("v351 original probability sample SHA mismatch")
    if frame["source_input_QA"]["Natural_Earth"]["source_file_sha256"]!=NATURAL_EARTH_SHA:
        raise ValueError("v351 actual Natural Earth map-units SHA mismatch")
    if frame["source_input_QA"]["WorldClim_BIO1"]["source_zip_sha256"]!=WORLDCLIM_SHA:
        raise ValueError("v351 WorldClim BIO1 source SHA mismatch")
    if frame["source_full_cell_grid_QA"]["source_climate_supported_land_cells"]!=13437:
        raise ValueError("v351 world cell population drift")
    selected=[]
    for cont in CONTINENTS:
        arr=[r for r in frame["strata"] if r["continent"]==cont and
             r["source_macrotemperature_band_C"]=="temperate"]
        if len(arr)!=1 or len(arr[0]["sampled_cells"])!=2:
            raise ValueError("Frozen v351 worldwide temperate stratum missing")
        cells=arr[0]["sampled_cells"]
        chosen=min(cells,key=lambda x:x["source_grid_cell_sha256"])
        selected.append({"continent":cont,"source_parent_population_n":arr[0]["n_equal_area_source_land_climate_cells"],
          "source_wave_A_cell_inclusion_probability":1/arr[0]["n_equal_area_source_land_climate_cells"],
          **chosen})
    if len({c["source_original_fixed_grid_index"] for c in selected})!=6:
        raise ValueError("Frozen wave A original sampled cells not distinct")
    return selected

def projected_cell_square(cell, margin_deg=0.001):
    from pyproj import Transformer
    ix,iy=map(int,cell["source_original_fixed_grid_index"].split(":"))
    x0=X_START+ix*CELL_M
    y0=Y_START+iy*CELL_M
    inv=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    from shapely.geometry import box
    geom=box(x0,y0,x0+CELL_M,y0+CELL_M)
    corners=[inv.transform(x,y) for x in (x0,x0+CELL_M/2,x0+CELL_M)
                             for y in (y0,y0+CELL_M/2,y0+CELL_M)]
    lons=[v[0] for v in corners];lats=[v[1] for v in corners]
    if max(lons)-min(lons)>=180:raise ValueError("Original selected cell crosses dateline")
    b=(min(lons)-margin_deg,min(lats)-margin_deg,max(lons)+margin_deg,max(lats)+margin_deg)
    if not(-180<=b[0]<b[2]<=180 and -90<=b[1]<b[3]<=90):
        raise ValueError("Original selected cell WGS bbox invalid")
    return geom,b

def clipped_geodesic_rail_length(geom,source_poly,forward,inverse,geod):
    from shapely.ops import transform
    transformed=transform(forward.transform,geom)
    clipped=transformed.intersection(source_poly)
    if clipped.is_empty:return 0.0,0
    candidates=[]
    def walk(g):
        if g.geom_type=="LineString":
            if len(g.coords)>=2 and g.length>0:candidates.append(g)
        elif hasattr(g,"geoms"):
            for sub in g.geoms:walk(sub)
    walk(clipped)
    return sum(abs(geod.geometry_length(transform(inverse.transform,g))) for g in candidates),len(candidates)

def source_overture(cell):
    import pyarrow.compute as pc
    from overturemaps import record_batch_reader
    from shapely import wkb
    from pyproj import Geod,Transformer
    from overturemaps.core import StacIndexWarning
    import warnings
    poly,bbox=projected_cell_square(cell)
    forward=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    inverse=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    geod=Geod(ellps="WGS84")
    errors=[]
    for attempt in range(1,4):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error",StacIndexWarning)
                reader=record_batch_reader("segment",bbox=bbox,release=RAIL_RELEASE,stac=True)
                if reader is None:
                    return {"status":"complete_no_source_transportation_segments",
                        "source_rail_features_in_original_cell":0,
                        "source_rail_geodesic_km_in_original_cell":0.0,
                        "source_complete":True,
                        "source_rail_candidate_in_original_land_cell":False,
                        "source_total_transportation_records_examined":0,
                        "original_source_attempt":attempt}
                n=Counter();rail_m=0.0
                for batch in reader:
                    n["source_transportation_rows"]+=batch.num_rows
                    if n["source_transportation_rows"]>SOURCE_MAX_TRANSPORT_ROWS:
                        raise ValueError("Source transport row cap exceeded; cell unidentifiable")
                    idx=batch.schema.get_field_index("subtype")
                    if idx<0:raise ValueError("Original Overture subtype missing")
                    matched=batch.filter(pc.equal(batch.column(idx),"rail"))
                    n["source_rail_tagged_before_class_filters"]+=matched.num_rows
                    if not matched.num_rows:continue
                    for row in matched.to_pylist():
                        good,reason=eligible_segment(row)
                        if not good:
                            n["source_rail_filter_reject_"+str(reason)]+=1
                            continue
                        geom=wkb.loads(bytes(row["geometry"]))
                        metres,numsegments=clipped_geodesic_rail_length(
                            geom,poly,forward,inverse,geod)
                        if metres>0:
                            rail_m+=metres
                            n["source_rail_eligible_clipped_features"]+=1
                            n["source_rail_eligible_clipped_parts"]+=numsegments
            if rail_m<0:raise ValueError("Negative source geodesic track length")
            return {"status":"original_overture_source_complete",
                    "source_rail_features_in_original_cell":n["source_rail_eligible_clipped_features"],
                    "source_rail_geodesic_km_in_original_cell":round(rail_m/1000,3),
                    "source_complete":True,
                    "source_rail_candidate_in_original_land_cell":rail_m>0,
                    "source_total_transportation_records_examined":n["source_transportation_rows"],
                    "source_rail_tagged_before_class_filters":n["source_rail_tagged_before_class_filters"],
                    "source_rail_clipped_parts":n["source_rail_eligible_clipped_parts"],
                    "source_filter_rejections":{k:v for k,v in n.items()
                        if k.startswith("source_rail_filter_reject_")},
                    "original_source_attempt":attempt,
                    "source_WGS84_query_bbox_rounded":list(map(lambda v:round(v,3),bbox)),
                    "original_railway_source_not_confirmed_train_activity":True}
        except Exception as e:
            errors.append({"attempt":attempt,"error_type":type(e).__name__,
                           "short_error":str(e)[:190]})
            if attempt<=2:time.sleep(SOURCE_STAC_ERROR_BACKOFF[attempt-1])
    return {"status":"railway_source_unavailable_NOT_zero_rail_or_plants",
            "source_complete":False,"source_rail_candidate_in_original_land_cell":None,
            "source_rail_features_in_original_cell":None,
            "source_rail_geodesic_km_in_original_cell":None,
            "error_attempts":errors}

def original_inat_count(cell):
    # only public original research-grade wild plant census COUNTS from a
    # rectangular bbox containing the 100km projected cell, not inside it.
    poly,bbox=projected_cell_square(cell)
    west,south,east,north=bbox
    params={"taxon_id":47126,"quality_grade":"research","captive":"false",
       "geo":"true","d1":"2025-01-01","d2":"2025-12-31",
       "swlat":f"{south:.5f}","swlng":f"{west:.5f}",
       "nelat":f"{north:.5f}","nelng":f"{east:.5f}",
       "per_page":1}
    u=SOURCE_I_NAT+"?"+urllib.parse.urlencode(params)
    try:
        req=urllib.request.Request(u,headers={"User-Agent":"tetsu-v352-global-source-opportunity/1.0"})
        with urllib.request.urlopen(req,timeout=30) as resp:
            body=resp.read(30000)
            obj=json.loads(body)
            count=obj.get("total_results")
            if not isinstance(count,int) or count<0:
                raise ValueError("Original iNat count schema invalid")
            return {"source_2025_research_grade_public_plant_records_bbox_count":count,
                    "source_total_count_only_not_complete_plant_geolocation_census":True,
                    "source_rectangular_bbox_not_exact_equal_area_cell":True,
                    "source_HTTP":resp.status}
    except Exception as e:
        return {"source_2025_research_grade_public_plant_records_bbox_count":None,
           "status":"original_iNat_source_count_unavailable",
           "source_error_type":type(e).__name__,
           "source_short_message":str(e)[:155]}

def run(source_frame_path,checkpoint=None):
    j=json.loads(Path(source_frame_path).read_text(encoding="utf8"))
    chosen=wave_a(j)
    rows=[]
    for cell in chosen:
        info=source_overture(cell)
        time.sleep(1.25)
        observed=original_inat_count(cell)
        rows.append({"source_continent":cell["continent"],
            "source_original_temperature_band":"temperate_8_to_20C",
            "source_BIO1_original_climate_C":cell["original_macroclimate_BIO1_C"],
            "source_original_cell_grid_index":cell["source_original_fixed_grid_index"],
            "source_original_cell_sha256":cell["source_grid_cell_sha256"],
            "source_coarse_100km_center":[cell["coarse_sample_cell_center_lat"],
                                         cell["coarse_sample_cell_center_lon"]],
            "source_parent_landcell_population_n":cell["source_parent_population_n"],
            "source_design_cell_inclusion_probability":cell["source_wave_A_cell_inclusion_probability"],
            "original_Overture":info,"original_iNaturalist_source_count_only":observed})
        print(json.dumps({"source_progress_continent":cell["continent"],
            "source_rail_status":info["status"],
            "source_clipped_rail_km":info["source_rail_geodesic_km_in_original_cell"],
            "source_original_iNat_bbox_records":
              observed["source_2025_research_grade_public_plant_records_bbox_count"]}),
              flush=True)
        if checkpoint is not None:
            partial={"version":"v352-original-partial-source-availability",
                "n_source_cells_processed":len(rows),"n_source_cells_expected":len(chosen),
                "source_results_NOT_GLOBAL_WHEN_PARTIAL":True,
                "rows":rows}
            Path(checkpoint).parent.mkdir(parents=True,exist_ok=True)
            Path(checkpoint).write_text(json.dumps(partial,ensure_ascii=False,indent=2))
    return {"version":"v352-world-source-random-land-grid-waveA-rail-presence",
        "original_frame_sample_sha256":SOURCE_FRAME_SHA,
        "original_Overture_release":RAIL_RELEASE,
        "original_Overture_client_version":importlib.metadata.version("overturemaps"),
        "original_source_cells_preselected_without_using_rail_or_plants":True,
        "n_source_waveA_cells":len(rows),
        "n_original_rail_present_cells":sum(x["original_Overture"]["source_rail_candidate_in_original_land_cell"] is True for x in rows),
        "n_original_rail_absent_cells_source_verified":sum(x["original_Overture"]["source_rail_candidate_in_original_land_cell"] is False for x in rows),
        "n_rail_source_unavailable_cells":sum(not x["original_Overture"]["source_complete"] for x in rows),
        "rows":rows,"original_public_iNat_source_is_count_only_no_spatial_classes":True,
        "not_an_estimate_of_global_railway_network_length":True,
        "no_plant_warmth_or_trackbed_ecological_effect_estimated":True,
        "no_original_precise_plant_geolocation_observer_species_raw_or_railway_geometries_republished":True}

def main():
    a=argparse.ArgumentParser();a.add_argument("--frame",required=True);a.add_argument("--output",required=True)
    p=a.parse_args()
    result=run(p.frame,checkpoint=p.output)
    dst=Path(p.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
