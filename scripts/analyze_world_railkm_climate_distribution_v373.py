#!/usr/bin/env python3
"""v373 genuine global railway-km macroclimate exposure census.

Full Overture 2026-11 railway 1,536,226 normal candidate source segments
and 2,237,043km aggregated by true 100/200km EPSG6933 equal-area cells.
Original WorldClim2.1 BIO1 1970-2000 + NaturalEarth land definitions.
This is the WORLDWIDE PHYSICAL CLIMATE DISTRIBUTION OF RAIL INFRASTRUCTURE,
not an ecological species thermal niche, railway temperature sensor reading,
train seed dispersal or railroad-versus-control biotic effect.

Necessary confounding baseline for original H3 cross-climate plant analyses.
"""
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from shapely.geometry import Point
from pyproj import Transformer

from generate_global_equal_area_climate_land_grid_v351 import (
    source_continents,worldclim_tif,read_temp,climate_scale,
    WORLDCLIM_SHA,CONTINENTS,temp_band,sha_file
)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    load_json,RAIL_SOURCE_KM,RAIL_RELEASE,GRID_WIDTHS
)

WORLDMAP_SHA="718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b"
BANDS=("cool","temperate","warm")

def mapcell(ix,iy,width):
    x=X_START+(ix+.5)*width
    y=Y_START+(iy+.5)*width
    return x,y

def valid_source_exposure(j):
    if j.get("source_Oveture_release")!=RAIL_RELEASE or not j.get(
      "source_original_fullworld_rail_128_files_SHA_verified"):
        raise ValueError("Missing true original entire world 128 GERS source QA")
    if abs(j.get("source_fullworld_candidate_normal_rail_km",0)-RAIL_SOURCE_KM)>.05:
        raise ValueError("Not authentic original 2237043.481km railway centreline frame")
    for k in GRID_WIDTHS:
        data=j.get("source_global_gridded_rail_exposure_km_by_cell_width",{}).get(str(k))
        if not isinstance(data,list) or not data:
            raise ValueError("Original global railway 100km or 200km source missing")
        length=sum(float(x["source_real_original_Overture_normal_rail_km_clipped_into_cell"])
                   for x in data)
        if abs(length-RAIL_SOURCE_KM)>.15:
            raise ValueError("Global original grid railway km sum invalid")
    return True

def weighted_average(pairs):
    sw=sum(w for w,x in pairs)
    return sum(w*x for w,x in pairs)/sw if sw else None

def build(rail,land,raster,scale,width):
    geo=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    rmap={}
    for src in rail:
        index=(int(src["source_equal_area_grid_ix"]),int(src["source_equal_area_grid_iy"]))
        if index in rmap:raise ValueError("Rail source same grid cell duplicated")
        km=float(src["source_real_original_Overture_normal_rail_km_clipped_into_cell"])
        if km<=0 or not math.isfinite(km):raise ValueError("Invalid original rail length km")
        rmap[index]=km
    totals=defaultdict(lambda:defaultdict(lambda:{"landcells":0,
                  "rail_positive_landcells":0,"railkm":0.,
                  "temperatures":[],"weighted_rail_temperatures":[]}))
    globalQA=Counter()
    lat_km=defaultdict(float)
    classified_railkm=0.
    # Every possible source equal-area land-cell center is counted, including
    # mapped rail=0. Do not make denominator only railroad cell opportunity.
    nx=int(34_800_000//width);ny=int(14_800_000//width)
    for ix in range(nx):
        for iy in range(ny):
            x,y=mapcell(ix,iy,width)
            lon,lat=geo.transform(x,y)
            globalQA["all_equal_area_grid_cell_centers"]+=1
            if not(math.isfinite(lon) and math.isfinite(lat)):
                globalQA["invalid_geo_inverse"]+=1;continue
            if not(-180<=lon<=180 and -90<=lat<=90):
                globalQA["outside_earth"]+=1;continue
            hits=[k for k in CONTINENTS if land[k].covers(Point(lon,lat))]
            if len(hits)!=1:
                globalQA["landcenter_ocean_or_uncertain"]+=1
                continue
            temp=read_temp(raster,lat,lon)
            if temp is None:
                globalQA["landcenter_original_BIO1_missing"]+=1
                continue
            physical=temp/scale
            if not -70<=physical<=50:
                raise ValueError("Original WorldClim climate units wrong")
            k=hits[0]
            band=temp_band(physical)
            cell=totals[k][band]
            globalQA["source_center_land_and_climate_supported"]+=1
            cell["landcells"]+=1
            cell["temperatures"].append(physical)
            km=rmap.get((ix,iy),0.)
            if km>0:
                globalQA["landcenter_supported_cells_with_real_original_rail"]+=1
                cell["rail_positive_landcells"]+=1
                cell["railkm"]+=km
                cell["weighted_rail_temperatures"].append((km,physical))
                classified_railkm+=km
    all_rail=sum(rmap.values())
    bycont={}
    for cont in CONTINENTS:
        data=totals.get(cont,{})
        row={}
        for band in BANDS:
            vals=data.get(band,{"landcells":0,"rail_positive_landcells":0,
                               "railkm":0.,"temperatures":[],"weighted_rail_temperatures":[]})
            row[band]={
              "source_total_original_land_cell_centers":vals["landcells"],
              "source_original_land_cell_centers_with_rail":vals["rail_positive_landcells"],
              "source_railkm_in_original_climate_band":round(vals["railkm"],4),
              "source_landcenter_equal_weighted_mean_BIO1_C":mean(vals["temperatures"]) if vals["temperatures"] else None,
              "source_original_railkm_weighted_mean_BIO1_C":weighted_average(vals["weighted_rail_temperatures"])}
        land_n=sum(v["source_total_original_land_cell_centers"] for v in row.values())
        km_n=sum(v["source_railkm_in_original_climate_band"] for v in row.values())
        for band,r in row.items():
            landshare=r["source_total_original_land_cell_centers"]/land_n if land_n else None
            railshare=r["source_railkm_in_original_climate_band"]/km_n if km_n else None
            r["fraction_of_continent_landcenter_cells"]=landshare
            r["fraction_of_continent_railkm"]=railshare
            r["railkm_share_divided_by_landcell_share"]=(railshare/landshare
              if landshare is not None and railshare is not None and landshare>0 else None)
        bycont[cont]=row
    bands={}
    for band in BANDS:
        cells=[totals[cont][band] for cont in CONTINENTS]
        km=sum(z["railkm"] for z in cells)
        n=sum(z["landcells"] for z in cells)
        bands[band]={"source_global_landcenter_cell_n":n,
                     "source_railpositive_landcenter_cell_n":sum(z["rail_positive_landcells"] for z in cells),
                     "source_railkm_in_temperature_band":round(km,4),
                     "railkm_weighted_BIO1_C":weighted_average(
                         [v for z in cells for v in z["weighted_rail_temperatures"]])}
    all_land_n=sum(x["source_global_landcenter_cell_n"] for x in bands.values())
    for band,row in bands.items():
        row["fraction_all_source_landcenters"]=row["source_global_landcenter_cell_n"]/all_land_n if all_land_n else None
        row["fraction_all_climate_classifiable_railkm"]=(
           row["source_railkm_in_temperature_band"]/classified_railkm if classified_railkm else None)
    globalQA["original_full_global_rail_geodesic_km"]=all_rail
    globalQA["source_landcenter_climate_classified_railkm"]=classified_railkm
    globalQA["source_original_railkm_unclassifiable_as_land_center_or_BIO1"]=all_rail-classified_railkm
    globalQA["fraction_global_railkm_with_landcenter_climate"]=classified_railkm/all_rail
    return {"source_grid_width_m":width,"source_full_global_climate_context_QA":dict(globalQA),
        "source_each_of_six_continents_climate_band":bycont,
        "global_temperature_band_rail_vs_land_distributions":bands,
        "source_original_railkm_weighted_BIO1_all_valid_landcenter_cells_C":
             weighted_average([v for c in CONTINENTS for band in BANDS
                for v in totals[c][band]["weighted_rail_temperatures"]]),
        "source_landcell_equalweighted_BIO1_all_valid_landcenter_cells_C":mean([
             x for c in CONTINENTS for band in BANDS for x in totals[c][band]["temperatures"]])}

def analyze(source_rail,land_source,climate_source):
    orig=load_json(source_rail)
    valid_source_exposure(orig)
    if sha_file(land_source)!=WORLDMAP_SHA or sha_file(climate_source)!=WORLDCLIM_SHA:
        raise ValueError("Source worldwide climate or land original checksum differs")
    land,landqa=source_continents(land_source)
    raster,mem,rasqa=worldclim_tif(climate_source)
    try:
        scale,phys=climate_scale(raster)
        result=[build(orig["source_global_gridded_rail_exposure_km_by_cell_width"][str(width)],
                      land,raster,scale,width) for width in GRID_WIDTHS]
    finally:
        raster.close();mem.close()
    return {"version":"v373-original-global-infrastructure-railkm-climate-bias-population-audit",
       "original_GERS_global_release":RAIL_RELEASE,"original_railkm_WGS84_geodesic_km":RAIL_SOURCE_KM,
       "source_worldclim_BIO1_1970_2000_is_original_site_macroclimate_not_species_affinity":True,
       "source_NaturalEarth_map_110m_checksum":WORLDMAP_SHA,
       "source_WorldClim_zip_SHA256":WORLDCLIM_SHA,
       "source_NaturalEarth_QA":landqa,"source_WorldClim_QA":rasqa,
       "source_original_temperature_unit_verification":phys,
       "global_source_grid_rail_climate_results":result,
       "limits":[
         "100/200km source railway length assigned to original cell center BIO1 rather than fine resolved track temperature; coastal cells and islands with ocean center excluded",
         "Overture 2026 conventional candidate railway mapped infrastructure not operated railway schedules/traffic or train-transported seeds",
         "Railkm-weighted vs land-cell weighted temperature difference is a descriptive unequal physical exposure frame, not vegetation functional response",
         "Six continents are not 6 independently sampled flora experiments; geographic temperature and landuse differences govern railway construction",
         "WorldClim mean temperature is 1970–2000 macroscale rather than observed rail ballast microclimate in 2025",
         "Separate 100km/200km spatial supports are sensitivity not independent replication",
         "No plant community dataset processed at this source-only physical infrastructure baseline"]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-v370-rail-grid",required=True)
    p.add_argument("--natural-earth",required=True)
    p.add_argument("--worldclim-zip",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    obj=analyze(a.original_v370_rail_grid,a.natural_earth,a.worldclim_zip)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({"version":obj["version"],"results":obj["global_source_grid_rail_climate_results"]},
                     ensure_ascii=False),flush=True)
if __name__=="__main__":main()
