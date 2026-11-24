#!/usr/bin/env python3
"""v351 finite-world land-grid probability frame, selected BEFORE railway/plant results.

Original 100 km × 100 km EPSG:6933 grid-center source population, only
Natural Earth admin-0 land continents and official WorldClim BIO1 climate.
2 without replacement within each of 18 continent × macroclimate strata.
No railway-biological or citizen-source observation is consulted.
See docs/tetsu_v351_global_equal_area_worldclim_land_frame_frozen_ja.md.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,zipfile
from collections import Counter,defaultdict
from pathlib import Path

WORLDCLIM_SHA="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
NATURAL_EARTH_COMMIT="ca96624a56bd078437bca8184e78163e5039ad19"
NATURAL_EARTH_URL=("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
                   +NATURAL_EARTH_COMMIT+"/geojson/ne_110m_admin_0_countries.geojson")
WORLDCLIM_URL="https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip"
CELL_M=100_000
X_START=-17_000_000
X_END=17_000_000
Y_START=-7_300_000
Y_END=7_300_000
SEED="tetsu-v351-20261124-global-landcell-temperature-before-rail-plant-observations"
CONTINENTS=("Africa","Asia","Europe","North America","South America","Oceania")
TEMP_BANDS=("cool","temperate","warm")
DRAW_PER_STRATUM=2

def sha_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()

def temp_band(celsius):
    if not math.isfinite(celsius):return None
    return "cool" if celsius<8 else "temperate" if celsius<20 else "warm"

def sample_exact_without_replacement(frame,k,seed):
    """Source frame independent of all rail or plant detections.

    An individual cell has exact pi=min(k,N)/N under this finite SRSWOR.
    """
    if len(frame)!=len(set(c["cell_ix_iy"] for c in frame)):
        raise ValueError("Duplicate cells in original sampling frame")
    ordered=sorted(frame,key=lambda c:c["cell_ix_iy"])
    n=len(ordered); drawn=min(k,n)
    if n==0:return [],None
    rng=random.Random(hashlib.sha256(seed.encode("utf8")).hexdigest())
    return rng.sample(ordered,drawn),drawn/n

def source_continents(path):
    from shapely.geometry import shape
    from shapely.ops import unary_union
    from shapely.prepared import prep
    from shapely.validation import make_valid
    from shapely.geometry import Polygon, MultiPolygon
    doc=json.loads(Path(path).read_text(encoding="utf8"))
    if doc.get("type")!="FeatureCollection" or not isinstance(doc.get("features"),list):
        raise ValueError("Not original Natural Earth country GeoJSON")
    geom=defaultdict(list);qa=Counter()
    for feat in doc["features"]:
        prop=feat.get("properties") or {}
        cont=prop.get("CONTINENT")
        qa["original_countries"]+=1
        if cont not in CONTINENTS:
            qa["excluded_antarctica_or_other"]+=1
            continue
        g=shape(feat.get("geometry"))
        if g.is_empty:raise ValueError("Empty original Natural Earth country geometry")
        if not g.is_valid:
            repaired=make_valid(g)
            polygons=[]
            def collect_valid_polys(part):
                if isinstance(part,Polygon):polygons.append(part)
                elif isinstance(part,MultiPolygon):
                    polygons.extend(part.geoms)
                elif hasattr(part,"geoms"):
                    for sub in part.geoms:collect_valid_polys(sub)
            collect_valid_polys(repaired)
            if not polygons:raise ValueError("Unrepairable invalid Natural Earth land polygon")
            new=unary_union(polygons)
            if not new.is_valid:raise ValueError("Repaired Natural Earth polygon is still invalid")
            # A topological repair can change small coastal edge area but cannot
            # be chosen based on Overture rail or plant observations.
            area_drift=abs(new.area-g.area)/max(abs(g.area),1e-15)
            if area_drift>0.01:
                raise ValueError("Natural Earth repair changed polygon area by >1%")
            qa["source_invalid_country_polygons_repaired"]+=1
            qa["source_max_polygon_area_relative_drift_ppm"] = max(
                qa["source_max_polygon_area_relative_drift_ppm"],
                round(area_drift*1e6,3))
            g=new
        geom[cont].append(g)
        qa["country_"+cont]+=1
    if set(geom)!=set(CONTINENTS):raise ValueError("Not all six continents represented by source")
    return {k:prep(unary_union(v)) for k,v in geom.items()},dict(qa)

def worldclim_tif(zip_path):
    f=zipfile.ZipFile(zip_path)
    all_names=[x.filename for x in f.infolist() if x.filename.lower().endswith(".tif")]
    if len(all_names)!=19:
        f.close();raise ValueError("Original WorldClim must contain exactly 19 TIFF BIO values")
    bio1=[x for x in all_names if x.endswith("bio_1.tif")]
    if len(bio1)!=1:
        f.close();raise ValueError("Cannot uniquely identify official WorldClim BIO1 source")
    from rasterio.io import MemoryFile
    raw=f.read(bio1[0]); f.close()
    mem=MemoryFile(raw)
    raster=mem.open()
    if raster.crs is None or raster.count!=1:
        raster.close();mem.close()
        raise ValueError("Original WorldClim BIO1 missing spatial units or unexpected channels")
    if raster.crs.to_epsg()!=4326:
        raster.close();mem.close()
        raise ValueError("WorldClim BIO1 expected EPSG:4326 latlon source raster")
    return raster,mem,{"original_bioclimatic_TIFF_count":19,"source_BIO1_archive_member":bio1[0],
                      "source_BIO1_bytes":len(raw),"source_BIO1_sha256":hashlib.sha256(raw).hexdigest(),
                      "worldclim_crs":str(raster.crs),"worldclim_raster_shape":list(raster.shape),
                      "worldclim_original_nodata":raster.nodata}

def read_temp(src,lat,lon):
    value=next(src.sample([(lon,lat)],masked=True))[0]
    try:
        if bool(getattr(value,"mask",False)):return None
        v=float(value)
    except (TypeError,ValueError):return None
    if not math.isfinite(v) or (src.nodata is not None and math.isclose(v,src.nodata)):return None
    return v

def climate_scale(src):
    # Physical BIO1 units checked against independent distant geographic controls.
    controls={"Tokyo":read_temp(src,35.6812,139.7671),
              "Prague":read_temp(src,50.0839,14.4350)}
    t,p=controls["Tokyo"],controls["Prague"]
    if t is None or p is None:raise ValueError("BIO1 source controls at Tokyo and Prague missing")
    source_in_c=(5.0<=t<=35.0 and -5.0<=p<=20.0)
    source_x10=(50.0<=t<=350.0 and -50.0<=p<=200.0)
    if source_in_c==source_x10:
        raise ValueError("BIO1 original unit scaling unresolved; cannot classify climate bands")
    scale=1.0 if source_in_c else 10.0
    return scale,{"original_Tokyo_BIO1":t,"original_Prague_BIO1":p,
                  "BI01_physical_T_conversion_divisor":scale,
                  "two_independent_control_locations_source_C":
                     {"Tokyo":t/scale,"Prague":p/scale}}

def cell_population(land,src,scale):
    from shapely.geometry import Point
    from pyproj import Transformer
    transform=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    strata=defaultdict(list);qa=Counter()
    ix_count=(X_END-X_START)//CELL_M
    iy_count=(Y_END-Y_START)//CELL_M
    for ix in range(ix_count):
        x=X_START+(ix+0.5)*CELL_M
        for iy in range(iy_count):
            y=Y_START+(iy+0.5)*CELL_M
            lon,lat=transform.transform(x,y)
            qa["total_global_equal_area_cell_centers"]+=1
            if not (math.isfinite(lon) and math.isfinite(lat) and -90<=lat<=90 and -180<=lon<=180):
                qa["invalid_global_proj"]+=1;continue
            p=Point(lon,lat)
            hits=[continent for continent in CONTINENTS if land[continent].covers(p)]
            if len(hits)!=1:
                qa["center_ocean_or_continent_ambiguous"]+=1;continue
            continent=hits[0];qa["source_continent_on_land_cells"]+=1
            t=read_temp(src,lat,lon)
            if t is None:
                qa["source_land_without_WorldClim_BIO1"]+=1;continue
            temp=t/scale
            if not -70<=temp<=50:raise ValueError("WorldClim temperature outside physical plausible bounds")
            band=temp_band(temp)
            qa["source_climate_supported_land_cells"]+=1
            strata[(continent,band)].append({
              "cell_ix_iy":f"{ix}:{iy}",
              "cell_center_lat":lat,"cell_center_lon":lon,
              "source_BIO1_C":temp,
              "continent":continent,"climate_band":band
            })
    if qa["source_climate_supported_land_cells"]<2000:
        raise ValueError("Expected >=2000 global source land-grid cells; partial world frame rejected")
    return strata,dict(qa)

def summarize(strata,qa,meta):
    out=[]
    for continent in CONTINENTS:
        for band in TEMP_BANDS:
            available=strata.get((continent,band),[])
            selected,pi=sample_exact_without_replacement(
                available,DRAW_PER_STRATUM,SEED+"|"+continent+"|"+band)
            rec={
              "continent":continent,"source_macrotemperature_band_C":band,
              "n_equal_area_source_land_climate_cells":len(available),
              "draws_source_cells":len(selected),
              "source_finite_population_inclusion_probability_per_cell":pi,
              "source_sampling_eligible_only_if_land_center_has_BIO1":True,
              "sampled_cells":[{
                "source_grid_cell_sha256":hashlib.sha256(x["cell_ix_iy"].encode()).hexdigest(),
                "source_original_fixed_grid_index":x["cell_ix_iy"],
                "coarse_sample_cell_center_lat":round(x["cell_center_lat"],2),
                "coarse_sample_cell_center_lon":round(x["cell_center_lon"],2),
                "original_macroclimate_BIO1_C":round(x["source_BIO1_C"],3),
                "source_has_any_railway":None,
                "source_iNaturalist_plant_occurrences":None
              } for x in selected]
            }
            out.append(rec)
    if sum(row["n_equal_area_source_land_climate_cells"] for row in out)!=qa["source_climate_supported_land_cells"]:
        raise AssertionError("Stratum population count mismatch")
    centers=[x["source_original_fixed_grid_index"] for row in out for x in row["sampled_cells"]]
    if len(centers)!=len(set(centers)):raise AssertionError("Sampled cells duplicated")
    return {"version":"v351-real-global-land-grid-worldclim-population-source",
            "source_input_QA":meta,"source_full_cell_grid_QA":qa,
            "strata":out,"source_strata_n":len(out),
            "n_sampled_unique_land_cells":len(centers),
            "n_original_continents":len(CONTINENTS),
            "sampled_cell_unique_hash":hashlib.sha256("\n".join(sorted(centers)).encode()).hexdigest(),
            "original_railway_geometries_not_yet_tested":True,
            "original_plant_observations_not_yet_tested":True,
            "not_an_actual_global_rail_length_probability_frame":True,
            "ecological_rail_plant_thermal_effects_not_estimated":True,
            "no_original_precise_plant_coordinates_or_observers_republished":True,
            "sampled_cell_centres_are_100km_grid_centroids_not_rare_species_observations":True}

def analyze(naturalearth,worldclim):
    land,landqa=source_continents(naturalearth)
    checksum=sha_file(worldclim)
    if checksum!=WORLDCLIM_SHA:raise ValueError("WorldClim official source SHA256 mismatch")
    raster,mem,worldqa=worldclim_tif(worldclim)
    try:
        scale,unitqa=climate_scale(raster)
        cells,frameqa=cell_population(land,raster,scale)
    finally:raster.close();mem.close()
    meta={"Natural_Earth":{
            "original_repository_commit":NATURAL_EARTH_COMMIT,
            "source_url":NATURAL_EARTH_URL,
            "source_file_sha256":sha_file(naturalearth),**landqa},
          "WorldClim_BIO1":{"source_url":WORLDCLIM_URL,
            "source_zip_sha256":checksum,**worldqa,**unitqa},
          "equal_area_crs":"EPSG:6933",
          "sampling_frame_cell_width_m":CELL_M,
          "original_frame_x_m":[X_START,X_END],
          "original_frame_y_m":[Y_START,Y_END],
          "sampling_draws_per_available_continent_temperature_stratum":DRAW_PER_STRATUM,
          "source_seed":SEED}
    return summarize(cells,frameqa,meta)

def main():
    p=argparse.ArgumentParser()
    for key in ("natural-earth","worldclim-zip","output"):
        p.add_argument("--"+key,required=True)
    a=p.parse_args()
    result=analyze(a.natural_earth,a.worldclim_zip)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
