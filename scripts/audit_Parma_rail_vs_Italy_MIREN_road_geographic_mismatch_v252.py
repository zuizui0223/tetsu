#!/usr/bin/env python3
"""v25.2 original Italy country-reference geography and era falsification.

Revisit latest 2026-11-17 source-authored country_validation.json from
research/validation_20261117 Git commit a8138a50, Git blob SHA pinned.
The only 10-country-reference contrast surviving all eight exploratory
ranking/calibration/sensitivity gates in v25.1 was IT_road.
But the original railway flora is 28 stations **in the province of
Parma** (Adorni & Ghillani 2020 published PDF field survey 2012-13),
whereas the source MIREN 'domestic' road plots are in central Italy
2022, bbox 42.149692..42.47389002 N, 12.89728797..14.13079 E.
Country equal DOES NOT mean region matched or field-year matched.

UN Geographic Names Experts Group official city representative Parma
WGS84: latitude 44.8018611111111, longitude 10.32735
https://ungegn.un.org/dashboard/cities/details?id=1538
This is only PARMA CITY CENTER, NOT coordinates of all 28 original
provincial station plots. Closest MIREN road_bbox corner gives a lower
bound on city-center-to-source-control distance (~331 km); do NOT
claim all 28 railway stations are >=331 km away, their precise field
locations have not been joined in this audit.

Also primary botanical field survey mainly 2012–2013 vs 2022 MIREN
road plots, around 9–10 years temporally different. No direct
contemporary railway/station vs roadside natural baseline.
Original railway author study may have occasional supplemental observations;
the original sampling years do not imply all individuals observed in
a single 2013 survey. All country prediction results remain source
observational data, not independent causal adaptation.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path

BASE_SOURCE_GIT_BLOB="df8b84ba40c5eb6ee70536b188e2a4a5cfabeb9a"
CITY_PARMA_UNGEGN=(44.8018611111111,10.32735)
CITY_SOURCE_URL="https://ungegn.un.org/dashboard/cities/details?id=1538"
SOURCE_ITALY_STUDY="Adorni & Ghillani 2020; Floristic survey in railway stations of the province of Parma"
SOURCE_STUDY_PDF="https://www.museocivico.rovereto.tn.it/UploadDocs/21188_Art08_Adorni_C_agg_min.pdf"
SOURCE_RAIL_STATIONS_N=28
SOURCE_RAIL_MAIN_YEARS=(2012,2013)

def haversine_km(a,b):
    lat1,lon1=a;lat2,lon2=b
    r=math.pi/180
    delta_lat=(lat2-lat1)*r;delta_lon=(lon2-lon1)*r
    x=math.sin(delta_lat/2)**2+math.cos(lat1*r)*math.cos(lat2*r)*math.sin(delta_lon/2)**2
    return 6371.0088*2*math.asin(min(1,math.sqrt(x)))

def min_to_lonlat_bbox(city,lats,lons):
    if any(len(x)!=2 or x[0]>x[1] for x in (lats,lons)):
        raise ValueError("Original control data bounding box inconsistent")
    p=(max(lats[0],min(city[0],lats[1])),max(lons[0],min(city[1],lons[1])))
    return haversine_km(city,p),p

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+b"\0"+raw).hexdigest()

def analyze(country):
    if country["sites"].get("IT")!=SOURCE_RAIL_STATIONS_N:
        raise ValueError("Original Parma 28 real station species survey cohort changed")
    original=country["controls"]["IT_road"]
    if original.get("roads")!=3 or original.get("plots")!=60 or original.get("years")!=[2022]:
        raise ValueError("Original 2022 Italian MIREN 3-road 60-plot domestic comparator changed")
    lat,lon=original["lat"],original["lon"]
    citydist,corner=min_to_lonlat_bbox(CITY_PARMA_UNGEGN,lat,lon)
    if not 320<citydist<345:raise ValueError("Geographic source bbox comparison produced impossible Parma-to-road distance")
    r=country["comparisons"]["IT_road"]
    if abs(r["common"]["delta"]-.30809816690971126)>1e-9:
        raise ValueError("Original latest Italy source country rank comparison drifted")
    return {"version":"v25.2",
      "status":"EXTERNAL_GEO_SOURCE_PARMA_VS_DOMESTIC_ITALIAN_ROAD_CONTROL_NOT_LOCAL",
      "original_2026_Nov17_country_validation_source_commit":"a8138a50c7db3abad7ea2f717dc25f374a34e87c",
      "original_CZ_IT_HR_BA_SK_SI_IN_AU_CN_comparison_data_quality":"secondary published original survey study results",
      "original_rail_study":SOURCE_ITALY_STUDY,
      "source_original_rail_study_PDF":SOURCE_STUDY_PDF,
      "source_original_rail_station_28_from_Parma_PROVINCE":SOURCE_RAIL_STATIONS_N,
      "source_Preslia_Italy_local_station_coordinates_all_28_verified":False,
      "rail_original_main_survey_2012_2013":list(SOURCE_RAIL_MAIN_YEARS),
      "domestic_reference_MIREN_mountain_road_source":{
          "n_true_2022_road_roads":original["roads"],
          "n_true_source_road_side_plot_observations":original["plots"],
          "survey_years":original["years"],
          "latitude_min_max_N":lat,"longitude_min_max_E":lon,
          "source_2022_minus_latest_original_main_rail_year":2022-SOURCE_RAIL_MAIN_YEARS[1],
          "source_2022_minus_earliest_original_main_rail_year":2022-SOURCE_RAIL_MAIN_YEARS[0]},
      "Parma_city_UNGEGN_official_point":{
         "latitude_N":CITY_PARMA_UNGEGN[0],"longitude_E":CITY_PARMA_UNGEGN[1],
         "source_url":CITY_SOURCE_URL},
      "distance_from_city_center_to_nearest_point_inside_Italian_MIREN_source_observation_bbox_km":citydist,
      "nearest_geographic_corner_not_a_confirmed_observed_MIREN_field_plot":{
          "latitude_N":corner[0],"longitude_E":corner[1]},
      "source_holdout_IT_foreign_rail_vs_domestic_road_rank_advantage_common_taxa":
          r["common"]["delta"],
      "source_holdout_IT_genus_resampled_common_rank_interval":r["common"]["genus_interval"],
      "source_holdout_IT_rank_advantage_extended_taxon_panel":r["extended"]["delta"],
      "source_external_Italy_road_reference_2022_does_not_match_Parma_rail_2012_13_geography_or_era":True,
      "not_claim_each_of_the_28_original_Parma_railway_stations_is_331km_from_road_reference":True,
      "not_a_causal_railway_specific_flora_preference_or_mowed_plant_adaptation_test":True,
      "interpretation":"The one observational country comparator passing all eight post hoc source gate checks is NOT geographically/temporally matched to Parma provincial railway station sites. Region-specific railway habitat filtering cannot be inferred from comparing country-matched but distant different-biome and era railway vs mountain road surveys."}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",default="research/validation_20261117/country_validation.json")
    p.add_argument("--output",required=True)
    a=p.parse_args()
    raw=Path(a.source).read_bytes()
    if git_blob(raw)!=BASE_SOURCE_GIT_BLOB:
        raise ValueError("Original latest Nov17 country_validation.json Git blob changed, stop")
    result=analyze(json.loads(raw))
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("V252_ITALY_DOMESTIC_ROAD_NOT_LOCAL_SOURCE="+json.dumps(result,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
