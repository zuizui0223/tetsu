#!/usr/bin/env python3
"""v23.9 source-original spatial local nontransport plot availability for true Czech 39 rail stations.

ORIGINALS:
 - Kutlvašr et al 2024 Preslia 39 station flora + corrected S2 original
   (2020–22 station-wide 3 habitat zones, 761 unique currently recorded taxa)
 - Sap1k published Git-derived Czech official SR70 source station point
   effective 2026-08-15; exact original 39 full station names -> 35
   unambiguous 2026 station REPRESENTATIVE points (not 2020 plant GPS)
 - Chytrý et al 2026 CzechVeg-Open v1.1 source full 117739 plot
   headers (MD5 pinned), fields DEG_LAT/DEG_LON WGS84 decimal,
   PRECISION source location uncertainty meters, original survey YEAR,
   original original-publication REFERENCE_TITLE, original HABITAT/
   LOCALITY broad text labels.

This is an OPPORTUNITY AUDIT, NOT proof of railway effects.
Station whole-yard regional flora vs 1-1000 m² CzechVeg plots have
very different sampling grains. Matching a 2026 station representative
point to a surrounding 2000/2015+ natural plot does not create an
independent, equal-survey-effort plot. Text with no "transport" label is
NOT necessarily true natural non-rail habitat. There is NO direct
plant sample GPS for original station-only species.
No original source exact site GPS or plot ID is published here.

Criteria: only original source nontransport-text plots with source
WGS84 coordinates and coordinate position uncertainty<=1000 meters,
separately <=250 m. Plot age 2015–2023 or 2000–2023, and air distance
<=5km / <=10km / <=25km of exact matching SR70 station centroids.
A station is "locally supported" if >=5 different plot IDs and >=2
different original source publications; separately >=10/3 studies.
This assesses source feasibility only. For stronger study, source
CzechVeg species incidence and physical environment plus plot vs
station grain standardization must be solved.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    get_meta
)
from audit_CzechVeg_abiotic_regional_reference_opportunity_v194 import (
    number
)
from audit_Czech_2024_S2_to_SR70_original_station_points_v153 import (
    run as czech_true_station_identity
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    haversine_km
)

REQUIRED_GEOGRAPHIC_SOURCE_COLUMNS=("DEG_LAT","DEG_LON","PRECISION","RELEVE_NR")
YEAR_COHORTS=(2000,2015)
PRECISION_GATES=(250,1000)
RADII_KM=(5,10,25)
CZECH_BOUNDS=(48.0,52.0,12.0,20.0)
EXPECTED_STATION_MATCHED=35
EXPECTED_SOURCE_PLOTS=117739

def original_reference_geo(source_header,meta):
    source={};qa=Counter();quality=defaultdict(Counter)
    with open(source_header,encoding="utf-8-sig",newline="") as fh:
        rows=csv.DictReader(fh)
        if set(REQUIRED_GEOGRAPHIC_SOURCE_COLUMNS)-set(rows.fieldnames or []):
            raise ValueError("CzechVeg original v1.1 missing source original WGS84 decimal degree coords or precision")
        for r in rows:
            qa["publisher_original_total_headers"]+=1
            pid=str(r.get("RELEVE_NR") or "").strip()
            if pid not in meta:
                qa["meta_ineligible_missing_year_publication_region_area"]+=1
                continue
            m=meta[pid]
            if m["label"] is not None or m["strict"] is not None:
                qa["original_plot_has_transport_text_not_nontransport_control"]+=1
                continue
            lat=number(r.get("DEG_LAT"),-90,90)
            lon=number(r.get("DEG_LON"),-180,180)
            precision=number(r.get("PRECISION"),0,1_000_000)
            if lat is None or lon is None:
                qa["source_original_wgs84_coordinate_missing"]+=1;continue
            if not(CZECH_BOUNDS[0]<=lat<=CZECH_BOUNDS[1]
                   and CZECH_BOUNDS[2]<=lon<=CZECH_BOUNDS[3]):
                qa["original_coordinate_outside_czech_bounds"]+=1;continue
            if precision is None:
                qa["source_original_position_precision_missing"]+=1;continue
            if precision>max(PRECISION_GATES):
                qa["source_original_position_precision_above_1km"]+=1;continue
            if m["year"]>2023:
                qa["original_survey_after_field_date_2023"]+=1;continue
            record={
              "lat":lat,"lon":lon,"precision_m":precision,
              "year":m["year"],"study":m["study"],
              "source_region":m["region"],"source_area_class":m["area"],
            }
            source[pid]=record
            for cohort in YEAR_COHORTS:
                if cohort<=m["year"]<=2023:
                    qa["source_eligible_coordinates_since_"+str(cohort)]+=1
            for gate in PRECISION_GATES:
                if precision<=gate:
                    qa["source_coordinate_precision_le_"+str(gate)+"m"]+=1
    if qa["publisher_original_total_headers"]!=EXPECTED_SOURCE_PLOTS:
        raise ValueError("CzechVeg full original 117739 source header count mismatch")
    return source,dict(qa)

def spatial_opportunity(matched_stations,source_control_plots):
    # Keep derived station coordinates and exact original plot coords
    # in CI memory only; never export station- or relevé-level coordinates.
    active=[]
    for row in matched_stations:
        xy=row.get("matched_czech_2026_source_station_centroid_WGS84")
        if xy is not None:
            active.append({
               "code":row["original_locality_code"],
               "source_original_rail_stretch":row["source_section"],
               "lat":xy["lat"],"lon":xy["lon"]})
    if len(active)!=EXPECTED_STATION_MATCHED:
        raise ValueError("v15.3 true source registry exactly 35/39 stations not recovered")
    points=list(source_control_plots.values())
    bands=defaultdict(lambda:defaultdict(list))
    site_qc={}
    for src in active:
        for cohort in YEAR_COHORTS:
            for precision in PRECISION_GATES:
                for km in RADII_KM:
                    vals=[]
                    for plot in points:
                        if plot["year"]<cohort or plot["precision_m"]>precision:continue
                        # Latitude pre-screen only to speed exact original
                        # WGS84 distance, never accept a plot on this alone.
                        if abs(src["lat"]-plot["lat"])>km/110.+.01:continue
                        d=haversine_km(src["lat"],src["lon"],plot["lat"],plot["lon"])
                        if d<=km:vals.append((d,plot))
                    studies=len({p["study"] for d,p in vals})
                    # The key holds the ORIGINAL 5 rail stretches, but
                    # no exact station names or plot IDs are output.
                    bands[(cohort,precision,km)][src["source_original_rail_stretch"]].append({
                      "n_plots":len(vals),"n_distinct_original_studies":studies,
                      "min_source_distance_km":min((d for d,p in vals),default=None)})
    result={}
    for (cohort,precision,km),slices in sorted(bands.items()):
        flat=[v for arr in slices.values() for v in arr]
        if len(flat)!=EXPECTED_STATION_MATCHED:raise ValueError("Source station support frame changed")
        eligible5=[v for v in flat if v["n_plots"]>=5 and v["n_distinct_original_studies"]>=2]
        eligible10=[v for v in flat if v["n_plots"]>=10 and v["n_distinct_original_studies"]>=3]
        result[f"year_since_{cohort}_precision_{precision}m_radius_{km}km"]={
            "n_original_matched_author_stations":len(flat),
            "original_stations_ge5_local_reference_plots_from_ge2_publications":len(eligible5),
            "original_stations_ge10_local_reference_plots_from_ge3_publications":len(eligible10),
            "median_control_plots_per_matched_station":median(x["n_plots"] for x in flat),
            "median_distinct_source_publications_per_station":median(x["n_distinct_original_studies"] for x in flat),
            "n_original_5_rail_stretches_with_one_or_more_supported_stations":{
               str(stretch):sum(v["n_plots"]>=5 and v["n_distinct_original_studies"]>=2 for v in vals)
               for stretch,vals in sorted(slices.items())},
            "surveyed_2015plus_status_not_equal_to_source_rail_2020_22_taxon_record_time":True,
        }
    return result,{
       "n_original_Czech_rail_stations_matched_registry_centroid":len(active),
       "n_original_pub_rail_stretches":len({x["source_original_rail_stretch"] for x in active}),
       "source_2026_rail_centroids_NOT_botanical_subplot_locations":True,
       "plot_and_station_exact_coordinates_not_published":True}

def run(s1,s2,sr70,headers):
    stations=czech_true_station_identity(s1,s2,sr70)
    meta,qc=get_meta(headers)
    nontransport,geoqa=original_reference_geo(headers,meta)
    matched=stations["all_39_station_original_source_identity_and_candidate_provenance"]
    opportunities,quality=spatial_opportunity(matched,nontransport)
    return {"version":"v23.9",
       "status":"GENUINE_CZECH_35RAIL_STATION_VS_CZECHVEG_INDEPENDENT_LOCAL_NONTRANSPORT_PLOT_OPPORTUNITY",
       "source_Czech_botanical_39stations_DOI":"10.23855/preslia.2024.247",
       "source_current_rail_station_registry_2026_provenance":
             stations["third_party_official_rail_registry_derivative"],
       "original_Preslia_botanical_source_provenance":
             stations["publisher_S1_and_corrected_S2_authentic_provenance"],
       "source_CzechVeg_publication_DOI":"10.5281/zenodo.19231152",
       "source_CzechVeg_plot_original_MD5":qc["publisher_original_headers_md5"],
       "source_original_nontransport_plot_geographic_quality":geoqa,
       "source_exact_name_matched_geographic_station_frame":quality,
       "local_original_nontransport_publication_reference_opportunity":opportunities,
       "source_nontransport_text_unverified_true_natural":True,
       "no_individual_station_plant_or_plot_coordinates_exported":True,
       "interpretation_limits":[
           "No plant species from CzechVeg were read or interpreted as field richness, true railway site beta not evaluated",
           "SR70 official derivative is effective 2026 vs author original 2020–22 field plots; point is station representative NOT original sampled plant polygon",
           "CzechVeg field plots and whole original Preslia station species lists are NOT the same site area or botanical survey protocol",
           "No transportation word in original source HABITAT or LOCALITY cannot guarantee a plot physically outside actual railway land",
           "Exact coordinate uncertainty <=250/1000m and distance radius are necessary but insufficient to match gravel/soil or original source vegetation type",
           "CzechVeg source plots may be from old years or multiple habitat types, sampling bias persists",
           "Nonrail plot availability from 5 original railway stretches does not create five experimental study systems",
           "Potential comparison must model nonrail natural species pool and unequal effort, not subtract station 100+ species from typical plot 10-20 species"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cz_s1","cz_s2","sr70","headers","output"):
        p.add_argument("--"+k.replace("_","-"),required=True)
    a=p.parse_args()
    ans=run(a.cz_s1,a.cz_s2,a.sr70,a.headers)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(ans,ensure_ascii=False,indent=2),encoding="utf8")
    print("V239_TRUE_CZECH_STATION_NONRAIL_SPATIAL_REFERENCE="+json.dumps({
      "status":ans["status"],
      "source_exact_name_matched_geographic_station_frame":ans["source_exact_name_matched_geographic_station_frame"],
      "source_original_nontransport_plot_geographic_quality":ans["source_original_nontransport_plot_geographic_quality"],
      "local_original_nontransport_publication_reference_opportunity":ans["local_original_nontransport_publication_reference_opportunity"]
    },ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
