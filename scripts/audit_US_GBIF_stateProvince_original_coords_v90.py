#!/usr/bin/env python3
"""v9.0: independent geographic QA of original GBIF stateProvince text
against ACTUAL coordinate point in U.S. Census 2025 1:500k state polygons.

Eight PRE-FROZEN U.S. states from v8.1, BEFORE and AFTER 2022 GBIF
occurrence search original row pages, two deterministic API offsets
per state/period. Source sample is BOUNDED first <=50k result offsets,
not a random or representative sample of ALL GBIF records.
All-Plantae original source and separately U.S.-country-GRIIS-listed
plant record subset (when enough sampled source rows).

No precise occurrence coordinates or individual occurrence identifiers
persist in resulting artifact. This is pure source quality audit, NOT
verification of actual weed establishment or freight seed carriage.
"""
from __future__ import annotations
import argparse,io,json,math,zipfile
from collections import Counter
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import requests,shapefile
from shapely.geometry import Point,shape
from shapely.ops import transform as shp_transform
from shapely.validation import make_valid
from pyproj import Transformer
from analyze_griis_alien_rail_source_v25 import valid_griis
from harvest_US_state_GRIIS_GBIF_plant_source_v81 import (
    CENSUS_URL,parse_census_state_table,PILOT_FIPS,source_get,GBIF_URL)

CENSUS_SHAPE_URL="https://www2.census.gov/geo/tiger/GENZ2025/shp/cb_2025_us_state_500k.zip"
PERIODS={"before_2022":"2018,2021","after_2022":"2023,2025"}
REQUEST_OFFSETS=(0,1000)
ROWS_PER_PAGE=100
MAX_SOURCE_GBIF_POINTS=8*2*len(REQUEST_OFFSETS)*ROWS_PER_PAGE
MAX_ZIP_BYTES=12000000
MAX_SOURCE_COUNT_DRIFT=.03
GEO_BORDER_TOLERANCE_METERS=5000
MIN_GRIIS_GEOPROBED_ROWS_PER_STATE_PERIOD=10
EQUAL_DISTANCE_CRS="EPSG:5070"
SOURCE_SEED_NOT_RANDOM_ALL_RECORDS=True


def census_original_polygons(zipraw):
    """Official Census 2025 cartographic geometry in WGS84, match EXACT
    selected FIPS state source frame. Not a site vegetation polygon."""
    if not (500000<len(zipraw)<MAX_ZIP_BYTES) or not zipfile.is_zipfile(io.BytesIO(zipraw)):
        raise ValueError("Original Census 2025 STATE boundary source missing/truncated")
    with zipfile.ZipFile(io.BytesIO(zipraw)) as archive:
        names={name.lower():name for name in archive.namelist()}
        opts=[name[:-4] for name in archive.namelist()
              if name.lower().endswith("_us_state_500k.shp")]
        if len(opts)!=1:raise ValueError("No single USA Census 2025 states SHP source")
        stem=opts[0]
        def member(ext):
            path=names.get((stem+ext).lower())
            if not path:raise ValueError("Census polygon source file missing "+ext)
            return io.BytesIO(archive.read(path))
        with shapefile.Reader(shp=member(".shp"),shx=member(".shx"),dbf=member(".dbf")) as reader:
            fieldnames=[x[0] for x in reader.fields[1:]]
            if not {"STATEFP","NAME","STUSPS"}.issubset(fieldnames):
                raise ValueError("Original 2025 USA state shapefile FIPS/name fields changed")
            polys={}
            for sr in reader.iterShapeRecords():
                attrs=dict(zip(fieldnames,sr.record))
                code=str(attrs["STATEFP"])
                if code not in PILOT_FIPS:continue
                if code in polys:raise ValueError("Duplicate FIPS source state geometry")
                poly=shape(sr.shape.__geo_interface__)
                if not poly.is_valid:poly=make_valid(poly)
                if poly.is_empty or poly.area<=0:
                    raise ValueError("Source Census state polygon empty")
                polys[code]={
                    "name":attrs["NAME"],"abbr":attrs["STUSPS"],"geom":poly}
    if set(polys)!=set(PILOT_FIPS):
        raise ValueError("At least one original 8-state FIPS boundary source missing")
    return polys


def occurrence_params(name,period,offset):
    if not name or period not in PERIODS or offset<0 or offset>50000:
        raise ValueError("Predeclared source state, temporal bin or source offset invalid")
    return {
        "country":"US","stateProvince":name,
        "kingdomKey":6,"taxonRank":"SPECIES",
        "year":PERIODS[period],"hasCoordinate":"true",
        "hasGeospatialIssue":"false",
        "limit":ROWS_PER_PAGE,"offset":offset,
    }


def chosen_offsets(n):
    if type(n)!=int or n<0:raise ValueError("Malformed GBIF state source count")
    if n<ROWS_PER_PAGE:return [0] if n else []
    out=[]
    for offset in REQUEST_OFFSETS:
        value=min(offset,n-ROWS_PER_PAGE)
        if value not in out:out.append(value)
    return out


def bounded_source_points(sess,code,state_name,period):
    if code not in PILOT_FIPS:raise ValueError("Not one of the 8 frozen source states")
    first=source_get(sess,occurrence_params(state_name,period,0))
    count=first.get("count")
    if type(count)!=int or count<0 or not isinstance(first.get("results"),list):
        raise ValueError("GBIF state source count/results missing")
    offsets=chosen_offsets(count)
    page0=first["results"]
    allrows=[]
    counts=[]
    seen=set()
    for off in offsets:
        page=first if off==0 else source_get(sess,occurrence_params(state_name,period,off))
        n=page.get("count")
        if type(n)!=int or not isinstance(page.get("results"),list):
            raise ValueError("Malformed authentic GBIF original source results")
        if count>0 and abs(n-count)/count>MAX_SOURCE_COUNT_DRIFT:
            raise ValueError("GBIF state species source count drift during source probes")
        counts.append(n)
        for row in page["results"]:
            ident=str(row.get("gbifID") or row.get("key") or "")
            if not ident or ident in seen:
                raise ValueError("Duplicate or invalid GBIF original row ID across source probe offsets")
            seen.add(ident)
            year=row.get("year")
            low,high=tuple(int(s) for s in PERIODS[period].split(","))
            if type(year)!=int or not low<=year<=high:
                raise ValueError("GBIF original occurrence event YEAR inconsistent with requested period")
            allrows.append(row)
    if len(allrows)>len(offsets)*ROWS_PER_PAGE or len(allrows)>len(REQUEST_OFFSETS)*ROWS_PER_PAGE:
        raise ValueError("GBIF original source exceeded fixed page budget")
    return allrows,{
      "source_GBIF_count_in_stateProvince_query":count,
      "original_source_offsets_probed":offsets,
      "original_source_row_pages_probed":len(offsets),
      "n_original_rows_returned":len(allrows),
      "original_source_full_count_drift_fraction":
          round((max(counts)-min(counts))/max(1,count),7),
      "NOT_random_probability_sample_of_all_GBIF_state_records":True,
    }


def assess_geographic_points(rows,poly,selected_name,accepted_keys):
    """Only save aggregate QC counts; original precise GBIF coordinates
    and gbifID never included in result."""
    to_dist=Transformer.from_crs("EPSG:4326",EQUAL_DISTANCE_CRS,always_xy=True)
    metric=shp_transform(to_dist.transform,poly)
    counts=Counter()
    by_gr={"all_plantae":Counter(),"USA_country_GRIIS_listed_only":Counter()}
    publishers=Counter()
    for row in rows:
        code=str(row.get("speciesKey") or "")
        gr=("USA_country_GRIIS_listed_only" if code in accepted_keys else None)
        groups=("all_plantae",gr) if gr else ("all_plantae",)
        dataset=str(row.get("datasetKey") or "unknown")
        publishers[dataset]+=1
        longitude=row.get("decimalLongitude")
        latitude=row.get("decimalLatitude")
        if not (type(longitude) in (float,int) and type(latitude) in (float,int)
                and math.isfinite(longitude) and math.isfinite(latitude) and
                -180<=longitude<=180 and -90<=latitude<=90):
            outcome="invalid_or_missing_reported_coordinates"
        else:
            point=Point(longitude,latitude)
            if poly.covers(point):
                outcome="coordinate_strictly_within_source_state_polygon"
            elif metric.distance(shp_transform(to_dist.transform,point))<=GEO_BORDER_TOLERANCE_METERS:
                outcome="coordinate_within_5km_cartographic_boundary_tolerance"
            else:
                outcome="coordinate_more_than_5km_outside_named_state_boundary"
        for group in groups:
            by_gr[group][outcome]+=1
        reported=str(row.get("stateProvince") or "")
        if reported.strip().lower()!=selected_name.strip().lower():
            counts["returned_stateProvince_text_differs_from_requested"]+=1
        status=str(row.get("occurrenceStatus") or "UNKNOWN")
        counts["raw_occurrenceStatus_"+status]+=1
    def summary(items):
        n=sum(items.values())
        included=items.get("coordinate_strictly_within_source_state_polygon",0)
        near=items.get("coordinate_within_5km_cartographic_boundary_tolerance",0)
        far=items.get("coordinate_more_than_5km_outside_named_state_boundary",0)
        return {
            "n_raw_GBIF_source_looked_up":n,
            "source_coordinate_QA_category_counts":dict(items),
            "fraction_strictly_in_Census_state_polygon":round(included/n,6) if n else None,
            "fraction_in_or_within_5km_Census_state_polygon":round((included+near)/n,6) if n else None,
            "fraction_more_than_5km_outside_Census_state_polygon":round(far/n,6) if n else None,
            "n_independent_source_rows_not_same_as_independent_observers":n
        }
    return {
        "all_Plantae":summary(by_gr["all_plantae"]),
        "USA_country_GRIIS_listed_records":summary(by_gr["USA_country_GRIIS_listed_only"]),
        "n_source_publisher_dataset_UUIDs_in_probe":len(publishers),
        "top_source_dataset_fraction_in_probe":round(max(publishers.values())/len(rows),6) if rows else None,
        "source_state_text_and_status":dict(counts),
        "no_exact_coordinate_or_GBIF_record_ID_in_artifact":True,
    }


def run(census_text,shape_zip,griis,sess):
    states=parse_census_state_table(census_text)
    polygons=census_original_polygons(shape_zip)
    listed=valid_griis(griis).get("US")
    if listed is None or len(listed)<100:
        raise ValueError("Original country-GRIIS U.S. checklist source missing")
    out={
        "complete":True,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "source_GBIF_periods":PERIODS,
        "frozen_8_state_FIPS_source":list(PILOT_FIPS),
        "2025_Census_state_shapefile_source_url":CENSUS_SHAPE_URL,
        "original_Census_state_polygon_zip_SHA256":sha256(shape_zip).hexdigest(),
        "original_Census_state_FIPS_txt_SHA256":sha256(census_text.encode()).hexdigest(),
        "sample_design":"2 deterministic stateProvince occurrence/search offsets (0, 1000) per period; bounded ≤200 raw source records per window; NOT random nationwide",
        "quality_boundary_tolerance_meters":GEO_BORDER_TOLERANCE_METERS,
        "states":{},"errors":{},
        "NO_evidence_for_rail_seed_distribution_or_state_nativity":True}
    for code in PILOT_FIPS:
        name=states[code]["name"]
        poly=polygons[code]
        if name!=poly["name"] or states[code]["abbreviation"]!=poly["abbr"]:
            raise ValueError("Original Census STATE FIPS name/abbreviation source mismatch")
        state={"name":name,"periods":{}}
        for period in PERIODS:
            try:
                records,audit=bounded_source_points(sess,code,name,period)
                qa=assess_geographic_points(records,poly["geom"],name,listed)
                state["periods"][period]={**audit,**qa,"source_complete":True}
                print(json.dumps({
                    "FIPS":code,"period":period,
                    "n_probed_rows":audit["n_original_rows_returned"],
                    "all_Plantae_far_fraction":qa["all_Plantae"][
                        "fraction_more_than_5km_outside_Census_state_polygon"],
                    "GRIIS_listed_probed":qa["USA_country_GRIIS_listed_records"][
                        "n_raw_GBIF_source_looked_up"]},ensure_ascii=False),flush=True)
            except Exception as exc:
                out["errors"][code+"_"+period]=str(type(exc).__name__)+": "+str(exc)[:150]
                state["periods"][period]={"source_complete":False}
        out["states"][code]=state
    out["complete"]=len(out["states"])==len(PILOT_FIPS) and not out["errors"]
    allvals=[v["periods"][p] for v in out["states"].values() for p in PERIODS]
    def summarize(key):
        src=[entry[key] for entry in allvals if entry.get("source_complete")]
        total=sum(z["n_raw_GBIF_source_looked_up"] for z in src)
        correct=sum(z["source_coordinate_QA_category_counts"].get(
             "coordinate_strictly_within_source_state_polygon",0) for z in src)
        near=sum(z["source_coordinate_QA_category_counts"].get(
             "coordinate_within_5km_cartographic_boundary_tolerance",0) for z in src)
        far=sum(z["source_coordinate_QA_category_counts"].get(
             "coordinate_more_than_5km_outside_named_state_boundary",0) for z in src)
        return {
            "n_original_GBIF_occurrence_rows_probed":total,
            "n_within_or_5km_cartographic_border":correct+near,
            "n_more_than_5km_outside_named_Census_state_polygon":far,
            "far_outside_fraction_in_bounded_source_probe":round(far/total,6) if total else None,
            "source_probe_denominator_not_random_unbiased_state_estimate":True,
        }
    out["summary"]={
        "n_8states_source_complete_both_periods":sum(all(
            v["periods"][p].get("source_complete") for p in PERIODS)
            for v in out["states"].values()),
        "all_Plantae":summarize("all_Plantae"),
        "USA_country_GRIIS_listed_records":summarize("USA_country_GRIIS_listed_records")
    }
    out["limitations"]=[
        "The GBIF occurrence/search offset probe is deterministic, not a random sample of all US GBIF records or all source GRIIS taxa; records beyond the first 50,000 offsets are not sampled.",
        "Cartographic Census 1:500k state polygons simplify coastlines and may omit tiny islands; +5km metric tolerance separates obvious state misplacement from border generalization, not true coordinate error.",
        "GBIF stateProvince text may refer to other geographies or project metadata. An outside polygon does not prove an actual plant was in the other state.",
        "Country GRIIS-listed taxa are not necessarily state-invasive; original GBIF positive indexed speciesKey facets are not independently georeferenced surveys.",
        "Provider/source datasets can concentrate sampled pages; source sample should not be used directly to correct state totals without a documented inclusion design.",
        "No rail freight route, viable weed seed dispersal, mowing, wind, regenerative phenotype plasticity or genetic adaptation measured."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--griis",required=True);p.add_argument("--output",required=True)
    a=p.parse_args()
    r1=requests.get(CENSUS_URL,timeout=(20,90))
    r1.raise_for_status()
    r2=requests.get(CENSUS_SHAPE_URL,timeout=(20,120))
    r2.raise_for_status()
    aliases=json.loads(Path(a.griis).read_text(encoding="utf-8"))
    with requests.Session() as session:
        session.headers.update({"User-Agent":"tetsu-v90-GBIF-stateProvince-2025-Census-polygon-sourceQA"})
        result=run(r1.text,r2.content,aliases,session)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":result["complete"],
       "summary":result["summary"],
       "errors":result["errors"],
       "not_verified_state_species_establishment":True},ensure_ascii=False),flush=True)
    if not result["complete"]:
        raise SystemExit("Original GBIF-geographical source QC incomplete; no ecological inference")


if __name__=="__main__":
    main()
