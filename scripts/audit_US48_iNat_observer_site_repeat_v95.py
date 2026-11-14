#!/usr/bin/env python3
"""v9.5 US48 REAL original iNaturalist-GBIF occurrence rows: cross-period
observer identity and coarse geographic-site overlap.

Source: GBIF dataset UUID 50c9509d-22c7-4a22-a47d-8c48425ef4a7
for the SAME 48 original state FIPS sources verified v9.3.
Early event years 2018–2021, late 2023–2025, NO 2022.
Deterministic shallow occurrence/search source API offsets (0,1000);
not representative random sampling of all iNaturalist observations.
Never report original usernames, observer IDs, observation IDs, or
precise latitude/longitude in output or logs. No tracking individuals.

Repeat observer/coarse site in these source records is NOT repeated
standardized biological field sampling. NO actual railway freight seeds,
new plant colonization, cutting/wind filtering or evolved adaptation.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from statistics import median
import requests
from pyproj import Transformer

from audit_US48_iNaturalist_GBIF_source_dominance_v93 import (
    INAT_GBIF_DATASET_UUID,WINDOWS)
from harvest_US_state_GRIIS_GBIF_plant_source_v81 import (
    GBIF_URL,CENSUS_URL,parse_census_state_table,source_get)
from analyze_griis_alien_rail_source_v25 import valid_griis

PAGES=(0,1000)
PAGE_SIZE=100
CELL_KM=(10,50)
SRC_CRS="EPSG:4326"
GRID_CRS="EPSG:5070"
MAX_INDEX_DRIFT=.025
FIXED_48_STATES=48
MIN_KNOWN_OBSERVER_RECORDS=20
MIN_COORDINATE_VALID_RECORDS=20
EXPECTED_BINS={"before_2022":"2018,2021","after_2022":"2023,2025"}
UNINFORMATIVE_OBSERVER_VALUES=frozenset((
    "unknown","anonymous","not available","n/a","na","none",
    "various","multiple observers","citizen science","inat",
    "unrecorded","observer"))
NOT_RANDOM_SOURCE_SAMPLE=True


def normalize_observer(record):
    """Observer labels are used in volatile memory ONLY; never persist
    their strings, hashes, user URLs or exact observer-specific counts.
    Prioritize stable recordedByID, fall back to Darwin Core recordedBy.
    """
    def parse(value):
        if isinstance(value,list):
            values=[str(x).strip() for x in value if x is not None]
            value=values[0] if values else ""
        if not isinstance(value,str):
            return None
        key=" ".join(value.casefold().strip().split())
        if not key or key in UNINFORMATIVE_OBSERVER_VALUES:
            return None
        if len(key)>300:
            return None
        return key
    return parse(record.get("recordedByID")) or parse(record.get("recordedBy"))


def source_query(name,window,offset):
    if (window not in EXPECTED_BINS or type(offset)!=int or offset not in PAGES
        or not name):
        raise ValueError("Only independently frozen 2018–21 /2023–25 U.S. source site probe offsets permitted")
    return {
        "country":"US","stateProvince":name,"datasetKey":INAT_GBIF_DATASET_UUID,
        "kingdomKey":6,"taxonRank":"SPECIES",
        "year":EXPECTED_BINS[window],"hasCoordinate":"true",
        "hasGeospatialIssue":"false","limit":PAGE_SIZE,"offset":offset
    }


def fetch_original_rows(session,name,window,expected_indexed):
    if type(expected_indexed)!=int or expected_indexed<=0:
        raise ValueError("Missing original GBIF scoped iNaturalist index denominator")
    queried=[]
    seen=set()
    pages=[]
    for offset in PAGES:
        if offset>=expected_indexed and offset!=0:continue
        data=source_get(session,source_query(name,window,offset))
        total=data.get("count")
        records=data.get("results")
        if type(total)!=int or total<0 or not isinstance(records,list) or len(records)>PAGE_SIZE:
            raise ValueError("GBIF original iNaturalist occurrence/search response malformed")
        if abs(total-expected_indexed)/expected_indexed>MAX_INDEX_DRIFT:
            raise ValueError("Original iNaturalist dataset-filtered GBIF source count drift >2.5% since v9.3")
        if not records:
            raise ValueError("Source iNaturalist occurrence rows empty at requested recorded-data offset")
        for r in records:
            rid=str(r.get("gbifID") or r.get("key") or "").strip()
            if not rid or rid in seen:
                raise ValueError("Duplicate, absent original occurrence ID or overlapped source pages")
            seen.add(rid)
            yr=r.get("year")
            low,high=(int(s) for s in EXPECTED_BINS[window].split(","))
            if type(yr)!=int or not low<=yr<=high:
                raise ValueError("Original iNaturalist event YEAR outside frozen GBIF temporal window")
            if str(r.get("datasetKey","")).lower()!=INAT_GBIF_DATASET_UUID:
                raise ValueError("Actual original GBIF occurrence row not from scoped iNaturalist UUID")
            queried.append(r)
        pages.append({"requested_original_GBIF_offset":offset,"n_actual_original_rows":len(records)})
    return queried,{
        "source_indexed_total_records_from_v93":expected_indexed,
        "n_original_source_sampled_records":len(queried),
        "source_GBIF_offset_pages":pages,
        "nonrandom_ordered_search_API_probe_not_randomized_sampling":True
    }


def cell_of_record(row,transformer,km):
    if km not in CELL_KM:raise ValueError("Only frozen coarse 10/50 km source grid")
    lon=row.get("decimalLongitude")
    lat=row.get("decimalLatitude")
    if (type(lon) not in (int,float) or type(lat) not in (int,float)
        or not math.isfinite(lon) or not math.isfinite(lat)
        or not -180<=lon<=180 or not -90<=lat<=90):
        return None
    x,y=transformer.transform(lon,lat)
    if not(math.isfinite(x) and math.isfinite(y)):return None
    scale=km*1000
    return (math.floor(x/scale),math.floor(y/scale))


def summarize_source_window(rows,accepted):
    """Keep identifiable strings + exact coordinates ONLY in memory;
    archive deidentified AGGREGATES by state and year-window.
    """
    t=Transformer.from_crs(SRC_CRS,GRID_CRS,always_xy=True)
    observers=set()
    sites={km:set() for km in CELL_KM}
    observer_sites={km:set() for km in CELL_KM}
    event_summary=Counter()
    staged=[]
    for row in rows:
        obs=normalize_observer(row)
        geo={km:cell_of_record(row,t,km) for km in CELL_KM}
        is_griis=str(row.get("speciesKey") or "") in accepted
        if is_griis:event_summary["record_US_country_GRIIS_listed"]+=1
        if obs:
            observers.add(obs)
            event_summary["records_with_observer_id_or_name"]+=1
        else:
            event_summary["records_without_usable_observer_identity"]+=1
        for km in CELL_KM:
            tile=geo[km]
            if tile is None:continue
            sites[km].add(tile)
            if obs:observer_sites[km].add((obs,tile))
        if geo[CELL_KM[0]] is not None:
            event_summary["records_with_usable_coordinates"]+=1
        staged.append((obs,geo,is_griis))
    return {
        "n_original_rows":len(rows),
        "n_observer_labels_in_sample_not_user_population":len(observers),
        "n_source_country_GRIIS_listed_plant_records_in_sample":
            event_summary["record_US_country_GRIIS_listed"],
        "n_records_with_usable_observer_identity":
            event_summary["records_with_observer_id_or_name"],
        "n_records_with_valid_original_coordinates":
            event_summary["records_with_usable_coordinates"],
        "n_coarse_10km_source_cells":len(sites[10]),
        "n_coarse_50km_source_cells":len(sites[50]),
        "original_source_observer_labels_and_sites_not_exported":True,
        "_staged":staged,"_observers":observers,
        "_sites":sites,"_observer_sites":observer_sites
    }


def temporal_overlap(early,late):
    """Counts for origin observer/source label intersection and same
    coarse tile; no individual observer, geocell or personal metadata
    are returned to user artifacts.
    """
    es=early["_observers"]
    ls=late["_observers"]
    result={
       "n_early_distinct_observer_labels_sampled":len(es),
       "n_late_distinct_observer_labels_sampled":len(ls),
       "n_observer_labels_in_BOTH_bounded_samples":len(es&ls),
       "n_late_source_rows_from_observer_seen_in_early_sample":sum(
           observer in es for observer,_,_ in late["_staged"] if observer),
       "n_late_source_rows_with_observer_identity":sum(
           bool(observer) for observer,_,_ in late["_staged"]),
       "observer_repeat_source_supported":(
           early["n_records_with_usable_observer_identity"]>=MIN_KNOWN_OBSERVER_RECORDS and
           late["n_records_with_usable_observer_identity"]>=MIN_KNOWN_OBSERVER_RECORDS),
       "coarse_site_overlap":{},
       "repeated_observer_coarse_location_overlap":{},
       "source_records_NOT_repeated_standardized_plot_surveys":True,
    }
    for km in CELL_KM:
        early_tiles=early["_sites"][km]
        late_tiles=late["_sites"][km]
        common=early_tiles&late_tiles
        repeat=early["_observer_sites"][km]&late["_observer_sites"][km]
        late_valid=[(obs,geos[km],isg) for obs,geos,isg in late["_staged"]
                    if geos[km] is not None]
        late_obs=[(obs,tile,isg) for obs,tile,isg in late_valid if obs]
        out={
          "n_early_distinct_geographical_cells_in_sample":len(early_tiles),
          "n_late_distinct_geographical_cells_in_sample":len(late_tiles),
          "n_geographical_cells_repeated_in_both_sample_year_bins":len(common),
          "n_late_records_in_early_sampled_geographical_cells":
              sum(tile in early_tiles for _,tile,_ in late_valid),
          "n_late_records_with_original_valid_coordinates":len(late_valid),
          "n_late_GRIIS_records_in_repeated_geographical_cells":
              sum(tile in early_tiles for _,tile,isg in late_valid if isg),
          "n_distinct_early_and_late_shared_observer_and_geographical_cell":
              len(repeat),
          "n_late_records_with_same_observer_same_coarse_site_as_early":
              sum((obs,tile) in early["_observer_sites"][km]
                  for obs,tile,_ in late_obs),
          "n_late_records_with_both_observer_and_coarse_location":len(late_obs),
          "NOT_plant_individual_or_true_revisit":True
        }
        result["coarse_site_overlap"][str(km)+"km"]=out
    return result


def state_record_summary(early,late,matched):
    """Strip all internal raw observer/site labels and records."""
    e={key:value for key,value in early.items() if not key.startswith("_")}
    l={key:value for key,value in late.items() if not key.startswith("_")}
    obs_count=matched["n_late_source_rows_with_observer_identity"]
    repeated=matched["n_late_source_rows_from_observer_seen_in_early_sample"]
    result={
       "source_pre_2018_21":e,
       "source_post_2023_25":l,
       "overlap":matched,
       "fraction_late_original_rows_with_repeated_prior_observer_label":
            round(repeated/obs_count,6) if obs_count else None,
       "observer_and_10km_50km_sampled_overlap_are_NOT_local_plant_establishment":True
    }
    if any(str(u) in json.dumps(result) for u in ("inaturalist.org/people/","inaturalist.org/observations/")):
        raise ValueError("Confidential source individual observer/record leaked into output")
    return result


def run(temporal,inat,accepted,census,sess):
    if temporal.get("complete") is not True or inat.get("complete") is not True:
        raise ValueError("Independent GBIF 48-state source counts incomplete")
    if (inat.get("source_original_GBIF_qualified_states")!=FIXED_48_STATES or
        inat.get("iNaturalist_original_GBIF_dataset_UUID")!=INAT_GBIF_DATASET_UUID or
        temporal.get("source_periods")!=EXPECTED_BINS or len(census)!=51 or len(accepted)<100):
        raise ValueError("Original iNaturalist UUID, source years, country taxon or Census geography changed")
    keys=sorted(inat["states"])
    if len(keys)!=FIXED_48_STATES:
        raise ValueError("Actual iNaturalist full 48-state source records missing")
    result={
        "complete":True,
        "UTC":datetime.now(timezone.utc).isoformat(),
        "study":"Bounded original GBIF iNaturalist same-observer and 10/50km source site recurrent recording audit",
        "predefined_48_original_US_Census_state_codes":keys,
        "country_GRIIS_US_listed_plant_keys_not_state_intro_count":len(accepted),
        "original_record_source":"GBIF occurrence/search scoped iNaturalist Research-grade dataset UUID",
        "source_UUID":INAT_GBIF_DATASET_UUID,
        "record_years":EXPECTED_BINS,
        "bounded_deterministic_search_offsets":PAGES,
        "n_requested_source_pages_per_state_period":len(PAGES),
        "source_page_size":PAGE_SIZE,
        "nonrandom_source_probe_DO_NOT_EXTRAPOLATE_TO_ALL_INAT_USERS":True,
        "states":{},"errors":{},
        "no_raw_observer_names_ids_record_ids_or_latlon_saved":True,
        "NO_actual_train_seed_carrying_or_rail_mowing_gust_adaptation":True
    }
    for code in keys:
        item=inat["states"][code]
        if item.get("complete") is not True or census[code]["name"]!=item["name"]:
            raise ValueError("Original Census FIPS state name and actual 48 iNaturalist source mismatch")
        windows={}
        try:
            for win in WINDOWS:
                total=item["source_result"]["two_GBIF_year_source_periods"][win][
                    "iNaturalist_original_all_Plantae_record_source_count"]
                observed,meta=fetch_original_rows(sess,item["name"],win,total)
                row_summary=summarize_source_window(observed,accepted)
                windows[win]={"data":row_summary,"page_audit":meta}
            match=temporal_overlap(windows["before_2022"]["data"],windows["after_2022"]["data"])
            clean=state_record_summary(windows["before_2022"]["data"],
                                      windows["after_2022"]["data"],match)
            clean["source_page_details"]={win:windows[win]["page_audit"] for win in WINDOWS}
            result["states"][code]=clean
            print(json.dumps({
                "FIPS":code,
                "n_sampled_rows_before":clean["source_pre_2018_21"]["n_original_rows"],
                "n_sampled_rows_after":clean["source_post_2023_25"]["n_original_rows"],
                "observer_labels_repeat":match["n_observer_labels_in_BOTH_bounded_samples"],
                "late_records_repeated_observer":match[
                    "n_late_source_rows_from_observer_seen_in_early_sample"],
                "late_10km_same_observer_site":match[
                    "coarse_site_overlap"]["10km"]["n_late_records_with_same_observer_same_coarse_site_as_early"],
                "n_source_errors":len(result["errors"])},ensure_ascii=False),flush=True)
        except Exception as exc:
            result["errors"][code]=type(exc).__name__+": "+str(exc)[:175]
            result["states"][code]={"source_incomplete":True}
            print(json.dumps({"FIPS":code,"source_error":result["errors"][code],
                              "n_source_errors":len(result["errors"])}),flush=True)
    success=[r for r in result["states"].values() if not r.get("source_incomplete")]
    result["complete"]=len(success)==FIXED_48_STATES and not result["errors"]
    def med(values):
        xs=[v for v in values if v is not None]
        return round(median(xs),6) if xs else None
    result["summary"]={
        "n_Census_state_sources_completed":len(success),
        "n_actual_GBIF_iNaturalist_original_sampled_rows_before":
            sum(v["source_pre_2018_21"]["n_original_rows"] for v in success),
        "n_actual_GBIF_iNaturalist_original_sampled_rows_after":
            sum(v["source_post_2023_25"]["n_original_rows"] for v in success),
        "median_early_sampled_observer_label_count":
            med([v["source_pre_2018_21"]["n_observer_labels_in_sample_not_user_population"] for v in success]),
        "median_late_sampled_observer_label_count":
            med([v["source_post_2023_25"]["n_observer_labels_in_sample_not_user_population"] for v in success]),
        "n_states_with_at_least_one_same_observer_label_in_both_bounded_samples":
            sum(v["overlap"]["n_observer_labels_in_BOTH_bounded_samples"]>0 for v in success),
        "median_fraction_late_rows_from_early_sampled_observer_label":
            med([v["fraction_late_original_rows_with_repeated_prior_observer_label"] for v in success]),
        "n_states_at_least_one_repeated_10km_tile":
            sum(v["overlap"]["coarse_site_overlap"]["10km"][
                "n_geographical_cells_repeated_in_both_sample_year_bins"]>0 for v in success),
        "n_states_at_least_one_repeated_observer_and_10km_tile":
            sum(v["overlap"]["coarse_site_overlap"]["10km"][
                "n_distinct_early_and_late_shared_observer_and_geographical_cell"]>0 for v in success),
        "median_fraction_late_rows_in_repeated_50km_source_tile":
            med([q["overlap"]["coarse_site_overlap"]["50km"][
                "n_late_records_in_early_sampled_geographical_cells"]/
                 q["overlap"]["coarse_site_overlap"]["50km"][
                "n_late_records_with_original_valid_coordinates"]
                 if q["overlap"]["coarse_site_overlap"]["50km"][
                    "n_late_records_with_original_valid_coordinates"] else None
                 for q in success]),
        "n_states_with_unknown_source_observer_identity_in_any_sample_period":
            sum(v["source_pre_2018_21"]["n_records_with_usable_observer_identity"]<
                v["source_pre_2018_21"]["n_original_rows"] or
                v["source_post_2023_25"]["n_records_with_usable_observer_identity"]<
                v["source_post_2023_25"]["n_original_rows"]
                for v in success),
    }
    result["limits"]=[
        "Exactly two source API offset pages (0 and 1000) per state-period: deterministic nonrandom search order, NOT probability sample of iNaturalist observers, locations or all GBIF observations.",
        "The source observer string may be generic, absent or modified; recordedByID and recordedBy are not independently validated same physical humans.",
        "Same observer and 10/50km coarse geographic cell sampled in both bins are NOT same survey plots, systematic revisit, site occupancy, botanical species detection or genetic repeat measures.",
        "The 48 source states repeatedly share some observers and cities, so state medians are not independent randomized biological replicates.",
        "GBIF 2018–21 and 2023–25 event source years contain different numbers of years and different iNaturalist app/reporting intensity, so missing sample observers are not absent users.",
        "GRIIS country-listed alien taxa are not confirmed alien at state level, nor does a late observed species imply arrival, colonization or seed movement.",
        "FAF freight economic OD, actual cargo-carried viable seed flux, 2026 rail/road tracks, rail mowing/herbicide, slipstream wind, plasticity, adaptation and local evolutionary genetics are NOT contained.",
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    for key in ("griis","temporal85","inat93","output"):
        p.add_argument("--"+key,required=True)
    a=p.parse_args()
    read=lambda k:json.loads(Path(getattr(a,k)).read_text(encoding="utf-8"))
    census_r=requests.get(CENSUS_URL,timeout=(20,90),
        headers={"User-Agent":"tetsu-48state-INAT-observer-site-original-GBIF-source-v95"})
    census_r.raise_for_status()
    census=parse_census_state_table(census_r.text)
    listed=valid_griis(read("griis")).get("US")
    if not listed:raise ValueError("Original U.S. country GRIIS source unavailable")
    with requests.Session() as session:
        session.headers.update({"User-Agent":"tetsu-v95-GBIF-inat-state-observer-site-source-QA"})
        out=run(read("temporal85"),read("inat93"),listed,census,session)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":out["complete"],"summary":out["summary"],
         "errors":out["errors"],"no_individual_observer_or_location_records_saved":True},
         ensure_ascii=False),flush=True)
    if not out["complete"]:
        raise SystemExit("Incomplete original iNaturalist observed-source sample; no inference")

if __name__=="__main__":
    main()
