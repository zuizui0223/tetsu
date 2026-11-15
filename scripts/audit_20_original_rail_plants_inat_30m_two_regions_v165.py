#!/usr/bin/env python3
"""v16.5 REAL iNaturalist 2018–2025 GEOLOCATED PHOTO feasibility for 20
author-original railway flora candidates preselected in SUCCESS v16.4.

This is the FIRST SOURCE-REGION wave of a potentially global sample.
Actual geometric 0–30m nearest REAL railway POLYLINE (not nearest OSM
node) computed from original Geofabrik Slovenia-2024 / Czech-2023 PBFs.
For each of 20 originally ranked botanical species, audit complete
iNaturalist photographed/research-grade source occurrences separately in
two locked COUNTRY BOUNDING BOXES. Also query GLOBAL total photo counts
as a context-only upper-bound; no global railway distance claim.
Czech and Slovenian bbox photo records are not an exhaustive world census.

Must reject public-coordinate obscuration/private/captive records and missing
or >30m photo georeferencing uncertainty, query aliases and API truncation.
A public point within 30m is NOT a fully accurate mapped *rail habitat*
observation; 'guaranteed' 30m only if distance+error_radius<=30m.
Global totals cannot establish global near-track counts. Photographed
plants need human/AI plant-size/flower-visibility audit before phenology or
height; observation date != actual first flowering day. No downloaded
photos or raw observation IDs, usernames or original coordinates exposed.
Does not identify railway preference vs nonrail regional plant source pools,
rail freight origin-destination, train wind or selection.
"""
from __future__ import annotations
import argparse,datetime as dt,hashlib,json,math,time
from collections import Counter,defaultdict
from pathlib import Path

from select_two_original_rail_20_species_for_inat_v164 import SPECIES_RE

API="https://api.inaturalist.org/v1/observations"
DATES=("2018-01-01","2025-12-31")
STUDY_REGIONS={
 "SI":{"bbox":(45.38,13.32,46.91,16.66),
       "original_OSM_PBF_SHA256":"b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39",
       "original_OSM_PBF_BYTES":280842431,
       "utm_epsg":32633,
       "OSM_year":2024},
 "CZ":{"bbox":(48.45,12.04,51.11,18.95),
       "original_OSM_PBF_SHA256":"26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747",
       "original_OSM_PBF_BYTES":816159635,
       "utm_epsg":32633,
       "OSM_year":2023}
}
SOURCE_CANDIDATES="v16.4_AUTHOR_ORIGINAL_REAL_RAIL_FLORA_20_SPECIES_FIRST_PHOTO_FEASIBILITY"
PHOTO_LICENSE_ALLOW={"cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa"}
# Narrow synonym from earlier primary-source botanical taxonomic audit in
# tetsu, never merge other genus species or treat aliases as absences.
CURATED_SOURCE_ACCEPTED_ALIASES={"Conyza canadensis":{"Erigeron canadensis"}}
RAIL_TYPES={"rail","narrow_gauge"}
ACCURACY_MAX_M=30.
RAIL_DISTANCE_M=30.
DAYS_FIRST=DATES[0];DAYS_LAST=DATES[1]
MAX_COMPLETE_PER_SOURCE_BBOX_TAXON=8000
PER_PAGE=200
API_DELAY_SECONDS=1.10

def validate_source_candidates(report):
    if report.get("version")!=SOURCE_CANDIDATES or report.get("complete") is not True:
        raise ValueError("Cannot use an unverified original 20 railway flora species study")
    rows=report.get("original_selected_20_species") or []
    if len(rows)!=20 or report.get("n_distinct_original_species")!=20:
        raise ValueError("Original source photo candidate list not exactly 20")
    names=[x["scientific_name_original_exact_binomial"] for x in rows]
    if len(set(names))!=20 or not all(SPECIES_RE.fullmatch(name) for name in names):
        raise ValueError("Selected taxa not 20 unique original scientifically valid binomials")
    if not all(x.get("NOT_railway_affinity_vs_nonnear_railway") is True and
               x.get("selection_not_based_on_iNaturalist_photos") is True
               for x in rows):
        raise ValueError("Original rail recorded taxa must be selected before iNaturalist data, not as proven rail affinity")
    if set(x["original_source_selected_from"] for x in rows)!={"SI","CZ"}:
        raise ValueError("Two authentic author original country sources required")
    return rows

def original_pbf_authentic(path,study):
    original=STUDY_REGIONS[study]
    if Path(path).stat().st_size!=original["original_OSM_PBF_BYTES"]:
        raise ValueError("Original source national Geofabrik PBF size changed")
    digest=hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda:fh.read(8*1024*1024),b""):digest.update(block)
    if digest.hexdigest()!=original["original_OSM_PBF_SHA256"]:
        raise ValueError("Original source national Geofabrik PBF SHA256 mismatch")
    return digest.hexdigest()

class RailwayPolylineSource:
    def __init__(self,source_path,study):
        # Official derived historical OSM railway rights-of-way are point-to-
        # LINE nearest-distance exposure proxies. No freight rail operation,
        # no track access or train-induced local airflow measured.
        import osmium
        from shapely.geometry import LineString
        from shapely.strtree import STRtree
        from pyproj import Transformer
        orig=STUDY_REGIONS[study]
        self.transformer=Transformer.from_crs(4326,orig["utm_epsg"],always_xy=True)
        projection=self.transformer
        source_bbox=orig["bbox"]
        self.source_ways_by_original_railway_tag=Counter()
        lines=[]
        class Handler(osmium.SimpleHandler):
            def way(that,w):
                kind=w.tags.get("railway")
                if kind not in RAIL_TYPES:return
                that_parent=self
                that_parent.source_ways_by_original_railway_tag[kind]+=1
                # No fake junctions added across separated OSM node IDs;
                # for nearest spatial location, individual way geometry
                # and its actual edges suffice.
                raw=[]
                for n in w.nodes:
                    if n.location.valid():
                        a=float(n.lon);b=float(n.lat)
                        if (math.isfinite(a) and math.isfinite(b) and
                                -180<=a<=180 and -90<=b<=90):
                            raw.append((a,b))
                if len(raw)<2:return
                x,y=projection.transform([x for x,_ in raw],[y for _,y in raw])
                if any(not(math.isfinite(a) and math.isfinite(b)) for a,b in zip(x,y)):
                    return
                line=LineString(zip(x,y))
                if line.length>0:lines.append(line)
        Handler().apply_file(str(source_path),locations=True,idx="flex_mem")
        if len(lines)<100 or not self.source_ways_by_original_railway_tag:
            raise ValueError("Original country OSM source has no authentic usable railway LINE geometries")
        self.lines=lines
        self.tree=STRtree(lines)
        self.study=study
        self.n_lines=len(lines)

    def distance_and_cell(self,lon,lat):
        from shapely.geometry import Point
        x,y=self.transformer.transform(lon,lat)
        pt=Point(x,y)
        i=int(self.tree.nearest(pt))
        d=float(pt.distance(self.lines[i]))
        if not math.isfinite(d) or d<0:
            raise ValueError("Nearest original physical railway polyline geometry distance invalid")
        return d,(int(math.floor(x/100.)),int(math.floor(y/100.)))

class ConservativeINat:
    def __init__(self,delay=API_DELAY_SECONDS):
        import requests
        self.session=requests.Session()
        self.session.headers.update({"User-Agent":"tetsu-original-botanist-citizen-photo-near-rail-feasibility/16.5 (https://github.com/zuizui0223/tetsu)"})
        self.last=0.
        self.delay=delay
        self.call_count=0

    def get(self,params):
        for attempt in range(5):
            wait=self.delay-(time.monotonic()-self.last)
            if wait>0:time.sleep(wait)
            self.last=time.monotonic()
            self.call_count+=1
            try:r=self.session.get(API,params=params,timeout=(12,95))
            except Exception:
                if attempt==4:raise
                time.sleep(3*(attempt+1))
                continue
            if r.status_code in (429,500,502,503,504):
                time.sleep(4*(attempt+1))
                continue
            r.raise_for_status()
            data=r.json()
            if "total_results" not in data or not isinstance(data.get("results"),list):
                raise ValueError("iNaturalist API changed documented observations paging schema")
            return data
        raise RuntimeError("iNaturalist API retry exhausted; never count missing API response as 0 observations")

def fixed_query(taxon,bbox=None,per_page=PER_PAGE):
    q={
        "taxon_name":taxon,"quality_grade":"research","photos":"true",
        "d1":DATES[0],"d2":DATES[1],"per_page":per_page,
        "order_by":"id","order":"asc"
    }
    if bbox is not None:
        q.update(zip(("swlat","swlng","nelat","nelng"),bbox))
    return q

def load_complete_photo_observations(client,taxon,bbox):
    p=fixed_query(taxon,bbox)
    first=client.get({**p,"page":1})
    total=int(first["total_results"])
    if total<0:raise ValueError("iNaturalist total results cannot be negative")
    if total>MAX_COMPLETE_PER_SOURCE_BBOX_TAXON:
        return {"total":total,"complete":False,
                "reason":"API_SOURCE_BBOX_EXCEEDS_8000_medium_batch_contract_not_counted_as_zero"},None
    rows=list(first["results"])
    pages=math.ceil(total/PER_PAGE)
    if pages>50:
        raise ValueError("iNaturalist documented API 10k pagination window exceeded")
    for page in range(2,pages+1):
        reply=client.get({**p,"page":page})
        if int(reply["total_results"])!=total:
            raise ValueError("Observation API country species count drift during pagination")
        rows.extend(reply["results"])
    ids=[r.get("id") for r in rows]
    if (len(ids)!=total or len(set(ids))!=total or
        any(not isinstance(i,int) for i in ids)):
        raise ValueError("Incomplete public iNaturalist original observation census; no partial 30m count")
    return {"total":total,"complete":True,"n_source_API_records":len(rows),"n_API_pages":max(pages,1)},rows

def exact_taxon_or_descendant(name,original):
    if not isinstance(name,str):return False
    equivalents={original}|CURATED_SOURCE_ACCEPTED_ALIASES.get(original,set())
    return any(
        name==target or (
            name.startswith(target+" ") and
            name[len(target)+1:].startswith(("subsp.","var.","subsp ","var ","ssp.","f.","form "))
        ) for target in equivalents
    )

def observation_quality(ob,source_species):
    if ob.get("obscured") or ob.get("geoprivacy") in ("obscured","private") or ob.get("taxon_geoprivacy") in ("obscured","private"):
        return None,"obscured_or_private"
    if ob.get("captive") or ob.get("captive") is True:
        return None,"captive_or_cultivated"
    original_name=(ob.get("taxon") or {}).get("name")
    if not exact_taxon_or_descendant(original_name,source_species):
        return None,"accepted_taxon_mismatch_unresolved_synonym_or_complex"
    a=ob.get("positional_accuracy")
    if isinstance(a,bool) or not isinstance(a,(int,float)) or not math.isfinite(a) or not 0<a<=ACCURACY_MAX_M:
        return None,"missing_or_above_30m_public_coordinate_accuracy"
    if not ob.get("photos"):
        return None,"no_photos"
    if not ob.get("observed_on"):
        return None,"date_missing"
    coords=(ob.get("geojson") or {}).get("coordinates")
    if (not isinstance(coords,list) or len(coords)!=2 or
        any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in coords)):
        return None,"missing_or_invalid_public_coordinates"
    lng,lat=coords
    if not(-180<=lng<=180 and -90<=lat<=90):
        return None,"invalid_lat_lon_public"
    return (float(lng),float(lat),float(a)),"pass"

def summarize_species_geo(source_records,original_species,rail):
    if source_records is None:
        return {"completed_original_source_photo_census":False,
                "n_point_near_30m_strict_accuracy":None,
                "n_location_certain_within_30m":None}
    bins=Counter();rejection=Counter();observers=defaultdict(set);cells=defaultdict(set)
    with_licensed=Counter();months=Counter();flower=Counter()
    for ob in source_records:
        valid,why=observation_quality(ob,original_species)
        rejection[why]+=1
        if valid is None:continue
        lon,lat,accuracy=valid
        distance,cell=rail.distance_and_cell(lon,lat)
        category=("point_within_30m" if distance<=30 else
                  "point_30to100m" if distance<=100 else
                  "point_100to1000m" if distance<=1000 else "point_beyond_1000m")
        bins[category]+=1
        observer=(ob.get("user") or {}).get("id")
        if observer is not None:observers[category].add(str(observer))
        cells[category].add(cell)
        if any(str(p.get("license_code") or "").lower() in PHOTO_LICENSE_ALLOW for p in ob.get("photos",[])):
            with_licensed[category]+=1
        if category=="point_within_30m":
            if distance+accuracy<=30:bins["certain_within_30m_accuracy_circle"]+=1
            # presence of flower in the photo is OPTIONAL phenotype,
            # absence of annotation is NOT evidence of biological nonflowering.
            ann={a.get("controlled_value_id") for a in ob.get("annotations",[])
                 if a.get("controlled_attribute_id")==12}
            if 13 in ann:flower["flower_visible_annotation"]=flower["flower_visible_annotation"]+1
            else:flower["no_verified_flower_visible_annotation"]=flower["no_verified_flower_visible_annotation"]+1
            mo=str(ob["observed_on"])[5:7]
            if len(mo)==2 and mo.isdigit():months[mo]+=1
    near=bins["point_within_30m"]
    return {
      "completed_original_source_photo_census":True,
      "n_source_API_photo_observations_returned":len(source_records),
      "n_after_strict_accuracy_and_taxon_public_geoprivacy_gate":rejection["pass"],
      "n_point_near_30m_strict_accuracy":near,
      "n_location_certain_within_30m":bins["certain_within_30m_accuracy_circle"],
      "n_point_near_30m_photo_CC_licensed":with_licensed["point_within_30m"],
      "n_distinct_observers_public_within_30m":len(observers["point_within_30m"]),
      "n_distinct_100m_public_cells_near_30m":len(cells["point_within_30m"]),
      "n_point_30to100m_strict_accuracy":bins["point_30to100m"],
      "n_point_100to1000m_strict_accuracy":bins["point_100to1000m"],
      "n_point_beyond_1000m_strict_accuracy":bins["point_beyond_1000m"],
      "near_30m_month_of_observation_counts_NOT_flowering_onset":dict(months),
      "near_30m_explicit_flower_visibility_annotation_counts_NOT_full_phenophase":dict(flower),
      "rejected_original_photo_quality_reasons":dict(rejection),
      "n_point_near_30m_with_photo_SHOWING_whole_plant_HEIGHT_verified":0,
      "whole_plant_height_photo_validation_not_performed":True,
      "not_verified_rail_right_of_way_or_train_wind":True
    }

def run(source_manifest,pbf_by_study,client=None):
    selected=validate_source_candidates(source_manifest)
    client=client or ConservativeINat()
    spatial_sources={}
    result={}
    for study in ("SI","CZ"):
        pbf=pbf_by_study[study]
        sha=original_pbf_authentic(pbf,study)
        rail=RailwayPolylineSource(pbf,study)
        spatial_sources[study]={
            "original_Geofabrik_PBF_SHA256":sha,
            "original_railway_way_type_counts":dict(rail.source_ways_by_original_railway_tag),
            "original_OSM_railways_lines_indexed":rail.n_lines,
            "original_source_photo_geographic_bbox_NOT_country_boundary":STUDY_REGIONS[study]["bbox"],
            "UTM_EPSG":STUDY_REGIONS[study]["utm_epsg"],
            "NOT_rail_waypoint_node_only_distance":True
        }
        species=[]
        for row in selected:
            name=row["scientific_name_original_exact_binomial"]
            meta,observations=load_complete_photo_observations(client,name,STUDY_REGIONS[study]["bbox"])
            profile=summarize_species_geo(observations,name,rail)
            species.append({
              "taxon":name,"original_source_selected_from":row["original_source_selected_from"],
              "source_API_bbox_photo_total":meta["total"],
              "source_query_complete":meta["complete"],
              "source_unenumerated_reason":meta.get("reason"),
              "source_original_photo_quality_and_geo_rail_polyline_30m":profile
            })
            print(json.dumps({"source_region":study,"taxon":name,
              "API_total":meta["total"],"complete":meta["complete"],
              "n_quality_gated_point_near_rail_30m":profile["n_point_near_30m_strict_accuracy"]},ensure_ascii=False),flush=True)
        result[study]=species
        del rail
    global_totals=[]
    for row in selected:
        name=row["scientific_name_original_exact_binomial"]
        data=client.get({**fixed_query(name,None,per_page=1),"page":1})
        global_totals.append({"taxon":name,"global_iNaturalist_research_photo_total_2018_2025":int(data["total_results"]),
                             "global_original_observation_count_NOT_global_rail_near_30m":True})
    output={
       "complete":True,
       "version":"v16.5_iNaturalist_first_TWO_SOURCE_COUNTRY_BBOX_real_OSM_rail_30m_feasibility",
       "source_API":API,"iNaturalist_API_policy_approximately_one_request_per_second":True,
       "source_2018_2025_dates":DATES,
       "source_candidate_20_original_taxa_selection_manifest":source_manifest["version"],
       "source_original_field_studies":["10.1007/s44473-026-00155-z","10.23855/preslia.2024.247"],
       "data_scope":"TWO source-study-country BOUNDING BOXES ONLY; GLOBAL total photo API counts contextual; no global rail-near estimate",
       "country_bbox_rail_geometry_original_source":spatial_sources,
       "n_source_authentic_original_selected_distinct_species":len(selected),
       "two_country_bbox_20taxon_real_iNat_photo_source_30m":result,
       "global_species_research_photo_totals_2018_2025_NOT_global_railway":global_totals,
       "n_country_species_strata_completed_without_API_cap":sum(x["source_query_complete"] for arr in result.values() for x in arr),
       "n_country_species_strata_api_capped_and_NEAR30_UNESTIMATED_NOT_ZERO":sum(not x["source_query_complete"] for arr in result.values() for x in arr),
       "n_complete_strata_with_geometrically_point_within30m":sum(
            (x["source_original_photo_quality_and_geo_rail_polyline_30m"]["n_point_near_30m_strict_accuracy"] or 0)
            for arr in result.values() for x in arr if x["source_query_complete"]),
       "n_complete_strata_with_photo_license_near30m":sum(
            x["source_original_photo_quality_and_geo_rail_polyline_30m"]["n_point_near_30m_photo_CC_licensed"]
            for arr in result.values() for x in arr if x["source_query_complete"]),
       "NOT_exhaustive_global_rail_30m_count":True,
       "NOT_actual_rail_photographed_plant_phenology_or_height_verified":True,
       "NOT_original_site_specific_true_rail_preference_vs_nonrail":True,
       "source_quality_limitations":[
         "Czech and Slovenian study COUNTRY BOUNDING BOXES overlap their borders; an iNaturalist API bbox is not the national political boundary, so totals are bbox totals.",
         "Source railway geometry is OSM 2023 Czech / 2024 Slovenia, later than some observations 2018-25 and original vegetation surveys 2020-23, including abandoned/low-use source sections.",
         "For within30m counts, public point distance<=30m AND public coordinate accuracy<=30m; only a subset has error circle entirely within 30m.",
         "A photographer near a railway, platform or station garden can have photographed a plant elsewhere; rail 30m geometric proximity is not legal railway habitat nor train-slipstream exposure.",
         "Research-grade exact binomial taxonomy and infra descendant labels require source name checks; unrecognized synonyms not treated as absence or automatically merged.",
         "API above-8000 source bbox strata are explicitly censored and do not count toward near30m totals; no 10k endpoint paging assumption.",
         "Photo existence/Creative Commons license does not establish that a whole plant, measurable height, open flower or first flowering date is visible.",
         "Observer/site concentration and different density of roads, settlements, urban opportunities create powerful sampling bias; no rail-specific taxon preference or within-taxon adaptation inferred.",
         "No private/obscured iNaturalist coordinates or raw IDs/usernames/coordinates or licensed image copies exported.",
         "User-requested full GLOBAL within30m railway census NOT accomplished; requires global rail spatial data and proper iNaturalist bulk/export/DwC ingestion beyond API medium batches."
       ],
    }
    return output

def main():
    p=argparse.ArgumentParser()
    for k in ("candidates","si-osm-pbf","cz-osm-pbf","output"):p.add_argument("--"+k,required=True)
    p.add_argument("--api-delay",type=float,default=API_DELAY_SECONDS)
    a=p.parse_args()
    if a.api_delay<1.:p.error("iNaturalist API must remain below 1 request per second")
    manifest=json.loads(Path(a.candidates).read_text(encoding="utf-8"))
    summary=run(manifest,{"SI":a.si_osm_pbf,"CZ":a.cz_osm_pbf},
                ConservativeINat(delay=a.api_delay))
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "version":summary["version"],
        "n_source_exact_species":summary["n_source_authentic_original_selected_distinct_species"],
        "n_source_country_species_completed":summary["n_country_species_strata_completed_without_API_cap"],
        "n_capped_strata_unestimated_not_zero":summary["n_country_species_strata_api_capped_and_NEAR30_UNESTIMATED_NOT_ZERO"],
        "n_near30m_photo_records_among_complete_strata_not_global":
          summary["n_complete_strata_with_geometrically_point_within30m"],
        "n_near30m_CC_photo_records":
          summary["n_complete_strata_with_photo_license_near30m"],
        "NOT_GLOBAL_near30_count":True
    },ensure_ascii=False))

if __name__=="__main__":main()
