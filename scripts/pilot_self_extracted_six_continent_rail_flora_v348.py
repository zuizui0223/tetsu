#!/usr/bin/env python3
"""v348: first-party six-continent OSM rail line × iNat original plant observation audit.

A railway-*proximity* citizen-science audit; NOT trackbed vegetation survey,
railway-specific selection, same-area plots, plant temperature or adaptation.
Source design frozen in docs/tetsu_v348_self_extracted_global_rail_proximity_frozen_design_ja.md
Output: city-level aggregate counts only; never site GPS, raw species names,
iNaturalist IDs, usernames, photos, original precise rail lines or rare plants.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re,time,urllib.parse,urllib.request,xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError,URLError

CITIES=(
    ("Tokyo","Asia","Japan",35.6812,139.7671),
    ("Prague","Europe","Czechia",50.0839,14.4350),
    ("Cape Town","Africa","South Africa",-33.9258,18.4230),
    ("Sydney","Oceania","Australia",-33.8830,151.2070),
    ("San Francisco","North America","USA",37.7764,-122.3940),
    ("Sao Paulo","South America","Brazil",-23.5357,-46.6340),
)
D=0.012
YEAR="2025"
NEAR_M=80.0
FAR_MIN_M=300.0
FAR_MAX_M=1000.0
MAX_POSITIONAL_ACCURACY_M=30.0
MAX_OBSERVATIONS=3600
PAGE=200
RATE_DELAY_S=1.05
MAX_SOURCE_BYTES=28_000_000
USER_AGENT="tetsu-original-rail-ecology-research/1.0 (public open-data source QA; no photos)"
OVERPASS=(
    "https://overpass.kumi.systems/api/interpreter",
    "https://overpass-api.de/api/interpreter",
)
INAT="https://api.inaturalist.org/v1/observations"
OSM_API="https://api.openstreetmap.org/api/0.6/map"
TAXON_ID=47126
MIN_SUPPORT_NEAR=(10,5)
MIN_SUPPORT_BACKGROUND=(30,10)


def bounds(lat,lon):
    return (lat-D,lon-D,lat+D,lon+D)


def get_bytes(url,timeout=35):
    req=urllib.request.Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json,application/xml,text/xml,*/*"})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        b=r.read(MAX_SOURCE_BYTES+1)
        if len(b)>MAX_SOURCE_BYTES:
            raise ValueError("Source body exceeds fixed bound; not a valid complete source")
        return b,r.geturl(),r.status


def eligible_rail(tags):
    return (tags.get("railway")=="rail"
            and tags.get("tunnel") not in ("yes","true","1")
            and tags.get("bridge") not in ("yes","true","1")
            and tags.get("service") not in ("yard","siding","spur")
            and tags.get("disused") not in ("yes","true","1")
            and tags.get("abandoned") not in ("yes","true","1"))


def overpass_rail_bytes(blob):
    d=json.loads(blob)
    if not isinstance(d,dict) or not isinstance(d.get("elements"),list):
        raise ValueError("Overpass original JSON lacks elements")
    ways=[];total=0
    for elem in d["elements"]:
        if elem.get("type")!="way":continue
        tags=elem.get("tags") or {}
        if tags.get("railway")=="rail":total+=1
        if not eligible_rail(tags):continue
        vertices=[(float(p["lat"]),float(p["lon"])) for p in elem.get("geometry",())
                  if "lat" in p and "lon" in p]
        if len(vertices)>=2:ways.append(vertices)
    return ways,{"source_rail_way_count_before_tag_filters":total,
                 "source_rail_ways_eligible":len(ways),
                 "source_rail_segments_eligible":sum(len(v)-1 for v in ways)}


def osm_xml_rail_bytes(blob):
    if not blob.lstrip().startswith(b"<?xml") and not blob.lstrip().startswith(b"<osm"):
        raise ValueError("OSM source is not XML")
    root=ET.fromstring(blob)
    if root.tag!="osm":raise ValueError("Not original OSM map root")
    nodes={e.attrib["id"]:(float(e.attrib["lat"]),float(e.attrib["lon"]))
           for e in root.findall("node")}
    ways=[];total=0
    for e in root.findall("way"):
        tags={tag.get("k"):tag.get("v") for tag in e.findall("tag")}
        if tags.get("railway")=="rail":total+=1
        if not eligible_rail(tags):continue
        pts=[nodes[n.attrib["ref"]] for n in e.findall("nd") if n.attrib["ref"] in nodes]
        if len(pts)>=2:ways.append(pts)
    return ways,{"source_rail_way_count_before_tag_filters":total,
                 "source_rail_ways_eligible":len(ways),
                 "source_rail_segments_eligible":sum(len(v)-1 for v in ways)}


def source_rail(lat,lon):
    a,b,c,d=bounds(lat,lon)
    source_errors=[]
    query=f'[out:json][timeout:45];way["railway"="rail"]({a:.5f},{b:.5f},{c:.5f},{d:.5f});out geom;'
    for endpoint in OVERPASS:
        url=endpoint+"?"+urllib.parse.urlencode({"data":query})
        try:
            blob,final,status=get_bytes(url,timeout=55)
            rail,qa=overpass_rail_bytes(blob)
            if not rail:source_errors.append({"provider":endpoint,"error":"no eligible rail ways"});continue
            return rail,{"provider":"OSM_Overpass","source_endpoint":endpoint,
                         "source_http":status,"source_bytes":len(blob),
                         "source_response_sha256":hashlib.sha256(blob).hexdigest(),**qa,
                         "source_errors_before_success":source_errors}
        except (HTTPError,URLError,TimeoutError,OSError,ValueError,KeyError,ET.ParseError) as e:
            source_errors.append({"provider":endpoint,"error_type":type(e).__name__,
                                   "status":getattr(e,"code",None),"short_message":str(e)[:140]})
    # Official OSM primary API is not a bulk method; small predeclared bbox only.
    osm_url=OSM_API+"?"+urllib.parse.urlencode({"bbox":f"{b:.5f},{a:.5f},{d:.5f},{c:.5f}"})
    try:
        blob,final,status=get_bytes(osm_url,timeout=75)
        rail,qa=osm_xml_rail_bytes(blob)
        if rail:
            return rail,{"provider":"OSM_primary_map_XML_small_bbox_fallback",
                         "source_endpoint":OSM_API,
                         "source_http":status,"source_bytes":len(blob),
                         "source_response_sha256":hashlib.sha256(blob).hexdigest(),**qa,
                         "source_errors_before_success":source_errors}
        source_errors.append({"provider":OSM_API,"error":"no eligible railway=rail ways"})
    except (HTTPError,URLError,TimeoutError,OSError,ValueError,KeyError,ET.ParseError) as e:
        source_errors.append({"provider":OSM_API,"error_type":type(e).__name__,
                               "status":getattr(e,"code",None),"short_message":str(e)[:140]})
    return [],{"provider":"OSM_failed","source_errors":source_errors,
               "source_rail_ways_eligible":0,"source_rail_segments_eligible":0}


def project(lat,lon,lat0,lon0):
    radius=6371008.8
    return ((math.radians(lon-lon0)*math.cos(math.radians(lat0))*radius),
            math.radians(lat-lat0)*radius)


def point_segment_distance(p,a,b):
    px,py=p; ax,ay=a; bx,by=b
    dx=bx-ax;dy=by-ay
    length=dx*dx+dy*dy
    if length==0:return math.hypot(px-ax,py-ay)
    t=max(0.0,min(1.0,((px-ax)*dx+(py-ay)*dy)/length))
    return math.hypot(px-ax-t*dx,py-ay-t*dy)


def rail_segments(ways,lat0,lon0):
    seg=[]
    for way in ways:
        p=[project(a,b,lat0,lon0) for a,b in way]
        seg.extend(zip(p,p[1:]))
    return seg


def min_rail_distance(lat,lon,lat0,lon0,segments):
    if not segments:return None
    p=project(lat,lon,lat0,lon0)
    return min(point_segment_distance(p,a,b) for a,b in segments)


def classify(distance,uncertainty):
    if distance is None or uncertainty is None:return "unknown"
    if not math.isfinite(distance) or not math.isfinite(uncertainty) or distance<0 or uncertainty<0:return "unknown"
    if distance+uncertainty<=NEAR_M:return "near"
    if distance-uncertainty>=FAR_MIN_M and distance+uncertainty<=FAR_MAX_M:return "background"
    return "ambiguous_buffer"


def plant_binomial(taxon):
    if not isinstance(taxon,dict) or taxon.get("rank")!="species":return None
    name=" ".join(str(taxon.get("name") or "").strip().lower().split())
    return name if re.fullmatch(r"[a-z-]+ [a-z-]+",name) else None


def quality_record(r):
    if not isinstance(r,dict):return None,"malformed"
    if r.get("obscured") is True or r.get("geoprivacy") not in (None,"open"):
        return None,"obscured_or_geoprivate"
    if r.get("captive") is True:return None,"captive"
    if r.get("quality_grade")!="research":return None,"not_research"
    taxon=plant_binomial(r.get("taxon"))
    if not taxon:return None,"not_exact_species_rank"
    accuracy=r.get("positional_accuracy")
    if isinstance(accuracy,bool):return None,"accuracy_missing"
    try:unc=float(accuracy)
    except (TypeError,ValueError):return None,"accuracy_missing"
    if not math.isfinite(unc) or unc<0 or unc>MAX_POSITIONAL_ACCURACY_M:
        return None,"accuracy_above_30m_or_invalid"
    geo=r.get("geojson") or {}
    xy=geo.get("coordinates") if isinstance(geo,dict) else None
    if not isinstance(xy,list) or len(xy)<2:return None,"no_public_georeferenced_geometry"
    try:lon,lat=map(float,xy[:2])
    except (TypeError,ValueError):return None,"bad_geocoordinates"
    if not math.isfinite(lon) or not math.isfinite(lat) or not(-180<=lon<=180 and -90<=lat<=90):
        return None,"bad_geocoordinates"
    return (taxon,lat,lon,unc),None


def inat_parameters(lat,lon,page=1):
    a,b,c,d=bounds(lat,lon)
    return {"taxon_id":TAXON_ID,"quality_grade":"research","captive":"false",
            "geo":"true","d1":YEAR+"-01-01","d2":YEAR+"-12-31",
            "swlat":f"{a:.5f}","swlng":f"{b:.5f}","nelat":f"{c:.5f}",
            "nelng":f"{d:.5f}","per_page":PAGE,"page":page,"order":"asc",
            "order_by":"id"}


def source_observations(lat,lon):
    body_hash=hashlib.sha256()
    base=INAT+"?"+urllib.parse.urlencode(inat_parameters(lat,lon,1))
    try:
        b,_,status=get_bytes(base,timeout=55)
        body_hash.update(b)
        js=json.loads(b)
        n=int(js.get("total_results",-1))
        if n<0 or not isinstance(js.get("results"),list):
            raise ValueError("Original iNaturalist results/total_results schema absent")
        qa={"provider":"iNaturalist_public_observations","source_HTTP":status,
            "frozen_year":YEAR,"reported_total_results":n,
            "max_per_window_observations":MAX_OBSERVATIONS,
            "complete_source_pagination":n<=MAX_OBSERVATIONS,
            "source_query_without_2025_results":n==0,
            "API_query_parameters":{k:v for k,v in inat_parameters(lat,lon).items()
                                    if k not in ("swlat","swlng","nelat","nelng")},
            "source_first_page_prefix_hash":hashlib.sha256(b).hexdigest()}
        if n>MAX_OBSERVATIONS:
            qa["stop_reason"]="original_year_count_over_frozen_census_cap"
            return [],qa
        all_results=js["results"].copy()
        pages=math.ceil(n/PAGE)
        if len(all_results)!=min(n,PAGE):
            raise ValueError("iNat page1 result count inconsistent with API reported total")
        for page in range(2,pages+1):
            time.sleep(RATE_DELAY_S)
            url=INAT+"?"+urllib.parse.urlencode(inat_parameters(lat,lon,page))
            raw,_,_=get_bytes(url,timeout=60)
            body_hash.update(raw)
            data=json.loads(raw)
            if not isinstance(data.get("results"),list):raise ValueError("iNat pagination results absent")
            all_results.extend(data["results"])
        if len(all_results)!=n:raise ValueError("Source changed while being paged; fail closed")
        qa.update({"source_pages":max(1,pages),
                   "source_result_count_verified":len(all_results),
                   "aggregate_response_sha256":body_hash.hexdigest()})
        return all_results,qa
    except (HTTPError,URLError,TimeoutError,ValueError,KeyError,OSError) as e:
        return [],{"provider":"iNaturalist_public_observations","source_error_type":type(e).__name__,
                   "source_error_status":getattr(e,"code",None),
                   "short_error":str(e)[:160],
                   "complete_source_pagination":False}


def city_audit(city):
    cityname,continent,country,lat0,lon0=city
    ways,rail_qa=source_rail(lat0,lon0)
    # Respect public API traffic after Overpass; one request per second max.
    time.sleep(RATE_DELAY_S)
    observations,inat_qa=source_observations(lat0,lon0)
    summary={"city":cityname,"country":country,"continent":continent,
             "original_year":YEAR,"anchor_geo_provenance":"predeclared urban-center approx anchor; not site visit",
             "window_lat_lon_halfspan_degrees":D,
             "rail_source":rail_qa,"plant_observation_source":inat_qa}
    if not ways or not inat_qa.get("complete_source_pagination"):
        return {**summary,"status":"source_unavailable_not_negative_biology",
                "near_observations":None,"background_observations":None,
                "rail_proximity_effect_NOT_estimated":True}
    segments=rail_segments(ways,lat0,lon0)
    accepted=Counter()
    counts=Counter()
    observed_species={"near":set(),"background":set()}
    for rec in observations:
        item,error=quality_record(rec)
        if error:counts[error]+=1;continue
        species,lat,lon,u=item
        d=min_rail_distance(lat,lon,lat0,lon0,segments)
        status=classify(d,u)
        counts[status]+=1
        if status in observed_species:
            observed_species[status].add(species)
    near,far=counts["near"],counts["background"]
    minnear,minnearsp=MIN_SUPPORT_NEAR
    minfar,minfarsp=MIN_SUPPORT_BACKGROUND
    nearsp=len(observed_species["near"]);farsp=len(observed_species["background"])
    enough=(near>=minnear and nearsp>=minnearsp and far>=minfar and farsp>=minfarsp)
    summary.update({"status":"proximity_both_arms_observed_DESCRIPTIVE_ONLY" if enough
                   else "insufficient_proximity_support_NOT_a_biological_zero",
                   "near_observations":near,"near_distinct_observed_species":nearsp,
                   "background_observations":far,"background_distinct_observed_species":farsp,
                   "record_exclusion_counts":dict(sorted(counts.items())),
                   "any_rail_near_source_support":near>0,
                   "proximity_two_group_sample_support":enough,
                   "near_background_distance_guards_m":{
                       "near_d_plus_accuracy_max":NEAR_M,"background_d_minus_accuracy_min":FAR_MIN_M,
                       "background_d_plus_accuracy_max":FAR_MAX_M},
                   "no_full_community_inventory_no_absence_inferred":True,
                   "no_species_names_original_GPS_or_observer_ID_exported":True,
                   "rail_proximity_effect_NOT_estimated":True})
    return summary


def run():
    cities=[city_audit(c) for c in CITIES]
    return {"version":"v348-self-extracted-OSM-iNaturalist-six-continent-first-source-wave",
            "source_design_pinned_before_outcomes":True,
            "year":YEAR,"regions_n_predeclared":len(CITIES),"regions":cities,
            "successful_source_regions":sum(x["status"] not in (
              "source_unavailable_not_negative_biology",) for x in cities),
            "regions_with_observed_near_and_background_support":sum(
              x.get("proximity_two_group_sample_support",False) for x in cities),
            "none_of_these_results_are_actual_rail_vegetation_releves_or_causal_effects":True,
            "globally_consistent_species_thermal_preference_NOT_calculated":True,
            "no_original_geographic_records_or_taxon_names_republished":True,
            "scientific_next_step":"Sampling-stratified standardized true railway and locally paired nonrail vegetation needed before causal/thermal meta-analysis",
            "source_licenses":"OSM ODbL; iNaturalist public observation terms and image rights remain with contributors"}

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True)
    a=p.parse_args()
    result=run()
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
