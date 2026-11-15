#!/usr/bin/env python3
"""v16.8 genuinely retrieve photographed iNaturalist records for the four
review-prioritized Czech ORIGINAL railway flora taxa and audit near30-vs-
100–1000m local control opportunities.

This is SOURCE FEASIBILITY for human image QA, NOT a photo-confirmed
within-taxon phenology, height, mowing or heritable adaptation effect.

Original taxon selection: publisher Czech/Slovenia railway flora 20 species
(v16.4); 40 COMPLETE geographically proximate photo source strata (v16.5);
v16.6 source-verified 11 basic photo gates; 10 angiosperm image-review gates
(v16.7); EXACT 4 of them with >=5 explicit Flowers photo annotations,
without treating the remaining annotations as "non-flowering" observations.
This selection occurred after knowing observational photo counts;
DO NOT treat as a confirmatory taxon selection for a railway effect.

Photographic metadata must be real iNaturalist source 2018–25, research
grade, visible photo, public non-obscured species-identified observations,
source coordinate uncertainty 1–30m. Original Czech national 2023 OSM
PBF true railway ways nearest LINE distance (NOT nearest OSM node).
Near for primary review requires point-line distance<=30m and recorded
coordinate-error circle wholly inside 30m line buffer. Off-rail proxy:
photo 100–1000m from nearest rail polyline, same strict accuracy.
Both near and controls have at least one Creative Commons-licensed photo.
This is NOT proof of actual railway track habitat for the near photos nor
a nonrail vegetation census for the control band.

One-to-one deterministic matching:
  * exact same species, calendar year AND month;
  * geographical nearest <=10km great-circle across records;
  * non-identical original iNaturalist observation and original 100m cell;
  * one unique near and unique off-rail record used once each;
  * prefer different observer as a documented sensitivity, not causality.
Source CC photos and observation precise IDs/coordinates never written
to GitHub or exported. Only aggregate number of matched near/control
photos, year-month coverage, unique public observers/cells, and the
difference between same-observer vs different-observer opportunities.
No actual pixels acquired or manually rated. Positive Flowers annotations
are carried only as optional review hints, not phenology labels or first
flowering observations. No image size/height without calibrated scale.
"""
from __future__ import annotations
import argparse,datetime as dt,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

from audit_20_original_rail_plants_inat_30m_two_regions_v165 import (
   RailwayPolylineSource,ConservativeINat,STUDY_REGIONS,
   fixed_query,load_complete_photo_observations,observation_quality,
   PHOTO_LICENSE_ALLOW,original_pbf_authentic,
)
from audit_iNat_rail30m_flower_annotation_review_gate_v167 import (
   V166 as ORIGINAL_166,MIN_ANNOTATED_FLOWER,
)

PRIORITY_VERSION="v16.7_SOURCE_ORIGINAL_20_RAIL_TAXA_2BBOX_FLOWER_ANNOTATION_AND_IMAGE_REVIEW_GATE"
CZ_GIS_SHA="26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"
OBSERVATION_REGION="CZ"
N_SOURCE_PRIORITY=4
CONTROL_MIN_M=100.0
CONTROL_MAX_M=1000.0
MAX_LOCAL_DISTANCE_KM=10.0
DATE_START="2018-01-01"
DATE_END="2025-12-31"
SAMPLING_ALGORITHM="same calendar year+month, <=10km, 1:1, each original observation used only once"

def review_priority(source):
    if (source.get("complete") is not True or source.get("version")!=PRIORITY_VERSION or
        source.get("source_v166_author_photo_gis_ID")!="43693580409" or
        source.get("source_n_real_eligible_photo_gates_v166")!=11 or
        source.get("n_species_country_cells_passing_exploratory_ANGIOSPERM_image_review_gate")!=10):
        raise ValueError("Original source species/photo/CC review index v16.7 has drifted")
    ranked=source["exploratory_ANGIOSPERM_image_review_priority"]
    selected=[r for r in ranked
              if (r["near30_original_iNat_positive_FLOWERS_annotation"]>=MIN_ANNOTATED_FLOWER
                 and r["region"]=="CZ")]
    if (len(selected)!=N_SOURCE_PRIORITY or len({r["species"] for r in selected})!=N_SOURCE_PRIORITY
        or any(r["species"]=="Equisetum arvense" for r in selected)):
        raise ValueError("Cannot safely select original four non-horsetail angiosperms without preexisting independent metadata review gates")
    return selected

def month_from_original_record(ob):
    value=str(ob.get("observed_on") or "")
    try:date=dt.date.fromisoformat(value[:10])
    except ValueError:raise ValueError("Original iNaturalist observation date invalid")
    if not(DATE_START<=date.isoformat()<=DATE_END):
        raise ValueError("Original photographed iNaturalist observation outside fixed 2018-2025 time period")
    return f"{date.year:04}-{date.month:02}"

def has_licensed_photo(ob):
    return any((p.get("license_code") or "").lower() in PHOTO_LICENSE_ALLOW
               for p in (ob.get("photos") or []))

def has_explicit_flowers_tag(ob):
    return any(a.get("controlled_attribute_id")==12 and a.get("controlled_value_id")==13
               for a in (ob.get("annotations") or []))

def haversine_km(a,b):
    lat1,lon1=a;lat2,lon2=b
    p1,p2=math.radians(lat1),math.radians(lat2)
    dlat=p2-p1;dlon=math.radians(lon2-lon1)
    s=math.sin(dlat/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dlon/2)**2
    return 6371.0088*2*math.asin(min(1.0,math.sqrt(s)))

def _private_real_record(ob,rail,taxon):
    # Used transiently ONLY in memory; never serialize ID, observer, lon/lat,
    # photo URL or coarse 100m grid cell to GitHub artifacts.
    coord,reason=observation_quality(ob,taxon)
    if coord is None:return None,reason
    if not has_licensed_photo(ob):
        return None,"original_photo_no_allowed_CC_licence"
    lon,lat,accuracy=coord
    distance,cell=rail.distance_and_cell(lon,lat)
    if distance<=30 and distance+accuracy<=30:
        group="near30_certain"
    elif CONTROL_MIN_M<=distance<=CONTROL_MAX_M:
        group="control_100to1000"
    elif distance<=30:
        return None,"near30_point_only_uncertainty_circle_crosses_buffer"
    else:
        return None,"outside_predeclared_primary_near_or_100to1000m_control_band"
    user=(ob.get("user") or {}).get("id")
    # An anonymous user is still a unique photo opportunity, not counted
    # as an independent observer. That is not zero field sampling effort.
    return {
       "internal_original_id":int(ob["id"]),
       "internal_original_user_id":str(user) if user is not None else None,
       "lat":lat,"lon":lon,"internal_source_100m_cell":cell,
       "year_month":month_from_original_record(ob),
       "group":group,
       "annotated_flowers_positive":has_explicit_flowers_tag(ob),
       "original_source_distance_to_nearest_railway_line_m":distance,
       "original_coordinate_uncertainty_m":accuracy
    },"pass"

def _hash_id(internal_observation_id):
    return hashlib.sha256(("tetsu-v168-photo-pair-logistics-"+
          str(internal_observation_id)).encode("utf-8")).hexdigest()

def one_to_one_matched(near,control,max_km=MAX_LOCAL_DISTANCE_KM,
                       restrict_different_observers=True):
    """For each candidate 1:1 pair, prefer locally nearest control that
    shares a calendar month and year; deterministic hash-sorted near order.
    This is sample availability, NOT balance on habitat/observer/rail use.
    """
    control_by_month=defaultdict(list)
    for c in control:
        control_by_month[c["year_month"]].append(c)
    for group in control_by_month.values():
        group.sort(key=lambda c:_hash_id(c["internal_original_id"]))
    matched=[];used=set()
    for n in sorted(near,key=lambda r:_hash_id(r["internal_original_id"])):
        options=[]
        for c in control_by_month[n["year_month"]]:
            if (c["internal_original_id"] in used or
                c["internal_original_id"]==n["internal_original_id"] or
                c["internal_source_100m_cell"]==n["internal_source_100m_cell"]):
                continue
            if (restrict_different_observers and n["internal_original_user_id"] is not None
                and n["internal_original_user_id"]==c["internal_original_user_id"]):
                continue
            d=haversine_km((n["lat"],n["lon"]),(c["lat"],c["lon"]))
            if d<=max_km:
                options.append((d,_hash_id(c["internal_original_id"]),c))
        if not options:continue
        _,_,c=min(options,key=lambda x:(x[0],x[1]))
        used.add(c["internal_original_id"])
        matched.append((n,c))
    if len(matched)>min(len(near),len(control)):
        raise ValueError("One-to-one source pair match exceeds real photographed records")
    if len({n["internal_original_id"] for n,_ in matched})!=len(matched):
        raise ValueError("Original near source image observation reused")
    if len({c["internal_original_id"] for _,c in matched})!=len(matched):
        raise ValueError("Original distant source image observation reused")
    return matched

def summarize_source_species(obs,rail,species,source_prior):
    original_source_total=int(source_prior["source_bbox_iNat_research_photo_total"])
    if original_source_total!=len(obs):
        raise ValueError("Original API source result count unexpectedly changed versus original source v16.7; do not mix invisible version drifts")
    quality=Counter();near=[];control=[]
    for o in obs:
        record,reason=_private_real_record(o,rail,species)
        quality[reason]+=1
        if record is not None:
            if record["group"]=="near30_certain":near.append(record)
            elif record["group"]=="control_100to1000":control.append(record)
    if len(near)>source_prior["original_public_accuracy_circle_guaranteed_near30"]:
        raise ValueError("Czech original accurate near30 photo count increased beyond original documented source; investigate source data drift")
    if len(near)>source_prior["original_CC_photo_near30"]:
        raise ValueError("Strict source CC photographed 30m group exceeds old source license totals")
    if len(control)>source_prior["original_photo_controls_100to1000m"]:
        raise ValueError("Strict source CC control sample exceeds old original full-photo control")
    matched_diff=one_to_one_matched(near,control,restrict_different_observers=True)
    matched_any=one_to_one_matched(near,control,restrict_different_observers=False)
    positive=sum(n["annotated_flowers_positive"] for n in near)
    pairs_with_positive=sum(n["annotated_flowers_positive"] for n,_ in matched_diff)
    count_observers=lambda records:len({r["internal_original_user_id"] for r in records
                                      if r["internal_original_user_id"] is not None})
    all_sites=lambda records:len({r["internal_source_100m_cell"] for r in records})
    months=Counter(n["year_month"] for n in near)
    return {
      "source_species":species,
      "source_2018_2025_full_real_API_photograph_observation_count":len(obs),
      "n_source_point_near30_and_errorcircle_certain_and_CC_licensed":len(near),
      "n_source_photo_control100to1000m_accuracy30_CC":len(control),
      "n_source_near30_unique_100m_cells":all_sites(near),
      "n_source_near30_unique_public_observers":count_observers(near),
      "n_source_control_unique_100m_cells":all_sites(control),
      "n_source_control_unique_public_observers":count_observers(control),
      "n_same_year_month_unique_photo_pairs_within10km_any_observers":len(matched_any),
      "n_same_year_month_unique_photo_pairs_within10km_different_observers":len(matched_diff),
      "n_near30_source_original_positive_Flowers_annotation_in_eligible_CC_accuracy_circle_group":positive,
      "n_different_observer_local_matched_pairs_with_original_positive_Flowers_annotation_near_rail":pairs_with_positive,
      "n_matched_source_unique_near_100m_cells":all_sites([n for n,_ in matched_diff]),
      "n_matched_source_unique_control_100m_cells":all_sites([c for _,c in matched_diff]),
      "n_matched_source_distinct_months":len({n["year_month"] for n,_ in matched_diff}),
      "n_matched_source_years":len({n["year_month"][:4] for n,_ in matched_diff}),
      "matched_pair_by_calendar_yearmonth":dict(sorted(Counter(n["year_month"] for n,_ in matched_diff).items())),
      "source_photo_rejection_counts":dict(quality),
      "image_pixels_original_photo_and_geo_or_observer_ID_NOT_exposed":True,
      "cannot_identify_flower_negative_trait_from_untagged":True,
      "cannot_infer_actual_flowering_onset_or_photo_plant_height":True
    }

def run(source,rail,client=None):
    selected=review_priority(source)
    client=client or ConservativeINat(delay=1.10)
    outputs=[]
    for rank,row in enumerate(selected,1):
        species=row["species"]
        meta,obs=load_complete_photo_observations(
           client,species,STUDY_REGIONS["CZ"]["bbox"])
        if not meta["complete"] or obs is None:
            raise ValueError("Actual original iNaturalist country taxon primary reviewer photo census capped, never count as 0")
        prior_by_name={x["species"]:x for x in source["original_species_country_40cells_annotation_and_month_audit"]
                       if x["region"]=="CZ"}
        if species not in prior_by_name:
            raise ValueError("New species not source verified in v16.7")
        result=summarize_source_species(obs,rail,species,prior_by_name[species])
        result["original_preselected_author_botanist_taxon_source_v167_rank"]=rank
        outputs.append(result)
        print(json.dumps({"source_rank":rank,
             "species":species,
             "source_API":meta["total"],
             "strict_CC_accuracy_circle30":result["n_source_point_near30_and_errorcircle_certain_and_CC_licensed"],
             "control_CC_accuracy30":result["n_source_photo_control100to1000m_accuracy30_CC"],
             "exact_yearmonth_geo10km_diff_observer_matched":result["n_same_year_month_unique_photo_pairs_within10km_different_observers"]},
             ensure_ascii=False),flush=True)
    return {
       "complete":True,
       "version":"v16.8_GENUINE_Czech_four_angiosperms_local_CC_photo_PAIRS_review_opportunity",
       "original_source_v167_eligibility":"v16.7_SOURCE_ORIGINAL_20_RAIL_TAXA_2BBOX_FLOWER_ANNOTATION_AND_IMAGE_REVIEW_GATE",
       "original_source_20_author_botanist_real_railway_flora_DOIs":[
            "10.1007/s44473-026-00155-z","10.23855/preslia.2024.247"],
       "source_Geofabrik_original_Czech2023_PBF_SHA":CZ_GIS_SHA,
       "original_inat_photo_years":[DATE_START,DATE_END],
       "historical_OSM_railway_nearest_LINE_not_node":True,
       "n_original_source_selected_review_angiosperms":len(outputs),
       "n_source_species_CC_photo_matched_local_calendarmonth":sum(
          x["n_same_year_month_unique_photo_pairs_within10km_different_observers"] for x in outputs),
       "source_restricted_original_4_taxa_matched_opportunity":outputs,
       "selected_AFTER_photo_annotation_counts_not_independent_confounder_control":True,
       "NO_raw_original_photos_observation_ID_public_coordinates_or_usernames_exported":True,
       "study_limits":[
          "This counts potential 1:1 photograph pairs, with same year/month and geodesic<=10km, not measured phenological/morphological differences.",
          "Pairing is deterministic and GREEDY, not optimal matching and not randomized exposure allocation; possible matches do not verify source habitat.",
          "The near primary requires 2023 physical OSM rail-line distance+coordinate error<=30m, yet the image itself may depict a planted specimen or photograph location not plant site.",
          "100-1000m off-rail control is not a road-exclusive or nonrail habitat census; roadside, railway station background and urban habitat may differ.",
          "Only CC photographs are counted to enable later ethical human QA; no images or usernames downloaded/copied in this script.",
          "Only taxa selected after initial photo annotation feasibility are retained. The four original taxa do not constitute a predefined representative sample of railway plants.",
          "An explicit Flowers annotation positive on one or two paired photos does NOT identify first flowering day or true absence in its unannotated partner.",
          "No site-level actual mowing, train wind, local lineage, plant size scale, genotype or fitness measured.",
       ],
       "NOT_causal_railway_phenology_or_genetic_adaptation":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--verified-v167",required=True)
    p.add_argument("--czech-original-osm",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--api-delay",type=float,default=1.10)
    a=p.parse_args()
    if a.api_delay<1.0:p.error("iNaturalist rate limit 1 request per second")
    source=json.loads(Path(a.verified_v167).read_text(encoding="utf-8"))
    original_pbf_authentic(a.czech_original_osm,"CZ")
    railway=RailwayPolylineSource(a.czech_original_osm,"CZ")
    result=run(source,railway,ConservativeINat(delay=a.api_delay))
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
       "complete":True,"version":result["version"],
       "n_source_4_reviewer_species":result["n_original_source_selected_review_angiosperms"],
       "n_same_year_month_local_10km_1to1_CC_photo_pairs_different_observer":
         result["n_source_species_CC_photo_matched_local_calendarmonth"],
       "source_species_control_photo_matches":[{
          "species":v["source_species"],
          "near30_CC_certain":v["n_source_point_near30_and_errorcircle_certain_and_CC_licensed"],
          "control100to1000_CC":v["n_source_photo_control100to1000m_accuracy30_CC"],
          "matched_one_to_one":v["n_same_year_month_unique_photo_pairs_within10km_different_observers"],
          "original_positive_FLOWERS_in_near_photo_matched":
               v["n_different_observer_local_matched_pairs_with_original_positive_Flowers_annotation_near_rail"],
          "n_matched_yearmonths":v["n_matched_source_distinct_months"],
          "n_matched_years":v["n_matched_source_years"],
        } for v in result["source_restricted_original_4_taxa_matched_opportunity"]],
       "no_actual_photo_image_pixels_or_flowering_shift_observed":True
    },ensure_ascii=False))

if __name__=="__main__":main()
