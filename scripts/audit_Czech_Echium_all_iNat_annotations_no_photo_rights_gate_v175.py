#!/usr/bin/env python3
"""v17.5 Authentic Czech Echium iNaturalist public photo ALL-RECORD phenology
annotation audit, with ZERO source image pixels downloaded or redistributed.

Independent candidate context: authenticated CzechVeg source v17.4, original
2,510,340 species×relevé from 117,739 plot headers, same source study,
original phytoregion, decade and plot area class RAIL-text vs ROAD-text
observational prevalence for Echium +0.0690920. This is a selection proxy
from DIFFERENT ORIGINAL botanist data, not genuine railway selection.

All 4,235 original 2018-25 Czech bbox research-grade PHOTO observations
of Echium vulgare, as independently authenticated v16.6: no CC0/BY or
other image license restrictions, no 1:1 year/month/10km matching.
Read PUBLIC metadata only; never GET an image, store URLs, photographer
identifiers, precise GPS, iNat original ID or sensitive observations.
Historic immutable 2023 Czech Geofabrik RAIL LINE PBF source verified SHA.

For original public un-obscured photo observations, report all high quality
coordinate and positional uncertainty strata including >30m/unknown
accuracies. SOURCE STRICT PRIMARY: reported public pos accuracy 1–30m,
photo location point <=30m rail line, error-circle within30m separately.
Source CONTROL: point 100–1000m, and entire public positional error
circle within [100,1000] separately. No physical rail right-of-way claim.

Official current iNaturalist annotations (Photos & Fruits group 12):
  13=Flowers (at least one OPEN flower VISIBLE)
  21=No Flowers or Fruits (no open flowers/fruits/seeds VISIBLE)
  15=Flower Buds, 14=Fruits or Seeds (not inferred as no open flowers).
  unannotated != no flowers and != plant not flowering.
  13+21 contradiction excluded. Pos and Neg refer to PHOTO evidence,
  NOT actual blooming timing, true nonflowering physiology or fitness.
Source: https://help.inaturalist.org/en/support/solutions/articles/151000191830
  and .../151000170738.

Gate of phenotype inferability PRE-ESTABLISHED in this version:
  >=50 exclusively positive- or negative-labelled observations in EACH
  exposure stratum, of which >=10 have Flowers and >=10 have explicit
  No Flowers or Fruits, >=5 distinct observational calendar months.
If not met, label logistic model is NOT IDENTIFIABLE; no imputation of
negative labels from absent annotation and no regression p-values.
Even if met, selection into annotation and local land-use confounding
prevent causal railway flowering inference, require external QA.

Export only country/category/annotation counts; NO plant visual scoring
claimed. This is NOT an "all Europe rail" regression; no global rail
geometry nor pan-European photograph denominator is yet built.
"""
from __future__ import annotations
import argparse,datetime,json,math
from collections import Counter,defaultdict
from pathlib import Path

from audit_20_original_rail_plants_inat_30m_two_regions_v165 import (
    RailwayPolylineSource, ConservativeINat, STUDY_REGIONS,
    load_complete_photo_observations, observation_quality,
    exact_taxon_or_descendant, original_pbf_authentic
)

TARGET="Echium vulgare"
EXPECTED_TOTAL=4235
SOURCE_V166="v16.6_AUTHENTIC_2COUNTRY_20_ORIGINAL_RAIL_TAXA_PUBLIC_PHOTO_30m_FEASIBILITY_QA"
PBF_SHA="26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"
ANNOTATED_MIN_EACH_GROUP=50
ANNOTATED_MIN_PER_CLASS=10
MIN_MONTHS=5
ANNO_TERM=12
ANNO_VALS={"photo_open_flower_visible":13,"photo_no_flower_fruit_seed_visible":21,
           "photo_buds_visible":15,"photo_fruits_or_seeds_visible":14}

def validate_prior(prior):
    if prior.get("version")!=SOURCE_V166 or prior.get("complete") is not True:
        raise ValueError("Original author 20 source photo rail 30m geography audit v16.6 missing")
    src=prior["source_species_country_counts"]
    x=[r for r in src if r["taxon"]==TARGET and r["region_source_bbox"]=="CZ"]
    if len(x)!=1 or x[0]["source_bbox_iNat_research_photo_total"]!=EXPECTED_TOTAL:
        raise ValueError("Original full 4,235 Czech source iNaturalist Echium photos changed")
    if prior["source_original_publications_and_Geofabrik_SHA"]["CZ"]["original_Geofabrik_PBF_SHA256"]!=PBF_SHA:
        raise ValueError("Source true original CZ2023 OSM PBF fingerprint mismatch")
    return x[0]

def annotation_class(ob):
    ann={a.get("controlled_value_id") for a in (ob.get("annotations") or [])
         if a.get("controlled_attribute_id")==ANNO_TERM}
    if 13 in ann and 21 in ann:return "contradictory_photo_annotations_exclude"
    if 13 in ann:return "photo_open_flower_visible"
    if 21 in ann:return "photo_no_flower_fruit_seed_visible"
    if 15 in ann and 14 in ann:return "photo_buds_and_fruits_visible_but_open_flower_unknown"
    if 15 in ann:return "photo_buds_visible_but_open_flower_unknown"
    if 14 in ann:return "photo_fruits_visible_but_open_flower_unknown"
    if ann:return "other_photo_annotation_unresolved"
    return "unannotated_open_flower_UNKNOWN_NOT_NONFLOWERING"

def public_metadata_gate(ob):
    """Do not require photo asset reuse rights or <30m accuracy.
    Strict prior source acceptance is checked SEPARATELY, unchanged.
    """
    if ob.get("obscured") or ob.get("geoprivacy") in ("obscured","private") or ob.get("taxon_geoprivacy") in ("obscured","private"):
        return None,"obscured_or_private"
    if ob.get("captive") is True:
        return None,"captive_or_cultivated"
    if not exact_taxon_or_descendant((ob.get("taxon") or {}).get("name"),TARGET):
        return None,"taxon_mismatch"
    if not ob.get("photos"):
        return None,"no_photos"
    date=ob.get("observed_on")
    try:dt=datetime.date.fromisoformat(str(date)[:10])
    except (ValueError,TypeError):
        return None,"date_missing_or_invalid"
    if not(2018<=dt.year<=2025):
        return None,"date_outside_2018to2025"
    coords=(ob.get("geojson") or {}).get("coordinates")
    if not isinstance(coords,list) or len(coords)!=2 or any(
            isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v)
            for v in coords):
        return None,"missing_or_invalid_public_coordinates"
    lon,lat=coords
    if not(-180<=lon<=180 and -90<=lat<=90):
        return None,"invalid_lat_lon_public"
    acc=ob.get("positional_accuracy")
    if isinstance(acc,bool) or not isinstance(acc,(int,float)) or not math.isfinite(acc) or acc<=0:
        acc=None
        acu="accuracy_unknown"
    elif acc<=30:acu="accuracy_1to30m"
    elif acc<=100:acu="accuracy_30to100m"
    elif acc<=1000:acu="accuracy_100to1000m"
    else:acu="accuracy_above1000m"
    return {"lon":float(lon),"lat":float(lat),
            "accuracy_m":float(acc) if acc is not None else None,
            "accuracy_class":acu,"month":f"{dt.month:02d}",
            "year":dt.year}, "pass"

def distance_zone(d):
    if d<=30:return "point_0to30m"
    if d<=100:return "point_30to100m"
    if d<=1000:return "point_100to1000m"
    return "point_over1000m"

def phenotype_gate(counts,months):
    def available(z):
        pos=counts[z]["photo_open_flower_visible"]
        neg=counts[z]["photo_no_flower_fruit_seed_visible"]
        n=pos+neg
        k=sum(v>=1 for v in months[z].values())
        return {"n_exclusive_positive":pos,"n_explicit_photo_negative":neg,
                "n_true_labeled_photo_outcomes":n,"n_distinct_months_any_photos":k,
                "metadata_photo_label_model_gate_passed":n>=ANNOTATED_MIN_EACH_GROUP and
                        pos>=ANNOTATED_MIN_PER_CLASS and neg>=ANNOTATED_MIN_PER_CLASS and
                        k>=MIN_MONTHS}
    return {"near30_accuracy30":available("strict_near30"),
            "control100to1000_accuracy30":available("strict_control100to1000")}

def process_records(obs,rail,prior_row):
    if len(obs)!=EXPECTED_TOTAL:
        raise ValueError("Full original source Czech 2018–25 photo cohort not 4235")
    all_public=Counter();quality=Counter()
    zones=defaultdict(Counter);months=defaultdict(Counter)
    precision=defaultdict(Counter)
    original_strict_near=0;original_strict_guaranteed=0;original_strict_control_point=0
    for o in obs:
        coords,q=public_metadata_gate(o)
        quality[q]+=1
        if coords is None:continue
        lon,lat,acc,acat=coords["lon"],coords["lat"],coords["accuracy_m"],coords["accuracy_class"]
        dist,_=rail.distance_and_cell(lon,lat)
        region=distance_zone(dist)
        anno=annotation_class(o)
        zones[region][anno]+=1
        zones["all_public_usable"][anno]+=1
        precision[region][acat]+=1
        precision["all_public_usable"][acat]+=1
        months[region][coords["month"]]+=1
        if acat=="accuracy_1to30m":
            strict,q_strict=observation_quality(o,TARGET)
            if strict is None or q_strict!="pass":
                raise ValueError("Original v16.6 photo strict GPS quality logic not recovered")
            if region=="point_0to30m":
                original_strict_near+=1
                zones["strict_near30"][anno]+=1
                months["strict_near30"][coords["month"]]+=1
                if dist+acc<=30:
                    original_strict_guaranteed+=1
                    zones["guaranteed_circle_near30"][anno]+=1
                    months["guaranteed_circle_near30"][coords["month"]]+=1
            if region=="point_100to1000m":
                original_strict_control_point+=1
                zones["strict_control100to1000"][anno]+=1
                months["strict_control100to1000"][coords["month"]]+=1
                if dist-acc>=100 and dist+acc<=1000:
                    zones["guaranteed_circle_control100to1000"][anno]+=1
                    months["guaranteed_circle_control100to1000"][coords["month"]]+=1
        all_public[anno]+=1
    expected=prior_row
    actual=(original_strict_near,original_strict_guaranteed,original_strict_control_point)
    known=(expected["n_rail_polyline_point_near_30m"],
           expected["n_rail_polyline_guaranteed_circle_within_30m"],
           expected["n_source_species_photo_100to1000m_controls"])
    if actual!=known:
        raise ValueError("Original v16.6 strict photo+OSM source counts changed: "+
                         str(actual)+" != "+str(known))
    if sum(quality.values())!=EXPECTED_TOTAL or sum(all_public.values())!=quality["pass"]:
        raise ValueError("Source full photo API quality/annotation accounting not conserved")
    gate=phenotype_gate(zones,months)
    result={
       "complete":True,"version":"v17.5_REAL_CZ_Echium_full_4235_iNat_metadata_NO_image_rights_filter",
       "source_original_record_DOI":"10.23855/preslia.2024.247",
       "source_independent_CzechVeg_context_DOI":"10.5281/zenodo.19231152",
       "source_verified_CzechVeg_context_rail_minus_road_per_plot_2018to2025_photo_COHORT_NOT_same_study":0.06909200397868588,
       "original_Geofabrik_2023_CZ_true_rail_LINE_SHA256":PBF_SHA,
       "source_2018_2025_Czech_bbox_research_grade_photo_API_records":EXPECTED_TOTAL,
       "full_original_API_photo_record_count_and_strict_geography_v166_reproduced":True,
       "photo_original_source_privacy_taxon_and_GPS_quality":dict(quality),
       "photo_annotation_codes_OFFICIAL_term12":{
          "13":"at least one OPEN flower is visible",
          "21":"NO flower fruit or seed visible in photo, NOT biological nonflowering",
          "15":"flower BUDS visible; not a true open-flower negative",
          "14":"fruit/seeds visible; not a true open-flower negative",
          "missing":"UNKNOWN and NOT NONFLOWERING",
       },
       "all_public_observations_annotation_outcomes":dict(all_public),
       "all_coordinate_accuracies_NOT_for_true_rail_exposure":{k:dict(v) for k,v in precision.items()},
       "original_railway_polyline_distance_zone_anthropogenic_observation_metadata":{
          k:dict(v) for k,v in zones.items()},
       "original_observation_MONTH_counts_NOT_FIRST_FLOWERING":{k:dict(sorted(v.items())) for k,v in months.items()},
       "fixed_2018_2025_label_coverage_thresholds":{
          "minimum_explicit_labeled_photos_per_group":ANNOTATED_MIN_EACH_GROUP,
          "minimum_positive_and_negative_each_per_group":ANNOTATED_MIN_PER_CLASS,
          "minimum_distinct_original_photo_observation_months":MIN_MONTHS
       },
       "source_observed_annotated_photo_phenotype_regression_eligibility":gate,
       "source_model_identifiable_for_two_photo_distance_groups":all(
            v["metadata_photo_label_model_gate_passed"] for v in gate.values()),
       "no_actual_image_pixels_downloaded":True,
       "NO_photo_image_copyright_restriction_used_for_METADATA":True,
       "no_1to1_matching_or_record_discard_for_yearmonth10km":True,
       "no_unannotated_photo_coded_negative":True,
       "no_true_plant_first_flowering_date_or_genetic_mowing_adaptation_estimated":True,
       "source_observation_limitations":[
         "All 4235 original public photo observation records are metadata census, with separate privacy, taxon and accuracy quality; cannot regress flower phenotype unless explicitly annotated positive/negative source support.",
         "The official negative annotation 21 says NO flower, fruit or seed visible in image, not the original plant biologically nonflowering or even necessarily devoid of buds.",
         "A photo with Flower Buds 15 or Fruits and Seeds 14 only cannot be treated as an independent no-open-flowers negative without actual pixel review.",
         "Distance is to fixed original 2023 Czech OSM rail LINE and 2018-25 observations; plants may be near platforms, roads, gardens, or rail embankment not directly under railway management.",
         "iNaturalist photography selection and annotation behavior differ with geography, season, observer, habitat and distance; no causally identifiable rail exposure.",
         "This examines only one source species and Czech bbox, NOT a global Europe flower phenotype regression.",
         "CzechVeg matched rail/road textual vegetation context is independent of these photographs, but different study/source plot grain and spatial mismatch preclude combining their 6.91pp contrast as a within-species effect.",
         "No raw original source image pixels, public observation IDs, usernames, original photo URLs or individual precise coordinates included."
       ]
    }
    return result

def run(v166,rail,client=None):
    prior=validate_prior(v166)
    client=client or ConservativeINat(delay=1.1)
    meta,obs=load_complete_photo_observations(client,TARGET,STUDY_REGIONS["CZ"]["bbox"])
    if not meta["complete"] or meta["total"]!=EXPECTED_TOTAL or obs is None:
        raise ValueError("Current real iNaturalist source cohort not fully retrievable; cannot pretend no photos")
    return process_records(obs,rail,prior)

def main():
    p=argparse.ArgumentParser()
    for x in ("original-v166","czech-2023-OSM","output"):
        p.add_argument("--"+x,required=True)
    a=p.parse_args()
    v166=json.loads(Path(a.original_v166).read_text(encoding="utf-8"))
    original_pbf_authentic(a.czech_2023_OSM,"CZ")
    rail=RailwayPolylineSource(a.czech_2023_OSM,"CZ")
    out=run(v166,rail)
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
      "version":out["version"],
      "n_source_original_iNat_photo_records":EXPECTED_TOTAL,
      "n_public_quality_gated":out["photo_original_source_privacy_taxon_and_GPS_quality"].get("pass",0),
      "annotated_source_full_public":out["all_public_observations_annotation_outcomes"],
      "near30_source_all_photo_annotations":out["original_railway_polyline_distance_zone_anthropogenic_observation_metadata"].get("strict_near30",{}),
      "control100to1000_source_all_photo_annotations":out["original_railway_polyline_distance_zone_anthropogenic_observation_metadata"].get("strict_control100to1000",{}),
      "both_group_identifiable":out["source_model_identifiable_for_two_photo_distance_groups"],
      "source_explicit_per_group_positive_negative":out["source_observed_annotated_photo_phenotype_regression_eligibility"],
      "no_image_copyright_filter_or_false_unannotated_negative":True
    },ensure_ascii=False))

if __name__=="__main__":main()
