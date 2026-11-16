#!/usr/bin/env python3
"""v17.6: genuine source-requery reproducible PUBLIC iNaturalist observation
review manifest for non-redistributed Echium Czech photo phenotyping.

User's revised design: DO NOT filter photos by CC0/CC BY because no photo
pixels are exported. DO NOT force one-to-one same-month<=10km observation
pair matching. Use all source 2018-25 research-grade Czech bbox 4,235
photographed observations, with source-valid public location/species/date
4,176 eligible; manual-photo labeling may happen later and is BLANK today.
Source authentic Czech 2023 OSM country PBF SHA verified, source v16.6
strict 178 near30/725 control100–1000 replicated, and source v17.5
full-API metadata/annotation counts checksum verified BEFORE exporting.

Deliverable is an ephemeral GitHub Actions ARTIFACT (not Git-tracked)
of public observation IDs (these link back to iNat, where original public
location and credited copyrighted photos remain with their authors) and
blank human visual labels. No copyrighted photos, protected geolocations,
source user handles or exact coordinate columns are exported. The reviewer
is NOT truly blind if the iNat page is opened and its map seen. Treat this
as distance/annotation metadata withheld from the REVIEW WORKSHEET, not
a fully blinded experimental assessor. Do not claim unseen pixels scored.

Two linked files by opaque REVIEW token:
- review_queue.csv: original PUBLIC observation id/permalink, fields for
  flower OPEN visible (YES/NO/UNCERTAIN), whether target whole plant,
  image visually readable, do NOT infer first flowering/height.
  Source annotation/exposure not present to reduce confirmation bias.
- analysis_context.csv: source original calendar year/month/DOY, original
  public GPS accuracy category, original OSM true rail LINE distance band,
  10km COARSENED locality grid, SHA-pseudonymous public observer token,
  explicit source iNat annotation 13/21/buds/fruits/UNKNOWN. This is a
  proxy for sampling and eventual post-reviewed trait modeling.
  Does NOT contain exact GPS, original observation ID or media URL.
- first_wave_review_queue.csv: all 178 public-coordinate-accurate <=30m
  records and deterministic random-like (SHA-based) up to 178 of
  source-accurate public observations 100–1000m from rail. Selection does
  not use annotation status, flowers, license, photo image or outcome,
  and does not force one-to-one pairs. More off-rail comparison records
  and all other strict and loose source records preserved in full queue.

The manifest is explicitly PUBLIC-OBSERVATION-ONLY: a public observation
URL may itself reveal an already-public observation locality; never
include geoprivacy obscured/private records. The original metadata
time and distance are from dynamic sources; version via verified
source run ID. No actual railway-right-of-way or mowing exposure claim.

No genetic adaptation, trait selection, open flower prevalence, fitness
or global Europe photo regression inferred before independent human
labels and spatial/seasonal and observer observation-effort modeling.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter
from pathlib import Path

from audit_20_original_rail_plants_inat_30m_two_regions_v165 import (
   RailwayPolylineSource, ConservativeINat, STUDY_REGIONS,
   load_complete_photo_observations, original_pbf_authentic
)
from audit_Czech_Echium_all_iNat_annotations_no_photo_rights_gate_v175 import (
   TARGET,EXPECTED_TOTAL,public_metadata_gate,annotation_class,
   distance_zone,process_records
)

PRIOR_VERSION="v17.5_REAL_CZ_Echium_full_4235_iNat_metadata_NO_image_rights_filter"
FIRST_WAVE_MAX_CONTROL=178
REVIEW_LABELS=("yes","no","uncertain")

def token_id(public_id):
    if isinstance(public_id,bool) or not isinstance(public_id,int) or public_id<=0:
        raise ValueError("Original iNaturalist public observation numeric identifier invalid")
    return hashlib.sha256(("tetsu-real-Echium-photo-visible-review-v176-"+str(public_id)).encode("utf-8")).hexdigest()[:20]

def verify_source(current,prior):
    if (prior.get("version")!=PRIOR_VERSION or prior.get("complete") is not True or
        prior.get("source_2018_2025_Czech_bbox_research_grade_photo_API_records")!=EXPECTED_TOTAL or
        prior.get("no_actual_image_pixels_downloaded") is not True or
        prior.get("no_1to1_matching_or_record_discard_for_yearmonth10km") is not True):
        raise ValueError("Original full Czech 4,235 source photo metadata v17.5 genuine audit missing")
    keys=[
       "photo_original_source_privacy_taxon_and_GPS_quality",
       "all_public_observations_annotation_outcomes",
       "original_railway_polyline_distance_zone_anthropogenic_observation_metadata",
       "source_model_identifiable_for_two_photo_distance_groups"]
    if any(current.get(k)!=prior.get(k) for k in keys):
        raise ValueError("Original real source photo and official Flowers/No Flowers annotation counts drifted from v17.5")
    if current["photo_original_source_privacy_taxon_and_GPS_quality"].get("pass")!=4176:
        raise ValueError("Original 4,176 source public geo/taxon/research-photo observations changed")
    return True

def create_queue(obs,rail,original_v175):
    if len(obs)!=EXPECTED_TOTAL:raise ValueError("Source photo API 4235 records incomplete")
    review=[];context=[];seen=set()
    quality=Counter()
    for ob in obs:
        meta,q=public_metadata_gate(ob)
        quality[q]+=1
        if meta is None:continue
        observation_id=ob.get("id")
        key=token_id(observation_id)
        if key in seen:raise ValueError("Repeated source public iNaturalist observation ID review key")
        seen.add(key)
        dist,cell=rail.distance_and_cell(meta["lon"],meta["lat"])
        zone=distance_zone(dist)
        calendar=str(ob["observed_on"])[:10]
        year=int(calendar[:4])
        # Prefer exact source observation date over invented date from UI.
        import datetime
        date=datetime.date.fromisoformat(calendar)
        if not(2018<=date.year<=2025):raise ValueError("Source original photo date changed")
        public_observer=(ob.get("user") or {}).get("id")
        observer_group=(hashlib.sha256(
          ("tetsu-public-observer-group-v176-"+str(public_observer)).encode("utf-8")).hexdigest()[:20]
          if public_observer is not None else "")
        accuracy=meta["accuracy_m"]
        precise=meta["accuracy_class"]=="accuracy_1to30m"
        review.append({
          "public_observation_id":observation_id,
          "review_token":key,
          "public_iNaturalist_observation_page":f"https://www.inaturalist.org/observations/{observation_id}",
          "human_REVIEW_open_flower_visible_yes_no_uncertain":"",
          "human_REVIEW_species_visually_plausible_yes_no_uncertain":"",
          "human_REVIEW_whole_plant_visible_yes_no_uncertain":"",
          "human_REVIEW_plant_scale_measurable_yes_no_uncertain":"",
          "human_REVIEW_image_target_unclear_or_multiple_plants":"",
          "human_REVIEW_notes":"",
        })
        context.append({
          "review_token":key,
          "source_species":TARGET,
          "source_original_public_photo_date":calendar,
          "source_original_observed_YEAR":year,
          "source_original_DOY":date.timetuple().tm_yday,
          "source_original_calendar_MONTH":date.month,
          "source_public_coord_accuracy_class":meta["accuracy_class"],
          "source_public_coord_accuracy_m":accuracy if accuracy is not None else "",
          "source_nearest_real_2023_OSM_rail_LINE_distance_m":round(dist,1),
          "source_rail_LINE_distance_bin":zone,
          "source_public_error_circle_within30m":bool(
              accuracy is not None and dist+accuracy<=30),
          "source_PUBLIC_observer_pseudonym_for_cluster_adjustment":observer_group,
          "source_10km_COOARSENED_location_grid_not_true_GPS":
              f"{math.floor(cell[0]/100)}_{math.floor(cell[1]/100)}",
          "source_existing_iNat_photo_annotation_EXPLICIT_only":annotation_class(ob),
          "source_original_photo_license_NOT_USED_FOR_METADATA":True,
          "original_2023_OSM_line_not_actual_management_or_legal_railway":True,
        })
    if len(review)!=4176 or len(context)!=4176 or len(seen)!=4176:
        raise ValueError("Original public accepted 4,176 observations not conserved in human QA queue")
    # Always make outcome-blind wave, and keep ALL original review rows
    # available without throwing away/1:1 matching the remaining sample.
    entries=list(zip(review,context))
    near=[(r,c) for r,c in entries if (
       c["source_public_coord_accuracy_class"]=="accuracy_1to30m" and
       c["source_rail_LINE_distance_bin"]=="point_0to30m")]
    far=[(r,c) for r,c in entries if (
       c["source_public_coord_accuracy_class"]=="accuracy_1to30m" and
       c["source_rail_LINE_distance_bin"]=="point_100to1000m")]
    if len(near)!=178 or len(far)!=725:
        raise ValueError("Source original v17.5 full strict 178 near/725 control records not recovered")
    # First-wave cap sampled from controls by hashed original id only,
    # NOT its photos, any annotation field, observation date, or original
    # geographic distance beyond being in declared exposure category.
    far.sort(key=lambda z:hashlib.sha256(
        ("tetsu-public-qa-control-firstwave-v176-"+z[0]["review_token"]).encode("utf-8")).hexdigest())
    chosen=[r for r,_ in near]+[r for r,_ in far[:FIRST_WAVE_MAX_CONTROL]]
    chosen.sort(key=lambda r:r["review_token"])
    summary={
      "complete":True,"version":"v17.6_Echium_Czech_true_public_observation_ID_NO_pixels_label_REVIEW_queue",
      "source_2018_2025_original_research_grade_photo_API_count":EXPECTED_TOTAL,
      "original_source_v175_exact_phenophase_annotation_comparison_reproduced":True,
      "n_non_private_observation_IDs_with_public_original_photo_date_and_geometry":len(review),
      "n_available_original_near30_error30_source_records":len(near),
      "n_available_original_control_100to1000_error30_source_records":len(far),
      "n_FIRSTWAVE_source_observations":len(chosen),
      "n_FIRSTWAVE_source_near30_original_observations":len(near),
      "n_FIRSTWAVE_source_control_100to1000_observations":len(chosen)-len(near),
      "source_firstwave_selection_USED_original_source_observation_ID_hash_NOT_photo_annotation_or_phenotype":True,
      "no_1to1_matching_and_all_other_source_records_preserved_for_regression":True,
      "photo_pixels_NEVER_downloaded_or_redistributed":True,
      "only_originally_PUBLIC_non_obscured_iNaturalist_observation_IDs_in_EPHEMERAL_artifact_NOT_git_repo":True,
      "links_to_original_PUBLIC_observation_pages_can_reveal_originally_PUBLIC_localities":True,
      "source_iNaturalist_observer_id_NOT_exported_only_group_hash":True,
      "precise_public_coordinates_NOT_exported":True,
      "review_queue_OMITS_original_distance_and_iNat_phenology_annotations_BUT_original_iNat_pages_may_show_map_and_tags":True,
      "source_REVIEW_fields_left_all_blank":True,
      "not_fully_blinded_human_evaluation_by_design":True,
      "not_direct_measurement_of_flowering_onset_plant_height_or_local_genetic_adaptation":True,
    }
    return review,context,chosen,summary

def export_csv(path,records):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not records:raise ValueError("No source public photo records eligible for review")
    with path.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def run(source_v166,source_v175,rail,client=None):
    client=client or ConservativeINat(delay=1.1)
    response,obs=load_complete_photo_observations(client,TARGET,STUDY_REGIONS["CZ"]["bbox"])
    if not response["complete"] or response["total"]!=EXPECTED_TOTAL or obs is None:
        raise ValueError("Public original all photo source full observation API unavailable")
    p=next(x for x in source_v166["source_species_country_counts"]
           if x["region_source_bbox"]=="CZ" and x["taxon"]==TARGET)
    current=process_records(obs,rail,p)
    verify_source(current,source_v175)
    return create_queue(obs,rail,source_v175)

def main():
    p=argparse.ArgumentParser()
    for arg in ("source-v166","source-v175","original-Czech-OSM","output-dir"):
        p.add_argument("--"+arg,required=True)
    a=p.parse_args()
    v166=json.loads(Path(a.source_v166).read_text(encoding="utf-8"))
    v175=json.loads(Path(a.source_v175).read_text(encoding="utf-8"))
    original_pbf_authentic(a.original_Czech_OSM,"CZ")
    rail=RailwayPolylineSource(a.original_Czech_OSM,"CZ")
    review,context,chosen,summary=run(v166,v175,rail)
    dest=Path(a.output_dir);dest.mkdir(parents=True,exist_ok=True)
    export_csv(dest/"original_public_observation_REVIEW_FIELDS_BLANK.csv",review)
    export_csv(dest/"source_original_2023_rail_LINE_context_NO_coords.csv",context)
    export_csv(dest/"FIRST_WAVE_356_public_obs_REVIEW_FIELDS_BLANK.csv",chosen)
    (dest/"source_qa_counts.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False))

if __name__=="__main__":main()
