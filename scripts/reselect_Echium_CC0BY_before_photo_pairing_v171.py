#!/usr/bin/env python3
"""v17.1 RIGHT-AWARE matched Echium original 2023 Czech railway photo review.

v17.0 established 29 original iNaturalist true-source photographic near30
vs rail100-1000m 1:1 pairs matched by year/month, within 10km, separate
observers, BEFORE filtering by CC0/CC BY derivative distribution rights.
Only ONE of those 29 pairs had BOTH CC0/CC BY attributed photos.

This independent v17.1 feasibility SAMPLE selects CC0/CC BY licensed photo
records FIRST, then uses the SAME FROZEN 1:1 identical study year/month,
geography within10km, 30m coordinate-error circle inside original 2023
physical OSM railroad LINE and distinct observers to rematch. This may
recover additional rights-cleared photo-pairs while changing WHICH ORIGINAL
iNaturalist observations are paired. This is an exploratory RESELECTION
cohort for image REVIEW, not the frozen v16.8 29 matched records.

No true flowering onset/plant height from a single photo; no selected taxon
railway adaptation claim. Images posted only to ephemeral GitHub Actions
artifact, licensed CC0/BY with original photographer attribution, no
original observation IDs, user IDs, private GPS or exact photo times exported.

A/B side orders are blinded independent of licensing or floral appearance.
Reviewer blank CSV and unblind key separated as in v17.0. No automated
semantic species or flower identification, no real measuring scale inference.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from collections import Counter
from pathlib import Path
import requests
from PIL import Image

from audit_20_original_rail_plants_inat_30m_two_regions_v165 import (
    RailwayPolylineSource, ConservativeINat, STUDY_REGIONS,
    original_pbf_authentic, load_complete_photo_observations
)
from audit_Czech_four_flower_taxa_local_near30_vs_control_photo_pairs_v168 import (
    _private_real_record, one_to_one_matched, summarize_source_species
)
from build_Echium_29_railphoto_blinded_CC_contact_sheets_v170 import (
    TARGET,EXPECTED_REAL_PAIRS,EXPECTED_ORIGINAL_API_TOTAL,
    source_gate,eligible_public_photo,blind_order,compose_blind_pair,
    safe_actual_image_download
)

V171="v17.1_Echium_RIGHTS_FIRST_Czech2023_rail_real_photo_REmatched_CC0BY"

def licensed_first_matching(original_obs,rail,prior,reference):
    if len(original_obs)!=EXPECTED_ORIGINAL_API_TOTAL:
        raise ValueError("Original 4,235 Czech iNat all research-grade photo cohort changed")
    before=summarize_source_species(original_obs,rail,TARGET,prior)
    if before["n_same_year_month_unique_photo_pairs_within10km_different_observers"]!=EXPECTED_REAL_PAIRS:
        raise ValueError("Original v16.8 29 source pairs no longer verified")
    if before["n_same_year_month_unique_photo_pairs_within10km_different_observers"]!=reference["n_same_year_month_unique_photo_pairs_within10km_different_observers"]:
        raise ValueError("Original v16.8 source matched pair result changed")
    near=[];controls=[];orig={}
    qc=Counter()
    for ob in original_obs:
        record,reason=_private_real_record(ob,rail,TARGET)
        if record is None:
            qc["not_in_original_true_location_or_allCC_photo_eligible"]+=1
            continue
        info=eligible_public_photo(ob)
        if info is None:
            qc["source_in_band_but_no_CC0BY_image_or_attribution"]+=1
            continue
        orig[record["internal_original_id"]]=info
        if record["group"]=="near30_certain":
            near.append(record)
        elif record["group"]=="control_100to1000":
            controls.append(record)
    matched=one_to_one_matched(near,controls,restrict_different_observers=True)
    matched_any_observer=one_to_one_matched(near,controls,restrict_different_observers=False)
    result=[(n,c,orig[n["internal_original_id"]],orig[c["internal_original_id"]])
            for n,c in matched]
    qc.update({
       "n_source_original_v168_allCC_pairs_reproduced":29,
       "n_source_original_near30_CC0_BY_observations":len(near),
       "n_source_original_rail100to1000_CC0_BY_observations":len(controls),
       "n_source_rights_first_CC0_BY_same_yearmonth_10km_different_observer_pairs":len(result),
       "n_source_rights_first_CC0_BY_any_observer_pairs":len(matched_any_observer)
    })
    return result,dict(qc)

def run(v168,v166,rail,output,client=None):
    reference,prior=source_gate(v168,v166)
    client=client or ConservativeINat(delay=1.10)
    api,records=load_complete_photo_observations(client,TARGET,STUDY_REGIONS["CZ"]["bbox"])
    if not api["complete"] or records is None:
        raise ValueError("iNaturalist 4,235 original source records incompletely enumerated")
    pairs,qc=licensed_first_matching(records,rail,prior,reference)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    attribution=[];templates=[];unblind=[]
    failures=0;reviewed=0
    session=requests.Session()
    session.headers.update({"User-Agent":"tetsu-CC0BY-original-Czech-Echium-photo-REmatching-17.1 (github.com/zuizui0223/tetsu)"})
    for n,c,ninfo,cinfo in pairs:
        images={};details={}
        try:
            for role,info in (("near",ninfo),("control",cinfo)):
                time.sleep(1.10)
                image,metrics=safe_actual_image_download(session,info)
                images[role]=image;details[role]=metrics
        except (requests.RequestException,OSError,ValueError):
            failures+=1
            continue
        reviewed+=1
        order=blind_order(n,c)
        file=f"echium_rightsfirst_blind_pair_{reviewed:02d}.jpg"
        compose_blind_pair(images[order[0]],images[order[1]],reviewed).save(output/file,quality=87)
        unblind.append({
            "pair_number":reviewed,"A_source_group":order[0],"B_source_group":order[1],
            "original_source_year_month_MATCHED":True,
            "original_source_geographic_distance_up_to10km":True
        })
        for side,role in (("A",order[0]),("B",order[1])):
            info=ninfo if role=="near" else cinfo
            attribution.append({
                "pair_number":reviewed,"blinded_side":side,
                "source_photo_url_for_required_attribution":info["source_url"],
                "source_photo_rights_CC0_or_CC_BY":info["license_code"],
                "source_original_author_photo_attribution":info["attribution"],
                **details[role]
            })
            templates.append({
                "pair_number":reviewed,"blinded_side":side,"blind_image_file":file,
                "target_Echium_visually_plausible_yes_no_uncertain":"",
                "visible_open_flower_yes_no_uncertain":"",
                "whole_single_plant_yes_no_uncertain":"",
                "base_visible_yes_no_uncertain":"",
                "apex_visible_yes_no_uncertain":"",
                "reliable_absolute_size_calibration_visible_yes_no":"",
                "railway_trackbed_or_rail_actually_visible_yes_no_uncertain":"",
                "relative_apparent_silhouette_reliably_comparable_yes_no":"",
                "reviewer_blinded_comments":""
            })
    def write_csv(filename,records):
        if not records:return
        with (output/filename).open("w",encoding="utf-8",newline="") as fp:
            w=csv.DictWriter(fp,fieldnames=list(records[0]))
            w.writeheader();w.writerows(records)
    write_csv("original_photo_image_attribution.csv",attribution)
    write_csv("blank_blinded_photo_human_QA.csv",templates)
    write_csv("UNBLIND_ONLY_AFTER_rating.csv",unblind)
    summary={
        "complete":True,"version":V171,
        "source_inat_original_Echium_Czech_full_research_grade_photo_observations":api["total"],
        "source_v168_fixed_29_pairs_reproduced_BEFORE_rights_filter":EXPECTED_REAL_PAIRS,
        "source_after_first_CC0_BY_then_rematch_study_cohort_DIFFERENT_from_v168":True,
        "source_CC0_BY_rights_first_matching_counts":qc,
        "n_CC0BY_already_licensed_both_sides_photo_pairs_matched_after_rights_select":len(pairs),
        "n_actual_both_sides_photographic_image_bytes_decoded":2*reviewed,
        "n_actual_blinded_photo_contact_pairs_for_HUMAN_QA":reviewed,
        "n_source_image_download_or_decode_failure_pairs":failures,
        "actual_image_pixels_have_been_retrieved_from_rights_approved_public_sources":bool(reviewed),
        "photo_image_content_OPEN_FLOWER_or_whole_PLANT_NOT_yet_blind_reviewed":True,
        "all_blind_reviewer_phenotype_ratings_are_blank":True,
        "source_CC0BY_attributions_and_original_public_source_links_in_separate_file":True,
        "precise_location_observationIDs_and_userIDs_NOT_exposed":True,
        "no_plant_height_flowering_onset_or_rail_specific_adaptation_inferred":True,
        "method_warning":"Matches source CC0/BY photo records BEFORE pairing, unlike original v16.8's 29 CC-all record-pairs, and could select different iNaturalist original photos; source opportunity, not an independent ecological or biological effect."
    }
    (output/"study_source_real_pixels_result.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser()
    for k in ("source-v168","source-v166","source-Czech-2023-OSM","output-dir"):
        p.add_argument("--"+k,required=True)
    a=p.parse_args()
    v168=json.loads(Path(a.source_v168).read_text(encoding="utf-8"))
    v166=json.loads(Path(a.source_v166).read_text(encoding="utf-8"))
    original_pbf_authentic(a.source_Czech_2023_OSM,"CZ")
    rail=RailwayPolylineSource(a.source_Czech_2023_OSM,"CZ")
    run(v168,v166,rail,a.output_dir)

if __name__=="__main__":main()
