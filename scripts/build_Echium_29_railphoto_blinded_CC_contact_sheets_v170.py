#!/usr/bin/env python3
"""v17.0 ACTUAL iNaturalist CC0/CC-BY photograph pixel and blinded paired
contact-sheet review cohort for Czech original railway-flora Echium vulgare.

It re-queries real public iNaturalist research-grade photographed botanical
observations (2018–2025) and exact historical 2023 Czech Geofabrik 816MB OSM
railway LINE geometry. It reproduces original source-v16.8 *29* matched
same-year/month, within-10km, different-observer near vs far image metadata
pairs BEFORE filtering by redistributable image license/actual HTTP pixels.

Only CC0 or CC-BY independently photo-licensed iNaturalist image URLs
and photo attribution are placed in GitHub Actions artifacts, NOT git source.
No image ever fetched from an arbitrary API URL host. Contact sheets have
blind left/right A/B labels; a separate sidecar unblinds source near/control
after independent human review. No original obs IDs, user IDs, precise
coordinates, exact spatial cells, photographs outside the license
allowlist or unlicensed URL are exported. Photo attribution and original
public photo source URLs must be preserved with the binaries.

Automated quality is limited to successful REAL image GET, JPEG decode,
dimensions, image aspect/orientation and repeat-byte hash. It does not
identify flowers, plants, railway tracks, full height, a camera scale or
measured phenology; optional pre-existing iNat Flowers tags are not
positive/negative visual ground-truth. All actual botanical assessments
in reviewer_blank.csv remain intentionally unfilled.
"""
from __future__ import annotations

import argparse, csv, hashlib, io, json, math, time
from collections import Counter
from pathlib import Path

import requests
from PIL import Image, ImageDraw, ImageOps

from audit_20_original_rail_plants_inat_30m_two_regions_v165 import (
    RailwayPolylineSource, ConservativeINat, STUDY_REGIONS,
    original_pbf_authentic, load_complete_photo_observations,
)
from audit_Czech_four_flower_taxa_local_near30_vs_control_photo_pairs_v168 import (
    _private_real_record, one_to_one_matched, summarize_source_species,
)
from sample_photo_roi_candidates import allowed_photo_url, image_fetch_url

SOURCE_V168="v16.8_GENUINE_Czech_four_angiosperms_local_CC_photo_PAIRS_review_opportunity"
SOURCE_V166="v16.6_AUTHENTIC_2COUNTRY_20_ORIGINAL_RAIL_TAXA_PUBLIC_PHOTO_30m_FEASIBILITY_QA"
TARGET="Echium vulgare"
EXPECTED_REAL_PAIRS=29
EXPECTED_ORIGINAL_API_TOTAL=4235
ALLOWED_LICENSES={"cc0","cc-by"}
CONTACT_SIZE=(450,450)
MAX_IMAGE_BYTES=5_000_000
MAX_IMAGE_SIDE=1400

def source_gate(v168,v166):
    if not(v168.get("complete") is True and v168.get("version")==SOURCE_V168
           and v168.get("source_Geofabrik_original_Czech2023_PBF_SHA")==
                  "26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"
           and v168.get("NO_raw_original_photos_observation_ID_public_coordinates_or_usernames_exported") is True):
        raise ValueError("Source verified four Czech original real CC photo matched study v16.8 lost")
    if v166.get("complete") is not True or v166.get("version")!=SOURCE_V166:
        raise ValueError("Original v16.6 full iNaturalist photo/genuine historic OSM source audit absent")
    original=[v for v in v168["source_restricted_original_4_taxa_matched_opportunity"]
              if v["source_species"]==TARGET]
    if len(original)!=1:
        raise ValueError("Original v16.8 Echium genuine published-flora target not uniquely present")
    reference=original[0]
    if (reference["n_same_year_month_unique_photo_pairs_within10km_different_observers"]!=EXPECTED_REAL_PAIRS
        or reference["source_2018_2025_full_real_API_photograph_observation_count"]!=EXPECTED_ORIGINAL_API_TOTAL):
        raise ValueError("Original v16.8 29 matched genuine Echium pairs/4,235 country photo results drifted")
    prior=[v for v in v166["source_species_country_counts"]
           if v["region_source_bbox"]=="CZ" and v["taxon"]==TARGET]
    if len(prior)!=1 or prior[0]["source_bbox_iNat_research_photo_total"]!=EXPECTED_ORIGINAL_API_TOTAL:
        raise ValueError("Original source 2023 CZ study Echium full iNat original photo count drift")
    return reference,prior[0]

def eligible_public_photo(ob):
    """Return FIRST source-order matched image that is explicitly CC0 or CC BY.

    All other CC licenses were valid for v16.8 *metadata* counts, but this
    separate real-image republishing review purpose strictly limits to
    source CC0 / CC-BY and preserves author attribution.
    """
    for p in ob.get("photos") or ():
        license_code=str(p.get("license_code") or "").lower()
        url=p.get("url") or ""
        if license_code not in ALLOWED_LICENSES:
            continue
        if not allowed_photo_url(url):
            continue
        original_source_url=image_fetch_url(p)
        if not allowed_photo_url(original_source_url):
            continue
        # CC BY author metadata attribution is mandatory for image copies.
        attribution=str(p.get("attribution") or "").strip()
        if license_code=="cc-by" and not attribution:
            continue
        return {
            "source_url":original_source_url,
            "attribution":attribution if attribution else "CC0",
            "license_code":license_code,
        }
    return None

def blind_order(n_record,c_record):
    """A/B ordering depends only on deterministic hash of private IDs,
    NEVER the species, image flower tag, distance, year or manual QA.
    """
    h=hashlib.sha256(
       ("tetsu-v170-blind-contact-sheet-"+str(n_record["internal_original_id"])+
        "-"+str(c_record["internal_original_id"])).encode("utf-8")).digest()
    return ("near","control") if (h[0] & 1)==0 else ("control","near")

def choose_review_pairs(original_obs,rail,reference,prior):
    """Reproduce original 29 source paired PHOTO opportunities exactly,
    then restrict to rights-cleared CC0/CC-BY **pixel review** subset.
    """
    if len(original_obs)!=EXPECTED_ORIGINAL_API_TOTAL:
        raise ValueError("Original source iNat API full Echium photo census drift; never silently replace 29")
    original=summarize_source_species(original_obs,rail,TARGET,prior)
    if original["n_same_year_month_unique_photo_pairs_within10km_different_observers"]!=EXPECTED_REAL_PAIRS:
        raise ValueError("Original real iNat source matched photo cohort not exactly original 29")
    if (original["n_source_point_near30_and_errorcircle_certain_and_CC_licensed"]!=
        reference["n_source_point_near30_and_errorcircle_certain_and_CC_licensed"]):
        raise ValueError("Original source error circle and CC geometry photo support drift")
    near=[];far=[]
    idmap={}
    for ob in original_obs:
        rec,why=_private_real_record(ob,rail,TARGET)
        if rec is None:
            continue
        idmap[rec["internal_original_id"]]=ob
        if rec["group"]=="near30_certain":near.append(rec)
        elif rec["group"]=="control_100to1000":far.append(rec)
    matches=one_to_one_matched(near,far,restrict_different_observers=True)
    if len(matches)!=EXPECTED_REAL_PAIRS:
        raise ValueError("Real original v16.8 1:1 near vs control source matches changed")
    qc=Counter()
    licensed=[]
    for n,c in matches:
        nphoto=eligible_public_photo(idmap[n["internal_original_id"]])
        cphoto=eligible_public_photo(idmap[c["internal_original_id"]])
        if not nphoto or not cphoto:
            qc["pair_missing_CC0_CC_BY_original_public_photo_or_photo_attribution"]+=1
            continue
        licensed.append((n,c,nphoto,cphoto))
    qc["n_original_source_photo_pairs"]=len(matches)
    qc["n_public_both_sides_CC0_CC_BY_attributed_photo_pairs"]=len(licensed)
    return licensed,dict(qc)

def safe_actual_image_download(session,info):
    """No untrusted photo host and capped bytes; actual bytes decoded only
    from rights-cleared public iNaturalist static image hosts."""
    u=info["source_url"]
    if not allowed_photo_url(u):
        raise ValueError("Blocked non-iNaturalist photo download host")
    with session.get(u,stream=True,timeout=(12,38),allow_redirects=True) as r:
        r.raise_for_status()
        if not allowed_photo_url(r.url):
            raise ValueError("Redirect out of allowed original iNaturalist image hosts")
        chunks=[];size=0
        for block in r.iter_content(64*1024):
            size+=len(block)
            if size>MAX_IMAGE_BYTES:
                raise ValueError("Original source image exceeds 5MB cap")
            chunks.append(block)
    raw=b"".join(chunks)
    with Image.open(io.BytesIO(raw)) as photo:
        photo=ImageOps.exif_transpose(photo).convert("RGB")
        full_size=photo.size
        photo.thumbnail((MAX_IMAGE_SIDE,MAX_IMAGE_SIDE),Image.Resampling.LANCZOS)
        image=photo.copy()
    return image,{
        "source_original_public_image_file_sha256":hashlib.sha256(raw).hexdigest(),
        "source_original_image_dimensions_px":list(full_size),
        "source_preview_decoded_dimensions_px":list(image.size),
        "source_image_HTTPS_decoded":True
    }

def compose_blind_pair(A,B,index):
    width,height=CONTACT_SIZE
    sheet=Image.new("RGB",(2*width,height+80),"white")
    draw=ImageDraw.Draw(sheet)
    for side,image in (("A",A),("B",B)):
        img=image.copy()
        img.thumbnail((width-20,height-22),Image.Resampling.LANCZOS)
        xbase=0 if side=="A" else width
        px=xbase+(width-img.width)//2
        py=(height-img.height)//2
        sheet.paste(img,(px,py))
        draw.rectangle((xbase,0,xbase+width-1,height),outline="gray")
        draw.text((xbase+12,height+12),f"{side}  (blinded source role)",fill="black")
        draw.text((xbase+12,height+32),"CC0 / CC BY author credit: attribution.csv",fill="black")
    draw.text((8,height+56),f"Echium vulgare  matched image review sample {index:02d}; not a phenotype label",fill="black")
    return sheet

def run(v168,v166,rail,output,client=None,max_pairs=50):
    reference,prior=source_gate(v168,v166)
    if not 1<=max_pairs<=50:
        raise ValueError("Invalid number of source pair review slots")
    client=client or ConservativeINat(delay=1.10)
    api,original_obs=load_complete_photo_observations(
        client,TARGET,STUDY_REGIONS["CZ"]["bbox"])
    if not api["complete"] or original_obs is None:
        raise ValueError("Original 4,235 Czech Echium iNaturalist records not fully enumerated")
    pairs,qc=choose_review_pairs(original_obs,rail,reference,prior)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    template=[]
    provenance=[]
    unblind=[]
    notes=[]
    valid=0
    sess=requests.Session()
    sess.headers.update({"User-Agent":"tetsu-original-Czech-rail-photo-noncausal-review/17.0 (github.com/zuizui0223/tetsu)"})
    pairs=pairs[:max_pairs]
    for n,c,ninfo,cinfo in pairs:
        images={}
        measured={}
        try:
            for role,info in (("near",ninfo),("control",cinfo)):
                # Limit request frequency; no programmatic bulk/reuse of
                # original image pixels outside this rights-audited contact
                # sheet sample.
                time.sleep(1.10)
                img,met=safe_actual_image_download(sess,info)
                images[role]=img;measured[role]=met
        except (requests.RequestException,OSError,ValueError) as exc:
            qc["original_source_photo_HTTPS_fetch_or_decode_failed"]+=1
            continue
        order=blind_order(n,c)
        valid+=1
        filename=f"echium_blind_pair_{valid:02d}.jpg"
        A=images[order[0]];B=images[order[1]]
        sheet=compose_blind_pair(A,B,valid)
        sheet.save(output/filename,quality=87)
        unblind.append({"pair_number":valid,"blinded_A_source_group":order[0],
                        "blinded_B_source_group":order[1],
                        "pair_same_calendar_year_month":True,
                        "source_nearest_true_rail_line_and_public_accuracy_guaranteed30":True,
                        "control_original_OSM_rail_line_distance_100to1000m":True})
        for side,role in (("A",order[0]),("B",order[1])):
            info=ninfo if role=="near" else cinfo
            meta=measured[role]
            # The attributions reference the original public photo asset
            # and its human author; no original observation ID, GPS, full
            # user account or exact photo time is saved here.
            provenance.append({
               "pair_number":valid,"side_A_or_B":side,
               "original_photo_source_public_URL":info["source_url"],
               "original_photo_license":info["license_code"],
               "original_photo_attribution_credit":info["attribution"],
               **meta
            })
            template.append({
               "pair_number":valid,"blinded_side":side,"blind_pair_image":filename,
               "source_photo_appears_to_show_correct_Echium_species":"",
               "open_flower_visible_yes_no_uncertain":"",
               "whole_single_plant_visible_yes_no_uncertain":"",
               "plant_base_visible_yes_no_uncertain":"",
               "plant_top_visible_yes_no_uncertain":"",
               "single_target_or_many_plants":"",
               "unambiguous_physical_scale_in_image_yes_no":"",
               "rail_track_or_ballast_ACTUALLY_visible_yes_no_uncertain":"",
               "photograph_suitable_for_same_view_relative_morphology_yes_no":"",
               "freeform_review_notes":""
            })
    def csvwrite(name,records):
        if not records:return
        with (output/name).open("w",encoding="utf-8",newline="") as fp:
            writer=csv.DictWriter(fp,fieldnames=list(records[0]))
            writer.writeheader();writer.writerows(records)
    csvwrite("image_attribution_and_verified_bytes.csv",provenance)
    csvwrite("blinded_image_human_QA_BLANK_template.csv",template)
    csvwrite("UNBLIND_AFTER_human_QA.csv",unblind)
    stats={
       "complete":True,
       "version":"v17.0_Echium_GENUINE_CC0_CC_BY_photo_pixels_blinded_PAIRED_review_cohort",
       "source_original_study":"Kutlvasr_2024 original Czech railway station flora; candidate was preselected as one of 20 source railway taxa",
       "source_original_real_iNat_2018to2025_CZ_bbox_photo_records":api["total"],
       "original_v168_within10km_sameYEAR_MONTH_different_observer_real_photo_pairs":EXPECTED_REAL_PAIRS,
       "n_both_sides_CC0_or_CC_BY_with_attribution":qc["n_public_both_sides_CC0_CC_BY_attributed_photo_pairs"],
       "n_licensed_paired_photo_assets_successfully_retrieved_and_decoded":valid*2,
       "n_licensed_blinded_true_source_image_pairs_for_HUMAN_REVIEW":valid,
       "original_CC_image_lisense_subset_NOT_all_previous_v168_CC_types":True,
       "source_license_rejection_and_photo_decode_QA":qc,
       "original_actual_image_pixels_downloaded_and_decoded":True if valid else False,
       "original_plant_part_semantic_visual_classifications_COMPLETED":False,
       "original_photo_blinded_review_annotations_are_ALL_BLANK":True,
       "original_CC_image_credit_required_and_in_sidecar":True,
       "actual_image_sheets_in_GitHub_Actions_ARTIFACT_ONLY_not_committed_to_source_repo":True,
       "photo_or_observation_precise_lat_lon_and_original_observation_ID_NOT_EXPORTED":True,
       "no_rail_right_of_way_or_visible_flower_verification":True,
       "no_phenological_first_flowering_or_physical_height_inference":True,
       "notes":[
         "Thumbnails are actual iNaturalist CC0/CC BY photos, independent photographers credited in file sidecar; original full image URLs included for legal attribution.",
         "Photos are blinded A/B inside each genuine matched pair. Source role unblinding is held in separate file, to be opened after manual review.",
         "No pretrained image recognizer or human reviewed the actual plant content; image readability is NOT confirmed merely from downloaded pixels.",
         "Original 29 pairs included broader Creative Commons licenses, of which only both-side CC0/CC BY can be placed in this redistribution artifact.",
         "Separate true iNaturalist observation identities, original coordinates and observer user IDs were used transiently for 1:1 matching but are never exported.",
         "No species-internal measurable trait has yet been demonstrated; original site-based incidence cannot support exact true rail/road CWM without matched flora.",
         "Never upload trackside precise sensitive observation geolocations, and do not collect additional pictures by trespassing on tracks."
       ]
    }
    (output/"source_review_counts.json").write_text(json.dumps(stats,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(stats,ensure_ascii=False),flush=True)
    return stats

def main():
    p=argparse.ArgumentParser()
    for k in ("source-v168","source-v166","original-Czech-2023-osm","out"):p.add_argument("--"+k,required=True)
    p.add_argument("--max-pairs",type=int,default=29)
    a=p.parse_args()
    v168=json.loads(Path(a.source_v168).read_text(encoding="utf-8"))
    v166=json.loads(Path(a.source_v166).read_text(encoding="utf-8"))
    original_pbf_authentic(a.original_Czech_2023_osm,"CZ")
    rail=RailwayPolylineSource(a.original_Czech_2023_osm,"CZ")
    run(v168,v166,rail,a.out,max_pairs=a.max_pairs)

if __name__=="__main__":main()
