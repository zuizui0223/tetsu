#!/usr/bin/env python3
"""v16.7: from two AUTHOR-SOURCE-VERIFIED v16.6 photo metadata artifacts,
quantify usable flowering-annotation and seasonal photo REVIEW opportunities.
Not a phenology analysis or an analysis of full-plant size.

v16.5 iNaturalist source observations 2018-2025 around authentic 2023/2024
country OSM railway polylines were fully enumerated, non-obscured,
non-captive, research-grade and with public coordinate accuracy <=30m.
v16.6 independently checked author-original exact 20 species, actual
historic Geofabrik map SHA, public photo quality, observer and 100m cell
counts, and fixed baseline eligibility. 11/40 30m botanical species-country
cells passed the basic photo-opportunity gate, all in Czech bbox.

This v16.7 step uses only source v16.6 audited metadata.
- Never confuse "Flowers" annotation, if present, with first flowering date,
  non-flowering negative phenotype, herbivory, selection or fitness.
- A missing photo annotation means NO annotation, NOT "no flowers".
- Historical photos have neither known image scale nor original local
  mowing status, train-frequency, wind stress or individual height.
- Equisetum arvense is a sporophyte with fertile strobilus, not flowers:
  it has no true angiosperm flowering response.
- Globally counted iNat research photos are NOT globally rail-proximate.
- Source is exact two disjoint country BBOXES, not worldwide survey.
- Do NOT infer adaptation or community-weighted trait changes from photos.

Exploratory stage-2 IMAGE REVIEW opportunity thresholds, predeclared here
before seeing source v16.7 flower-tag/month results (NOT predeclared
before v16.5 30m photo abundance):
  v16.6 basic photo gate; ≥20 CC photo observations near30m;
  ≥20 within30m with public accuracy circle fully inside the line buffer;
  ≥3 distinct observation calendar months with ≥5 near-track records each;
  ≥5 explicit "Flowers" annotated near-track photos for OPTIONAL label QC.
Those are review logistics only, NOT an empirical flowering-shift test.
We can reduce candidate list based on review cost; never select positive
railway effect size post hoc.
"""
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path

V166="v16.6_AUTHENTIC_2COUNTRY_20_ORIGINAL_RAIL_TAXA_PUBLIC_PHOTO_30m_FEASIBILITY_QA"
REGIONS=("SI","CZ")
MIN_CC=20
MIN_GUARANTEED30=20
MIN_MONTHLY_N=5
MIN_MONTHS=3
MIN_ANNOTATED_FLOWER=5
ANGIOSPERM_REVIEW_ONLY=True
FLOWERLESS_TAXA={"Equisetum arvense"}
ORIGINAL_BOTANICAL=(
  "10.1007/s44473-026-00155-z","10.23855/preslia.2024.247"
)

def verified_rows(source):
    if source.get("complete") is not True or source.get("version")!=V166:
        raise ValueError("The original 2023/2024 Geofabrik+2018-25 iNat photo study source v16.6 is not authenticated")
    if not(source.get("not_an_accurate_GLOBAL_near_railway_photo_count") is True
           and source.get("not_an_actual_measurement_of_photo_plant_height_or_first_flowering") is True):
        raise ValueError("Original global-rail/photo-height/phenophase source limitation lost")
    rows=source.get("source_species_country_counts")
    if not isinstance(rows,list) or len(rows)!=40:
        raise ValueError("Two original source-country region x 20 different plant species expected")
    keys={(r["region_source_bbox"],r["taxon"]) for r in rows}
    if len(keys)!=40 or {k[0] for k in keys}!=set(REGIONS):
        raise ValueError("A source taxon/region is missing, duplicated, or an unverified new country")
    if source["n_original_taxa"]!=20 or source["n_original_source_country_bbox_taxon_complete"]!=40:
        raise ValueError("Original source 20 unselected species x 2 country bbox drift")
    if sum(r["n_rail_polyline_point_near_30m"] for r in rows)!=source["n_source_photo_rail_point_within30m_both_bboxes"]:
        raise ValueError("Source actual near30m botanist photo count not conserved")
    if sum(r["n_rail_polyline_guaranteed_circle_within_30m"] for r in rows)!=source["n_location_within_30m_even_with_coordinate_uncertainty_circle"]:
        raise ValueError("Original 30m point coordinate error circle support changed")
    if sum(r["n_rail_polyline_30m_with_CC_photo_license"] for r in rows)!=source["n_photo_CC_licensed_within30m_both_bboxes"]:
        raise ValueError("Original source CC photo count not conserved")
    if sum(bool(r["minimum_original_photo_feasibility_gate_passed"]) for r in rows)!=source["n_country_taxon_cells_passing_predeclared_minimum_photographic_opportunity"]:
        raise ValueError("Original v16.6 author-source basic eligibility gate changed")
    return rows

def evaluate(source):
    rows=verified_rows(source)
    out=[]
    for r in rows:
        months=r.get("original_photo_observation_months_near30")
        ann=r.get("original_photo_explicit_flower_visibility_near30")
        if not isinstance(months,dict) or not isinstance(ann,dict):
            raise ValueError("Original author source real photo month/flower metadata not available")
        near=r["n_rail_polyline_point_near_30m"]
        if any(not isinstance(n,int) or n<0 for n in months.values()):
            raise ValueError("Actual photo original calendar-month event counts invalid")
        if any(not isinstance(n,int) or n<0 for n in ann.values()):
            raise ValueError("Actual original flower-annotation tag counts invalid")
        if any(not (isinstance(k,str) and len(k)==2 and k.isdigit() and 1<=int(k)<=12) for k in months):
            raise ValueError("Original observation calendar month invalid")
        if sum(months.values())>near:
            raise ValueError("Observation month photos exceed original geolocated near30")
        original_tagged=int(ann.get("flower_visible_annotation",0))
        original_untagged=int(ann.get("no_verified_flower_visible_annotation",0))
        if set(ann)-{"flower_visible_annotation","no_verified_flower_visible_annotation"}:
            raise ValueError("Unexpected source annotation type; do not interpret a negative phenotype")
        if original_tagged+original_untagged!=near:
            raise ValueError("Some original near30 photos missing annotation QA categorization")
        flowerless=(r["taxon"] in FLOWERLESS_TAXA)
        n_months=sum(n>=MIN_MONTHLY_N for n in months.values())
        basic=bool(r["minimum_original_photo_feasibility_gate_passed"])
        image_review_feasible=bool(
            basic and not flowerless
            and r["n_rail_polyline_30m_with_CC_photo_license"]>=MIN_CC
            and r["n_rail_polyline_guaranteed_circle_within_30m"]>=MIN_GUARANTEED30
            and n_months>=MIN_MONTHS
        )
        annotated_training=bool(image_review_feasible and original_tagged>=MIN_ANNOTATED_FLOWER)
        out.append({
            "region":r["region_source_bbox"],
            "species":r["taxon"],
            "original_photo_rail_30m_count":near,
            "original_CC_photo_near30":r["n_rail_polyline_30m_with_CC_photo_license"],
            "original_public_accuracy_circle_guaranteed_near30":r["n_rail_polyline_guaranteed_circle_within_30m"],
            "original_photo_controls_100to1000m":r["n_source_species_photo_100to1000m_controls"],
            "source_near30_unique_observers":r["n_unique_public_observers_near30m"],
            "source_near30_distinct_100m_cells":r["n_distinct_public_100m_cells_near30m"],
            "n_distinct_observed_months_with_at_least5_photos":n_months,
            "source_near30_observation_month_counts_NOT_first_flowering":dict(sorted(months.items())),
            "near30_original_iNat_positive_FLOWERS_annotation":original_tagged,
            "near30_without_VERIFIED_positive_annotation_NOT_a_negative":original_untagged,
            "optional_positive_Flowers_annotation_fraction_of_all_near30_NOT_flowering_prevalence":
                (original_tagged/near if near else None),
            "is_spore_bearing_NON_flowering_horsetail":flowerless,
            "v166_basic_original_photo_gate_passed":basic,
            "v167_ANGIOSPERM_positive_image_REVIEW_logistics_passed":image_review_feasible,
            "v167_annotated_photo_training_minimum_met":annotated_training,
            "not_actual_photo_review_or_plant_height_measured":True
        })
    if len(out)!=40:
        raise ValueError("Original selected 20 species in two source regions not retained")
    eligible=sorted((x for x in out if x["v167_ANGIOSPERM_positive_image_REVIEW_logistics_passed"]),
      key=lambda x:(-x["original_CC_photo_near30"],-x["near30_original_iNat_positive_FLOWERS_annotation"],x["species"]))
    annotate=sorted((x for x in out if x["v167_annotated_photo_training_minimum_met"]),
      key=lambda x:(-x["near30_original_iNat_positive_FLOWERS_annotation"],-x["original_CC_photo_near30"],x["species"]))
    return {
        "complete":True,
        "version":"v16.7_SOURCE_ORIGINAL_20_RAIL_TAXA_2BBOX_FLOWER_ANNOTATION_AND_IMAGE_REVIEW_GATE",
        "source_verified_v166_original_independent_original_DOIs":list(ORIGINAL_BOTANICAL),
        "source_v166_author_photo_gis_ID":"43693580409",
        "source_n_real_eligible_photo_gates_v166":sum(x["v166_basic_original_photo_gate_passed"] for x in out),
        "source_total_gis_near30_photo_records_original":sum(x["original_photo_rail_30m_count"] for x in out),
        "source_total_near30_tagged_FLOWERS_includes_inapplicable_horsetail":sum(x["near30_original_iNat_positive_FLOWERS_annotation"] for x in out),
        "source_total_near30_tagged_FLOWERS_for_angiosperms_only":sum(
          x["near30_original_iNat_positive_FLOWERS_annotation"] for x in out
          if not x["is_spore_bearing_NON_flowering_horsetail"]),
        "n_species_country_cells_passing_exploratory_ANGIOSPERM_image_review_gate":len(eligible),
        "n_species_country_cells_with_min_5_positive_Flowers_tags":len(annotate),
        "exploratory_ANGIOSPERM_image_review_priority":[
          {k:x[k] for k in (
            "region","species","original_photo_rail_30m_count","original_CC_photo_near30",
            "original_public_accuracy_circle_guaranteed_near30",
            "source_near30_unique_observers","source_near30_distinct_100m_cells",
            "n_distinct_observed_months_with_at_least5_photos",
            "near30_original_iNat_positive_FLOWERS_annotation",
            "original_photo_controls_100to1000m"
          )} for x in eligible],
        "original_species_country_40cells_annotation_and_month_audit":out,
        "exploratory_not_predeclared_before_source_photo_count":True,
        "strict_limits":[
            "Source only two country bounding boxes, not global near30m railway counts.",
            "iNaturalist observer non-annotation means unknown photo phenophase, NOT verified nonflowering.",
            "Explicit Flowers annotation can be incorrect/misidentification; original photographic pixels not manually inspected.",
            "Camera photo date and flower observation date do NOT estimate individual first flowering or flowering shift.",
            "Strict public GPS error circles improve geometric confidence, but not right-of-way/trackbed assignment nor railway management exposure.",
            "No image height scale, foreground ROI or same-genotype repeated rail/nonrail phenotyping; no claim of morphological shift or genetic adaptation.",
            "Equisetum arvense cannot flower, so no flowering phenology group despite source annotation metadata.",
            "The stage-2 image review thresholds were introduced after already seeing aggregate source v16.6 photo abundance; they are screening, NOT a confirmatory inference.",
        ],
        "none_of_new_source_counts_identifies_trait_selection":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-v166",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--csv",required=True)
    a=p.parse_args()
    report=evaluate(json.loads(Path(a.source_v166).read_text(encoding="utf-8")))
    o=Path(a.output);o.parent.mkdir(parents=True,exist_ok=True)
    o.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    csv_path=Path(a.csv);csv_path.parent.mkdir(parents=True,exist_ok=True)
    cols=("region","species","original_photo_rail_30m_count","original_CC_photo_near30",
          "original_public_accuracy_circle_guaranteed_near30","original_photo_controls_100to1000m",
          "source_near30_unique_observers","source_near30_distinct_100m_cells",
          "n_distinct_observed_months_with_at_least5_photos",
          "near30_original_iNat_positive_FLOWERS_annotation",
          "near30_without_VERIFIED_positive_annotation_NOT_a_negative",
          "is_spore_bearing_NON_flowering_horsetail",
          "v166_basic_original_photo_gate_passed",
          "v167_ANGIOSPERM_positive_image_REVIEW_logistics_passed",
          "v167_annotated_photo_training_minimum_met")
    with csv_path.open("w",encoding="utf-8",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=cols)
        writer.writeheader()
        for row in report["original_species_country_40cells_annotation_and_month_audit"]:
            writer.writerow({k:row[k] for k in cols})
    print(json.dumps({
      "version":report["version"],
      "n_gis_near30_source":report["source_total_gis_near30_photo_records_original"],
      "v166_original_predeclared_photo_basic_pass":report["source_n_real_eligible_photo_gates_v166"],
      "n_angiosperm_positive_image_review_opportunity":report["n_species_country_cells_passing_exploratory_ANGIOSPERM_image_review_gate"],
      "n_cells_with_at_least5_positive_Flowers_annotations":report["n_species_country_cells_with_min_5_positive_Flowers_tags"],
      "n_photos_with_explicit_Flowers_annotation_angiosperms":report["source_total_near30_tagged_FLOWERS_for_angiosperms_only"],
      "species_eligible_for_manual_image_review":report["exploratory_ANGIOSPERM_image_review_priority"],
      "source_not_phenology_height_fitness_or_genetic_adaptation":True
    },ensure_ascii=False))

if __name__=="__main__":main()
