#!/usr/bin/env python3
"""v16.6 independently audit SUCCESS v16.4 publisher-authenticated top20 taxa and
SUCCESS v16.5 real iNaturalist research-grade photo observations <=30m from
true original 2023/2024 OSM railway polylines, keeping country study BBOX
and GLOBAL photo context strictly separate.

No API/network calls; source GitHub Actions artifacts are pinned to exact
verified workflow run IDs. All 20 x 2 original-region queries must be fully
paginated, have source observations, OSM SHA, and no original exact taxon
selection drift. Report true near30 record count, public coordinate
uncertainty guaranteed30 subset, CC photo license, distinct observer+100m
place cells, 100-1000m comparison record opportunity, AND separate GLOBAL
iNaturalist 2018-25 photo TOTAL (NOT GLOBAL NEAR30 railway count).

Decision threshold fixed in v16.4/5 JAPANESE protocol before iNat result:
- 30 original photo observations satisfying <=30m coordinate accuracy and
  point-to-rail polyline <=30m, >=10 unique observers, >=10 100m spatial
  cells, >=30 same-species 100-1000m source-region photo controls.
This gate is necessary only for detailed phenological/morphology PHOTO
opportunity; it does NOT prove original photographed entire plant is
measurable. Absolute height cannot be reliably inferred from one
unscaled image, and flower visibility != flowering onset.
Equisetum arvense is a spore-bearing plant (NOT a flower-phase candidate).

Underrepresented/unresolved original taxonomic synonyms, obscuration,
uncertain coordinates, spatial mismatch, and area bias never equated with
biological absence or no actual railway adaptation.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path

SOURCE_CANDIDATE_VERSION="v16.4_AUTHOR_ORIGINAL_REAL_RAIL_FLORA_20_SPECIES_FIRST_PHOTO_FEASIBILITY"
PHOTO_SOURCE_VERSION="v16.5_iNaturalist_first_TWO_SOURCE_COUNTRY_BBOX_real_OSM_rail_30m_feasibility"
GEO_SHA={
 "SI":"b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39",
 "CZ":"26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"
}
COUNTRIES=("SI","CZ")
MIN_NEAR30=30
MIN_OBSERVERS=10
MIN_100M_CELLS=10
MIN_CONTROL_100TO1000=30
N_ORIGINAL_SPECIES=20
ORIGINAL_RESEARCH_DATE=("2018-01-01","2025-12-31")

def validate_original_source(candidate,photo):
    if candidate.get("version")!=SOURCE_CANDIDATE_VERSION or candidate.get("complete") is not True:
        raise ValueError("Original 20 railway botanist taxa must be source-authenticated")
    if (photo.get("version")!=PHOTO_SOURCE_VERSION or photo.get("complete") is not True or
        photo.get("NOT_exhaustive_global_rail_30m_count") is not True):
        raise ValueError("Original iNaturalist bbox 30m count source/version incorrectly represented as global")
    if tuple(photo.get("source_2018_2025_dates") or ())!=ORIGINAL_RESEARCH_DATE:
        raise ValueError("Original source iNaturalist 2018-25 photo dates changed")
    if photo.get("source_candidate_20_original_taxa_selection_manifest")!=SOURCE_CANDIDATE_VERSION:
        raise ValueError("Original source photo candidate provenance mismatch")
    selection=candidate.get("original_selected_20_species") or []
    names=[x.get("scientific_name_original_exact_binomial") for x in selection]
    if len(names)!=N_ORIGINAL_SPECIES or len(set(names))!=N_ORIGINAL_SPECIES:
        raise ValueError("Not all 20 original real botanist species from preselected source study")
    geo=photo["country_bbox_rail_geometry_original_source"]
    sources=photo["two_country_bbox_20taxon_real_iNat_photo_source_30m"]
    if set(sources)!=set(COUNTRIES) or set(geo)!=set(COUNTRIES):
        raise ValueError("Slovenia/Czech original photo bbox source study countries missing or invented")
    if not all(geo[c]["original_Geofabrik_PBF_SHA256"]==GEO_SHA[c] for c in COUNTRIES):
        raise ValueError("2023/2024 source authentic historic OSM rail geometry SHA changed")
    if not all(len(sources[c])==20 and
               [r["taxon"] for r in sources[c]]==names for c in COUNTRIES):
        raise ValueError("Original source species selected before photo count order drifted")
    if photo["n_country_species_strata_completed_without_API_cap"]!=40:
        raise ValueError("40 original species-region full iNaturalist source census not completely paginated")
    if photo["n_country_species_strata_api_capped_and_NEAR30_UNESTIMATED_NOT_ZERO"]!=0:
        raise ValueError("One original source species iNaturalist API stratum was capped; never count as 0")
    if not all(r["source_query_complete"] for c in COUNTRIES for r in sources[c]):
        raise ValueError("Original species-country photo observation returns include unresolved cap")
    return names

def run(candidate,photo):
    names=validate_original_source(candidate,photo)
    sources=photo["two_country_bbox_20taxon_real_iNat_photo_source_30m"]
    global_totals={
        r["taxon"]:r["global_iNaturalist_research_photo_total_2018_2025"]
        for r in photo["global_species_research_photo_totals_2018_2025_NOT_global_railway"]}
    if len(global_totals)!=20 or set(global_totals)!=set(names):
        raise ValueError("Missing true source global photo totals for some original species, NOT global near rail")
    results=[]
    for c in COUNTRIES:
        for r in sources[c]:
            m=r["source_original_photo_quality_and_geo_rail_polyline_30m"]
            totals={
              "taxon":r["taxon"],"region_source_bbox":c,
              "source_original_iNat_global_research_grade_photo_total_NOT_global_rail":global_totals[r["taxon"]],
              "source_bbox_iNat_research_photo_total":r["source_API_bbox_photo_total"],
              "n_public_source_bbox_research_photo_pass_30m_coordinate_accuracy":m["n_after_strict_accuracy_and_taxon_public_geoprivacy_gate"],
              "n_rail_polyline_point_near_30m":m["n_point_near_30m_strict_accuracy"],
              "n_rail_polyline_guaranteed_circle_within_30m":m["n_location_certain_within_30m"],
              "n_rail_polyline_30m_with_CC_photo_license":m["n_point_near_30m_photo_CC_licensed"],
              "n_unique_public_observers_near30m":m["n_distinct_observers_public_within_30m"],
              "n_distinct_public_100m_cells_near30m":m["n_distinct_100m_public_cells_near_30m"],
              "n_source_species_photo_100to1000m_controls":m["n_point_100to1000m_strict_accuracy"],
              "n_source_species_photo_30to100m":m["n_point_30to100m_strict_accuracy"],
              "original_photo_observation_months_near30":m["near_30m_month_of_observation_counts_NOT_flowering_onset"],
              "original_photo_explicit_flower_visibility_near30":m["near_30m_explicit_flower_visibility_annotation_counts_NOT_full_phenophase"],
              "is_horsetail_spore_not_flowering_species":r["taxon"]=="Equisetum arvense",
              "photograph_whole_plant_absolute_height_measurable_not_yet_audited":True,
            }
            total=totals["n_public_source_bbox_research_photo_pass_30m_coordinate_accuracy"]
            near=totals["n_rail_polyline_point_near_30m"]
            certain=totals["n_rail_polyline_guaranteed_circle_within_30m"]
            licensed=totals["n_rail_polyline_30m_with_CC_photo_license"]
            if not all(isinstance(x,int) and x>=0 for x in (total,near,certain,licensed)):
                raise ValueError("Original source photo aggregate count missing/non-integer")
            if not(licensed<=near<=total<=r["source_API_bbox_photo_total"]):
                raise ValueError("Original published photo source count > iNaturalist API source denominator")
            if not certain<=near:
                raise ValueError("Original 30m public accuracy-circle guaranteed count cannot exceed coordinate count")
            # Check that strict georeferenced categories exactly conserve
            # the source accepted records (bboxes may include off-country).
            a=near+totals["n_source_species_photo_30to100m"]+totals["n_source_species_photo_100to1000m_controls"]+m["n_point_beyond_1000m_strict_accuracy"]
            if a!=total:
                raise ValueError("Source original strict positional-accuracy photo counts do not partition by real track distance")
            totals["minimum_original_photo_feasibility_gate_passed"]=bool(
              near>=MIN_NEAR30 and
              totals["n_unique_public_observers_near30m"]>=MIN_OBSERVERS and
              totals["n_distinct_public_100m_cells_near30m"]>=MIN_100M_CELLS and
              totals["n_source_species_photo_100to1000m_controls"]>=MIN_CONTROL_100TO1000
            )
            results.append(totals)
    if len(results)!=40:
        raise ValueError("Original two source study bbox × 20 plant taxa source count drifted")
    summed=sum(x["n_rail_polyline_point_near_30m"] for x in results)
    cc=sum(x["n_rail_polyline_30m_with_CC_photo_license"] for x in results)
    if (summed!=photo["n_complete_strata_with_geometrically_point_within30m"] or
        cc!=photo["n_complete_strata_with_photo_license_near30m"]):
        raise ValueError("Source photo study global 30m count accidentally mixed with original SI/CZ bbox aggregates")
    eligible=[x for x in results if x["minimum_original_photo_feasibility_gate_passed"]]
    by_species=[]
    for name in names:
        p=[x for x in results if x["taxon"]==name]
        by_species.append({
          "taxon":name,
          "source_full_global_iNat_research_photo_total_NOT_global_near30":
              global_totals[name],
          "n_30m_railway_photo_across_two_nonoverlapping_study_BBOX":sum(x["n_rail_polyline_point_near_30m"] for x in p),
          "n_30m_geoprivacy_CC_photo_across_two_BBOX":sum(x["n_rail_polyline_30m_with_CC_photo_license"] for x in p),
          "n_public_accuracy_circle_certain_within_30m_across_two_BBOX":sum(x["n_rail_polyline_guaranteed_circle_within_30m"] for x in p),
          "n_country_study_bbox_eligible_for_next_stage":sum(x["minimum_original_photo_feasibility_gate_passed"] for x in p),
          "independent_observer_count_cannot_sum_to_exact_world_unique":True,
        })
    return {
      "complete":True,"version":"v16.6_AUTHENTIC_2COUNTRY_20_ORIGINAL_RAIL_TAXA_PUBLIC_PHOTO_30m_FEASIBILITY_QA",
      "source_2026_11_15_original_rail_botanist_top20_preselected":SOURCE_CANDIDATE_VERSION,
      "source_2018_to_2025_real_iNaturalist_vs_2023_2024_railway_PBF":PHOTO_SOURCE_VERSION,
      "source_original_publications_and_Geofabrik_SHA":photo["country_bbox_rail_geometry_original_source"],
      "source_species_country_counts":results,
      "source_exact_original_20_species_rank_and_two_region_totals":by_species,
      "n_original_taxa":20,
      "n_original_source_country_bbox_taxon_complete":40,
      "n_source_photo_rail_point_within30m_both_bboxes":summed,
      "n_photo_CC_licensed_within30m_both_bboxes":cc,
      "n_location_within_30m_even_with_coordinate_uncertainty_circle":sum(
          x["n_rail_polyline_guaranteed_circle_within_30m"] for x in results),
      "n_country_taxon_cells_passing_predeclared_minimum_photographic_opportunity":len(eligible),
      "original_country_taxon_cells_passing_minimum_photo_gate":[
          {"region_bbox":x["region_source_bbox"],"taxon":x["taxon"],
            "n_geometric_photo_near30m":x["n_rail_polyline_point_near_30m"],
            "n_public_observers":x["n_unique_public_observers_near30m"],
            "n_100m_public_cells":x["n_distinct_public_100m_cells_near30m"],
            "n_source_100to1000m_photo_controls":x["n_source_species_photo_100to1000m_controls"]
          } for x in eligible],
      "predeclared_required_per_country_taxon_photo":{
          "point_near30m_strict_public_accuracy":MIN_NEAR30,
          "distinct_public_observers":MIN_OBSERVERS,
          "distinct_100m_cells":MIN_100M_CELLS,
          "control_photo_100to1000m":MIN_CONTROL_100TO1000
      },
      "not_an_accurate_GLOBAL_near_railway_photo_count":True,
      "not_an_actual_measurement_of_photo_plant_height_or_first_flowering":True,
      "not_adaptation_or_local_railway_selection_evidence":True,
      "measurement_design_limitations":[
          "Two authentic European country source-study BOUNDING BOXES, NOT an exhaustive global railway-proximity analysis.",
          "Global context total is original research-grade photo source count and DOES NOT contain global distance to a railway.",
          "Strict public positional accuracy <=30m + point to original OSM rail LINE <=30m may still fall outside right-of-way. Guaranteed circle within30m is separately counted.",
          "Observer IDs were destroyed in source summary; observer totals across countries cannot be combined as an exact unique observer count.",
          "A CC-licensed photo can be macro of a leaf/flower with no full plant or scale; actual absolute height/photo ROI not measured.",
          "Horsetail Equisetum is spore-bearing; no angiosperm flowers or first-flowering phenology can be inferred.",
          "Date of documented open flower appearance does NOT equal actual first flowering day of an individual.",
          "Photographic sampling effort, season, observer and town proximity are confounded; road and comparable distant rail controls remain to be matched.",
          "A named plant species common in original rail plot is not an independently estimated rail-specialist taxon, and its modern photos do not prove train-carried dispersal."
      ]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-top20-v164",required=True)
    p.add_argument("--source-photo-v165",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--summary-csv",required=True)
    a=p.parse_args()
    candidate=json.loads(Path(a.source_top20_v164).read_text(encoding="utf-8"))
    photo=json.loads(Path(a.source_photo_v165).read_text(encoding="utf-8"))
    x=run(candidate,photo)
    p1=Path(a.output);p1.parent.mkdir(parents=True,exist_ok=True)
    p1.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding="utf-8")
    p2=Path(a.summary_csv);p2.parent.mkdir(parents=True,exist_ok=True)
    fields=["taxon","region_source_bbox","source_original_iNat_global_research_grade_photo_total_NOT_global_rail",
     "source_bbox_iNat_research_photo_total","n_public_source_bbox_research_photo_pass_30m_coordinate_accuracy",
     "n_rail_polyline_point_near_30m","n_rail_polyline_guaranteed_circle_within_30m",
     "n_rail_polyline_30m_with_CC_photo_license","n_unique_public_observers_near30m",
     "n_distinct_public_100m_cells_near30m","n_source_species_photo_100to1000m_controls",
     "minimum_original_photo_feasibility_gate_passed"]
    with p2.open("w",encoding="utf-8",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields)
        writer.writeheader()
        for rec in x["source_species_country_counts"]:
            writer.writerow({k:rec[k] for k in fields})
    print(json.dumps({
        "version":x["version"],
        "n_original_source_species":x["n_original_taxa"],
        "n_country_bbox_full_census_strata":x["n_original_source_country_bbox_taxon_complete"],
        "n_photo_near30m":x["n_source_photo_rail_point_within30m_both_bboxes"],
        "n_photo_CClicensed_near30m":x["n_photo_CC_licensed_within30m_both_bboxes"],
        "n_location_circle_certain_near30m":x["n_location_within_30m_even_with_coordinate_uncertainty_circle"],
        "n_country_species_pass_predeclared_opportunity":x["n_country_taxon_cells_passing_predeclared_minimum_photographic_opportunity"],
        "opportunity_species_country":x["original_country_taxon_cells_passing_minimum_photo_gate"],
        "all_original_20_species_two_region_totals":x["source_exact_original_20_species_rank_and_two_region_totals"],
        "all_country_taxon_quality_counts":[{
          "country":r["region_source_bbox"],"species":r["taxon"],
          "near30":r["n_rail_polyline_point_near_30m"],
          "certain30":r["n_rail_polyline_guaranteed_circle_within_30m"],
          "CC":r["n_rail_polyline_30m_with_CC_photo_license"],
          "observers":r["n_unique_public_observers_near30m"],
          "100m_cells":r["n_distinct_public_100m_cells_near30m"],
          "control100to1000":r["n_source_species_photo_100to1000m_controls"],
          "pass":r["minimum_original_photo_feasibility_gate_passed"]
        } for r in x["source_species_country_counts"]],
        "NOT_global_near30m":True,
    },ensure_ascii=False))

if __name__=="__main__":main()
